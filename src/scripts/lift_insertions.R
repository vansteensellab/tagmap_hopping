#!/usr/bin/env Rscript

library(data.table)
library(argparse)
library(stringr)

parser <- ArgumentParser(description='Process counts from two sides of tagmentation and score putative intergrations')

parser$add_argument('--insert', help='file with putative intergrations or bam file')
parser$add_argument('--chain', help='chain file')
parser$add_argument('--insert-out', help='output for insert file on reference genome')
parser$add_argument('--insert-left', help='left over file for insert file on reference genome')
parser$add_argument('--chrom-sizes', nargs='?', help='file with sizes of chromosomes')
parser$add_argument('--home', nargs='?', help='file with home insertions')
parser$add_argument('--dist', nargs='?', help='distance to home insertions')



argv <- parser$parse_args()

save(argv, file='test.Rdata')


if (grepl('.bam', argv$insert)){
    insert = fread(cmd=paste('samtools view', argv$insert), header=F, sep='\t',
                   fill=T)
    colnames(insert)[1:5] = c('name', 'flag', 'seqnames', 'start', 'mapq')
    coordinates = c('start')

    header = fread(cmd=paste('samtools view -H', argv$insert), header=F, sep='\t',
                             fill=T)
    header[,seqnames:=gsub('SN:','',V2)]

    chromsize = fread(argv$chrom_sizes, col.names=c('seqnames', 'size'),
                      key='seqnames')
    header[,size:=chromsize[seqnames,size]]

    header[!is.na(size), V3:= paste0('LN:', size)]

    this_cmd = paste(commandArgs())
    extra_line = list('@PG', 'ID:tagmap_hopping', 'PN:tagmap_hopping',
                      paste0('CL:', this_cmd), 'VN:0.1')
    new_header = rbind(header[,-c('seqnames', 'size')], extra_line)
} else {
    insert = fread(argv$insert, stringsAsFactors=F)
    insert_order = insert$name
    coordinates = c('start', 'end', 'start_gap', 'end_gap')
    dist = fread(argv$dist, stringsAsFactors=F, key='name',
                 col.names=c('seq', 'start', 'end', 'name',
                             'p_adj', 'strand', 'seq_home', 'start_home',
                             'end_home', 'name_home', 'score_home', 'strand_home',
                             'dist'))
    dist = dist[,.SD[1,], by='name']
    home = fread(argv$home, stringsAsFactors=F, key='name')
    insert[,name_home:=dist[name, name_home]]
    insert[,c('start_home', 'end_home'):=home[name_home,c('start_gap','end_gap')]]
}

setkey(insert, 'name')

lift_list = lapply(coordinates, function(x){
    bed_dt = insert[,c('seqnames', x, 'name'), with=F]
    colnames(bed_dt) = c('seqnames', 'end', 'name')
    bed_dt[,end:=as.numeric(end)]
    bed_dt[,start:=end-1]
    bed_dt = bed_dt[!is.na(end), c('seqnames', 'start', 'end', 'name'), with=F]
    out = system2(command=c('liftOver', '/dev/stdin', argv$chain, '/dev/stdout',
                            argv$insert_left),
                  input=do.call(paste, c(bed_dt, sep='\t')),
                  stdout=T)
    if (length(out) == 0){
        out_dt = data.table(NA, bed_dt$name)
        colnames(out_dt) = c(x, 'name')
        setkey(out_dt, 'name')
    } else {
        out_dt = data.table(do.call(rbind,str_split(out, '\t')))[,c(3,4)]
        colnames(out_dt) = c(x, 'name')
        setkey(out_dt, 'name')
    }
    return(out_dt)
})

join <- function(dt1, dt2){
    return(merge(dt1, dt2, all=T))
}

lift_dt = Reduce(join, lift_list)

out = merge(lift_dt, insert[,-coordinates,with=F])[, colnames(insert), with=F]

if (grepl('.bam', argv$insert)){
    sam_dt = rbind(new_header, out[order(seqnames, as.numeric(start))], fill=T)
    sam = do.call(paste, c(sam_dt, sep='\t'))
    sam_trim = str_trim(gsub('NA','',sam))
    system2(command=c('samtools', 'view', '-Sb', '/dev/stdin'),
            input=sam_trim, stdout=argv$insert_out)

    system2(command=c('samtools', 'index', argv$insert_out))
} else {
    out[is.na(start), start:=start_home]
    out[is.na(end), end:=end_home]
    fwrite(out[insert_order, -c('name_home', 'start_home', 'end_home')],
           file=argv$insert_out, sep='\t', quote=F, na='.')
}

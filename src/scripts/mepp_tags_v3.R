#!/usr/bin/env Rscript

library(rtracklayer)
library(data.table)
library(argparse)
library(stringr)
library(metap)

parser <- ArgumentParser(description='Process counts from two sides of tagmentation and score putative intergrations')

parser$add_argument('--bed', help='bed file with putative intergrations')
parser$add_argument('--bam', nargs="*", help='sorted bam files with mapping for both sides')
parser$add_argument('--overhang', help='expected overhang of the intergration (e.g. TTAA, TA)')
parser$add_argument('--fasta', help='fasta file of reference genome (has to be indexed first)')
parser$add_argument('--out', help='file for result output')
parser$add_argument('--mate', '-m',
                    help=paste('which mate of the mate pair comes from the',
                               'intergration (not Tn5).',
                               'Only used for single sided tagmap'))

argv <- parser$parse_args()

save(argv, file='test.Rdata')

call_gap <- function(dt){
    start_count = dt[,sort(table(start_gap),decreasing=TRUE)[1]]
    gap_concordance = start_count / nrow(dt)
    end_count = dt[,sort(table(end_gap),decreasing=TRUE)[1]]
    return(list(start_gap = as.numeric(names(start_count)),
                end_gap = as.numeric(names(end_count)),
                gap_concordance = gap_concordance))

}

concordance_side <- function(dt, this_gap){
    start_g = this_gap[,start_gap]
    end_g = this_gap[,end_gap]
    side = dt[,ifelse(start_g-start_pair > end_pair-end_g, 1, 2)]
    concordance = dt[,(side==1 & strand2=='-') | (side==2 & strand2=='+') ]
    return(list(side=side, is_concordant=concordance, mapq=dt$mapq,
                strand=dt$strand2, fragment_size = dt$fragment_size))
}

find_center <- function(dt, dist){
    delta_dt = data.table(start=unlist(dt[1:(nrow(dt)-(dist+1)), 'start']),
                          left=unlist(dt[1:(nrow(dt)-(dist+1)), 'difference']),
                          right=unlist(dt[(dist+2):nrow(dt), 'difference']))
    delta_dt[,delta:=left-right]
    gap_start = delta_dt[(left < 0 & right > 0) | (right < 0 & left > 0),
                         start[which.max(abs(delta))]]
    return(list(start_gap=gap_start,
                end_gap=gap_start + dist))
} ##TODO: find center with zero length overhang

find_single_center <- function(dt, start, end, dist){
    strand_table = dt[, table(strand)]
    strand_i = which.max(strand_table)
    concordance = strand_table[strand_i] / sum(strand_table)
    strand = names(strand_table)[strand_i]
    start_gap = ifelse(strand=="+", as.numeric(start) - 1,
                       as.numeric(end) - dist)
    return(list(strand=strand,
                concordance=concordance,
                start_gap = start_gap,
                end_gap = start_gap + ifelse(dist==0,1,dist)))
}

if (!length(argv$bam) %in% c(1,2)){
    stop('only 1 or 2 bam files allowed')
}


write_empty <- function(argv){
    if (length(argv$bam) == 2){
        cat(paste("seqnames", "start", "end", "start_gap", "end_gap",
                  "center_sequence", "count_1", "D_1", "mapq_1", "p_1",
                  "concordance.x", "count_2", "D_2", "mapq_2", "p_2",
                  "concordance.y", "sump", "p_adj", "strand", "name", sep="\t"),
            file=argv$out)
    } else {
        cat(paste("seqnames", "start", "end", "start_gap", "end_gap",
                  "center_sequence", "count_1", "D_1", "mapq_1", "p_1",
                  "concordance", "p_adj", "strand", "name", sep="\t"),
            file=argv$out)
    }
}

bed_gr = reduce(import.bed(argv$bed))
bed_dt = as.data.table(bed_gr)
bed_dt[,region:=1:nrow(bed_dt)]

if (length(bed_gr)==0){
    write_empty(argv)
} else {

    overhang_len = ifelse(argv$overhang%in%c('-', '0'), 0,
                          str_length(argv$overhang))


    n_regions = length(bed_gr)
    region_rep = rep(1:n_regions, each=2)
    side_rep = rep(c(1,2), n_regions)
    region_dt = data.table(region=region_rep,
                           side=side_rep)
    region_dt = region_dt[, list(strand=c('+','-')),
                        by=c('region','side')]
    region_dt[,is_concordant:=(side==1 & strand=='-') |
                              (side==2 & strand=='+')]




    if (length(argv$bam) == 2){
        pair_list = lapply(1:2, function(i){
            bam = argv$bam[i]
            bam_dt = fread(cmd=paste('samtools sort -n', bam, '|',
                                     'bamToBed -bedpe -mate1 -i -'),
                           col.names=c('seqnames1', 'start1', 'end1',
                                       'seqnames2', 'start2', 'end2',
                                       'pair_id', 'mapq', 'strand1', 'strand2'))
            bam_dt[,start_gap:=ifelse(strand2=='+', start2, end2-overhang_len)]
            bam_dt[,end_gap:=start_gap + overhang_len]

            bam_dt[,start_pair:=ifelse(strand2=='+', start2, start1)]
            bam_dt[,end_pair:=ifelse(strand2=='+', end1, end2)]

            bam_dt[,fragment_size := end_pair - start_pair]

            bam_gr=bam_dt[, GRanges(seqnames=seqnames1,
                                    IRanges(start_gap, width=overhang_len))]

            o = findOverlaps(bam_gr, bed_gr)
            bam_merge = cbind(bam_dt[queryHits(o), ], bed_dt[subjectHits(o), ])
            return(bam_merge)
        })

        gap_dt = do.call(rbind, pair_list)
        gap_call = gap_dt[,call_gap(.SD), by='region']
        colnames(gap_call)[1] = 'reg'
        setkey(gap_call, 'reg')

        overlap_list = lapply(1:2, function(i){
            pair_count = pair_list[[i]]
            if(is.null(pair_count)){
                overlap_count = data.table(strand=character(),
                                         region=integer(),
                                         side=integer(),
                                         is_concordant=logical(),
                                         count=numeric(),
                                         mapq=numeric(),
                                         max_fragment=numeric())
            } else {
                call_dt = pair_count[,concordance_side(.SD, gap_call[reg==region,]) , by='region']


                overlap_count = call_dt[,list(count=nrow(.SD), mapq=sum(as.numeric(mapq)),
                                              max_fragment=max(fragment_size)),
                                            by=c('strand', 'region', 'side', 'is_concordant')]
            }
            overlap_merge = merge(region_dt, overlap_count, all=T,
                                  by=c('region', 'side', 'strand', 'is_concordant'))
            overlap_merge[is.na(overlap_merge)] <- 0
            balance = overlap_merge[,list(count=sum(count),
                                          mapq=sum(mapq)/sum(count),
                                          max_fragment=max(max_fragment),
                                          D=(count[side==2]/sum(count)-.5) * 2,
                                          p=ifelse(sum(count)>0,
                                                   stats::binom.test(count[side==1],
                                                                     sum(count))$p.value,
                                                   1)),
                                    by=c('region', 'is_concordant')]
            concordance = balance[,list(count=sum(count),
                                        D=D[is_concordant],
                                        mapq=mapq[is_concordant],
                                        max_fragment=max_fragment[is_concordant],
                                        p=p[is_concordant],
                                        concordance=count[is_concordant] /
                                                    sum(count)),
                                  by=c('region')]

            names(concordance)[2:ncol(concordance)] =
                paste0(names(concordance)[2:ncol(concordance)], '_', i)
            return(concordance)
        })


        overlap_dt = merge(overlap_list[[1]],
                           overlap_list[[2]], by=c('region'),
                           all=T)

        overlap_dt[, sump:=sump(c(p_1, p_2)), by='region']
        overlap_dt[, p_adj:=p.adjust(sump)]
        overlap_dt[, strand:=ifelse(D_1==1, '+', '-')]
        overlap_dt[is.na(D_1), strand:=ifelse(D_2==1, '-', '+')]
        overlap_merge = merge(gap_call, overlap_dt, by.x='reg', by.y='region')
        full_result = merge(bed_dt[,-c('strand', 'width')],
                            overlap_merge, by.x='region', by.y='reg')

    } else if (length(argv$bam) == 1){
        bam = argv$bam[1]
        cmd = paste('bedtools intersect -wa -wb -a', argv$bed, '-b', bam)
        overlap = system(cmd,intern=T)
        overlap_dt = as.data.table(do.call(rbind, str_split(overlap, '\t')))
        names(overlap_dt) = c('seqnames_region', 'start_region', 'end_region', 'mapq_region',
                              'seqnames_read', 'start_read',
                              'end_read', 'read_id', 'mapq', 'strand')
        overlap_dt[,start_region:=as.numeric(start_region)+1]
        overlap_dt[,end_region:=as.numeric(end_region)]
        overlap_dt[,mate:=ifelse(grepl('/1', read_id), 1, 2)]
        overlap_dt[,pair_id:= gsub('/[12]', '', read_id)]
        overlap_dt[,adj_strand:= ifelse(mate==1, strand,
                                        ifelse(strand=='+', '-', '+'))]

        pair_count = overlap_dt[, list(start_pair=min(as.numeric(start_read)),
                                       end_pair=max(as.numeric(end_read)),
                                       mapq=max(mapq), strand_pair=adj_strand[1]),
                                by=c('seqnames_region', 'start_region', 'end_region',
                                     'pair_id')]

        pair_count[,start_gap:=ifelse(strand_pair=='+', end_pair - overhang_len,
                                      start_pair)]
        pair_count[,end_gap:=ifelse(strand_pair=='+', end_pair,
                                    start_pair + overhang_len)]

        pair_merge = merge(bed_dt, pair_count, by.x=c('seqnames', 'start', 'end'),
                           by.y=c('seqnames_region', 'start_region', 'end_region'))

        gap_call = pair_merge[,call_gap(.SD), by='region']
        colnames(gap_call)[1] = 'reg'
        setkey(gap_call, 'reg')


        call_dt = pair_merge[,concordance_side(.SD, gap_call[reg==region,]) , by='region']


        overlap_count = call_dt[,list(count=nrow(.SD), mapq=sum(as.numeric(mapq))),
                                    by=c('strand', 'region', 'side', 'is_concordant')]
        overlap_merge = merge(region_dt, overlap_count, all=T,
                              by=c('region', 'side', 'strand', 'is_concordant'))
        overlap_merge[is.na(overlap_merge)] <- 0
        balance = overlap_merge[,list(count=sum(count),
                                      mapq=sum(mapq)/sum(count),
                                      D=(count[side==1]/sum(count)-.5) * 2,
                                      p=ifelse(sum(count)>0,
                                               stats::binom.test(count[side==1],
                                                                 sum(count))$p.value,
                                               1)),
                                by=c('region', 'is_concordant')]
        concordance = balance[,list(count=sum(count),
                                    D=D[is_concordant],
                                    mapq=mapq[is_concordant],
                                    p=p[is_concordant],
                                    concordance=count[is_concordant] /
                                                sum(count)),
                              by=c('region')]

        concordance[, p_adj:=p.adjust(p)]
        concordance[, strand:=ifelse(D==1, '+', '-')]
        overlap_merge = merge(gap_call, concordance, by.x='reg', by.y='region')
        full_result = merge(bed_dt[,-c('strand', 'width')],
                            overlap_merge, by.x='region', by.y='reg')

    } else {
        stop('only 1 or 2 bam files allowed')
    }

    full_result[, name:=paste0('insert_', .I)]

    lines = full_result[,paste(seqnames, start_gap, end_gap, name,
                               sep='\t', collapse='\n')]

    cmd = paste('bedtools getfasta -name -tab -fi', argv$fasta, '-bed /dev/stdin')

    getfasta = system(cmd, input=lines, intern=T)

    seq_dt = data.table(str_match(getfasta, '(.*)::(.*):(.*)-(.*)\t(.*)'))
    names(seq_dt) = c('string', 'insert', 'seqnames', 'start_gap', 'end_gap',
                      'center_seq')
    setkey(seq_dt, 'insert')
    full_result[,center_sequence:=toupper(seq_dt[name,center_seq])]

    fwrite(full_result[,-c('region')], argv$out, sep='\t', na='.', quote=F)
}

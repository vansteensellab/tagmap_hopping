#!/usr/bin/env Rscript

library(rtracklayer)
library(data.table)
library(argparse)
library(stringr)
library(metap)

parser <- ArgumentParser(description='Process counts from two sides of tagmentation and score putative intergrations')

parser$add_argument('--bw', nargs="*",
                    help=paste('one or two bigWig files from both sides of',
                               'putative integrations counting the mate from',
                               'the integration only'))
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
    concordance = dt[,(side==1 & strand_pair=='+') | (side==2 & strand_pair=='-') ]
    return(list(side=side, is_concordant=concordance, mapq=dt$mapq,
                strand=dt$strand_pair))
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
    if (length(argv$bw) == 2){
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
    region_dt[,is_concordant:=(side==1 & strand=='+') |
                              (side==2 & strand=='-')]




    if (length(argv$bw) == 2){
        pair_list = lapply(1:2, function(i){
            bam = argv$bam[i]
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
            return(pair_merge)
        })

        gap_dt = do.call(rbind, pair_list)
        gap_call = gap_dt[,call_gap(.SD), by='region']
        colnames(gap_call)[1] = 'reg'
        setkey(gap_call, 'reg')

        overlap_list = lapply(1:2, function(i){
            pair_count = pair_list[[i]]
            call_dt = pair_count[,concordance_side(.SD, gap_call[reg==region,]) , by='region']


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
            concordance = balance[,list(count=count[is_concordant],
                                        D=D[is_concordant],
                                        mapq=mapq[is_concordant],
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

        overlap_dt[, sump:=metap::sump(c(p_1, p_2))$p, by='region']
        overlap_dt[, p_adj:=p.adjust(sump)]
        overlap_dt[, strand:=ifelse(D_1==1, '+', '-')]
        overlap_dt[is.na(D_1), strand:=ifelse(D_2==1, '-', '+')]
        overlap_merge = merge(gap_call, overlap_dt, by.x='reg', by.y='region')
        full_result = merge(bed_dt[,-c('strand', 'width')],
                            overlap_merge, by.x='region', by.y='reg')

    } else if (length(argv$bw) == 1){
        r1 = import.bw(argv$bw[1])

        o = findOverlaps(bed_gr, r1)

        score_dt = data.table(r1=r1[subjectHits(o)]$score,
                              region=queryHits(o),
                              seqname=as.character(seqnames(r1[subjectHits(o)])),
                              start=start(r1[subjectHits(o)]))
        region_dt = as.data.table(bed_gr)
        region_dt[,region:=1:nrow(region_dt)]
        lines = region_dt[,paste(seqnames, start, end, region, sep='\t')]
        cmd = paste('bedtools intersect -F 0.5 -wa -wb -a /dev/stdin -b', argv$bam)
        overlap = system(cmd, input = lines,intern=T)

        overlap_dt = data.table(t(as.data.frame(str_split(overlap, '[\t/]'))))
        colnames(overlap_dt) = c('seqnames', 'start', 'end', 'region',
                                 'seqnames_r', 'start_r', 'end_r', 'read_id',
                                 'mate', 'mapq', 'strand')

        center_dt = overlap_dt[mate==argv$mate,
                               find_single_center(.SD, start, end, overhang_len),
                               by=c('seqnames', 'start', 'end', 'region')]
        if (nrow(center_dt) == 0){
            write_empty(argv)
        } else {
            lines = center_dt[,paste(seqnames, start_gap, end_gap, region,
                                sep='\t', collapse='\n')]

            cmd = paste('bedtools getfasta -tab -fi', argv$fasta, '-bed /dev/stdin')

            getfasta = system(cmd, input=lines, intern=T)

            seq_dt = data.table(str_match(getfasta, '(.*):(.*)-(.*)\t(.*)'))
            names(seq_dt) = c('string', 'seqnames', 'start_gap', 'end_gap', 'center_sequence')

            seq_dt[,c('start_gap', 'end_gap'):=list(as.numeric(start_gap),
                                               as.numeric(end_gap))]
            seq_dt[,center_sequence:=toupper(center_sequence)]

            center_merge = merge(center_dt, seq_dt[,-c('string')],
                by=c('seqnames', 'start_gap', 'end_gap'))
            full_result = center_merge[,c('seqnames', 'start', 'end', 'start_gap',
                                          'end_gap', 'center_sequence', 'concordance',
                                          'strand', 'region')]
        }
    } else {
        stop('only 1 or 2 bigwigs allowed')
    }

    if (!(length(argv$bw) == 1 && nrow(center_dt) == 0)){
        full_result[, name:=paste0('insert_', .I)]

        fwrite(full_result[,-c('region')], argv$out, sep='\t', na='.', quote=F)

    }
}

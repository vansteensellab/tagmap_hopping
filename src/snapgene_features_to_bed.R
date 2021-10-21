#!/usr/bin/env Rscript
library(data.table)
library(argparse)

parser <- ArgumentParser(description='Convert snapGene file to in silico bed file')

parser$add_argument('--bed', '-b', help='output file')
parser$add_argument('--features', '-f', help='file with features from snapGene')
parser$add_argument('--chr', '-c', help='chromosome to transfer coordinates to')
parser$add_argument('--start', '-s', type="integer",
                    help='starting location of insertion')
argv <- parser$parse_args()

opt = commandArgs(trailingOnly = TRUE)
# dt = fread(cmd="sed 's/\\.\\./\t/' data/Features_from_Jessica_Zuin_FMI_Targetingvector.txt",
dt = fread(cmd=paste("sed 's/\\.\\./\t/'", argv$features),
           col.names=c('name', 'start_char', 'end_char'),
           colClasses=list(character=1:3))

dt[,start:=as.numeric(gsub('[.]','', start_char)) + argv$start - 1]
dt[,end:=as.numeric(gsub('[.]','', end_char)) + argv$start - 1]

dt[,name:=gsub("'", 'prime', name)]
dt[,name:=gsub(" ", '_', name)]
dt[,name:=gsub(" ", '_', name)]

dt[, chr:=argv$chr]

fwrite(dt[,c('chr', 'start', 'end', 'name')], sep='\t', col.names=F,
       file=argv$bed)

#!/bin/bash

# This script was written to detect (putative) transposon integration sites
# from TagMap paired-end sequencing reads.
# 15 November 2021
# Koen Rademaker, adapted from Christ Leemans
#
# The script exploits TagMap primer design, in which read R2 will contain the
# duplicated transposon overhang sequence (e.g. TTAA for PiggyBac) at the first
# or last base pairs on the forward strand of reverse strand, respectively.
#
# Reads have been filtered earlier in the pipeline to only include proper paired
# reads and non-PCR duplicates. In this script, reads are additionally filtered
# by mapping quality (mapq).
#
# Data of detected transposon integration sites (TIS) are outputted for
# downstream pipeline processing in the following tab-separated structure:
# |read name|chr|read start|read end|strand|mapq|overhang|TIS start|TIS end|sample|
#
# Developed using Sambamba version 0.6.6, Awk version 5.0.1.
# Example usage: find_insertions.sh -b FW.bam RV.bam -m 10 -i TA -o out.txt
#
#------------------------------------------------------------------------------#

usage="$(basename "$0") [-h] -b <...> -m <INT> -s <...) -o <...>

script to find transposon insertion sites.

where:
    -h show this help text
    -b 1 or 2 bam files
    -m minimum for highest mapq on each side
    -s expected transposon insertion site overhang sequence (e.g. TA, TTAA)
    -o output file"


while getopts ':h:b:m:i:o:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    b) echo $OPTARG
       bam=($OPTARG)
       ;;
    m) mapq="$OPTARG"
       ;;
    s) insertion="$OPTARG"
       ;;
    o) out="$OPTARG"
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1
        ;;
  esac
done
shift $((OPTIND - 1))

#------------------------------------------------------------------------------#

# (0) Declare variables and initialise files
BAM_QNAME=1
BAM_CHR=3
BAM_SEQSTART=4
BAM_MAPQ=5
BAM_SEQ=10
INSERTION_LENGTH=${#insertion}

touch "${out}"


# (1) Extract insertion site sequences, whilst accounting for transposon
#     characteristics, e.g. length of the sequence (PiggyBac, TTAA = 4), and
#     accepting alternative sequences as observed in the reads.


# (1a)   Single Tagmentation reaction (e.g. forward, single input .bam file):
# (1a.1) Paired-end read R2 on forward (+) strand.
sambamba view -f sam -F 'mapping_quality >= '${mapq}' and second_of_pair and mate_is_reverse_strand' "${bam[0]}" | \
  awk -vOFS="\t" -vFILE=$(basename -s .bam ${bam[0]}) -vQ="$BAM_QNAME" -vCHR="$BAM_CHR" -vSEQSTART="$BAM_SEQSTART" -vMAPQ="$BAM_MAPQ" -vSEQ="$BAM_SEQ" -vL="$INSERTION_LENGTH" \
  '{print $Q,$CHR,$SEQSTART,$SEQSTART+length($SEQ)-1,"+",$MAPQ,substr($SEQ,1,L),$SEQSTART,$SEQSTART+(L-1),FILE}' >> "${out}"
# (1a.2) Paired-end read R2 on reverse (-) strand.
sambamba view -f sam -F 'mapping_quality >= '${mapq}' and second_of_pair and reverse_strand' "${bam[0]}" | \
  awk -vOFS='\t' -vFILE=$(basename -s .bam ${bam[0]}) -vQ="$BAM_QNAME" -vCHR="$BAM_CHR" -vSEQSTART="$BAM_SEQSTART" -vMAPQ="$BAM_MAPQ" -vSEQ="$BAM_SEQ" -vL="$INSERTION_LENGTH" \
  '{print $Q,$CHR,$SEQSTART,$SEQSTART+length($SEQ)-1,"-",$MAPQ,substr($SEQ,length($SEQ)-(L-1),length($SEQ)),$SEQSTART+length($SEQ)-L,$SEQSTART+length($SEQ)-1,FILE}' >> "${out}"


# (1b)    Two Tagmentation reactions (two input .bam files):
if [ "${bam[1]}" != "" ]
then
  # (1b.1) Paired-end read R2 on forward (+) strand
  sambamba view -f sam -F 'mapping_quality >= '${mapq}' and second_of_pair and mate_is_reverse_strand' "${bam[1]}" | \
    awk -vOFS="\t" -vFILE=$(basename -s .bam ${bam[1]}) -vQ="$BAM_QNAME" -vCHR="$BAM_CHR" -vSEQSTART="$BAM_SEQSTART" -vMAPQ="$BAM_MAPQ" -vSEQ="$BAM_SEQ" -vL="$INSERTION_LENGTH" \
    '{print $Q,$CHR,$SEQSTART,$SEQSTART+length($SEQ)-1,"+",$MAPQ,substr($SEQ,1,L),$SEQSTART,$SEQSTART+(L-1),FILE}' >> "${out}"
  # (1b.2) Paired-end read R2 on reverse (-) strand
  sambamba view -f sam -F 'mapping_quality >= '${mapq}' and second_of_pair and reverse_strand' "${bam[1]}" | \
    awk -vOFS='\t' -vFILE=$(basename -s .bam ${bam[1]}) -vQ="$BAM_QNAME" -vCHR="$BAM_CHR" -vSEQSTART="$BAM_SEQSTART" -vMAPQ="$BAM_MAPQ" -vSEQ="$BAM_SEQ" -vL="$INSERTION_LENGTH" \
    '{print $Q,$CHR,$SEQSTART,$SEQSTART+length($SEQ)-1,"-",$MAPQ,substr($SEQ,length($SEQ)-(L-1),length($SEQ)),$SEQSTART+length($SEQ)-L,$SEQSTART+length($SEQ)-1,FILE}' >> "${out}"
fi

# END OF SCRIPT, THE PIPELINE WILL CONTINUE USING THIS OUTPUT.

#------------------------------------------------------------------------------#

#!/bin/bash

usage="$(basename "$0") [-h] -c <...> -b [...] -c <...> -d <...> -g <INT> -t <...>

script to find overlap between home alignments and alignments of the hopping.

where:
    -h show this help text
    -a homology arms
    -i bowtie2 index
    -o output
"

while getopts ':h:a:i:o:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    a) arms="$OPTARG"
       ;;
    i) index="$OPTARG"
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


bowtie2 -f -U $arms -x $index | \
    samtools view -Sb | \
    bamToBed | \
    awk -vOFS='\t' '
        BEGIN {
            print "seqnames", "start_gap", "end_gap", "strand"
        }
        {
            if (NR==1){
                end=0
                start=$2
                end=$3
            } else {
                seq=$1
                start=$2<start?$2:start
                end=$3>end?$3:end
        }}END{
            print seq, start, end, "+"
        }' > $out

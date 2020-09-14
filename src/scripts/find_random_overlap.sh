#!/bin/bash

usage="$(basename "$0") [-h] -c <...> -b [...] -c <...> -d <...> -g <INT> -t <...>

script to find overlap between home alignments and alignments of the hopping.

where:
    -h  show this help text
    -b bed file for insertions (will be created)
    -r file with home regions (both real insertions as well as random alignments)
    -i insertion table
    -n label name for name in insert file
    -s label name for score in insert file
    -o output with overlap counts
"

score="p_adj"
label="name"

while getopts ':h:r:b:o:n:m:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    m) mast="$OPTARG"
       ;;
    r) regions="$OPTARG"
       ;;
    b) insert_bed="$OPTARG"
       ;;
    o) overlap="$OPTARG"
       ;;
    n) N="$OPTARG"
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

zegrep -v "^#" $mast | \
  shuf -n $N |
  awk -vOFS='\t' '{print $1, $5, $6}' | \
  bedtools sort -i - > $insert_bed

bedtools closest -d -a $insert_bed -b <(bedtools sort -i $regions) > $overlap

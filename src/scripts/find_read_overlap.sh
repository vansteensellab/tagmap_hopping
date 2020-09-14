#!/bin/bash

usage="$(basename "$0") [-h] -c <...> -b [...] -c <...> -d <...> -g <INT> -t <...>

script to find overlap between home alignments and alignments of the hopping.

where:
    -h  show this help text
    -b bam file with alignments of hopping element.
    -r file with home regions (both real insertions as well as random alignments)
    -t temporary bed file
    -c amount of threads/cores
    -o output with overlap counts
"

while getopts ':h:b:r:t:o:c:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    b) bam="$OPTARG"
       ;;
    r) regions="$OPTARG"
       ;;
    t) bed="$OPTARG"
       ;;
    o) overlap="$OPTARG"
       ;;
    c) threads="$OPTARG"
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


samtools sort -@ $threads -n $bam | \
    bamToBed -bedpe -i /dev/stdin | \
    awk -vOFS='\t' '{if ($2 > 0){print $1, ($2<$5)?$2:$5, ($3>$6)?$3:$6, $7, $8, $9}}' | \
    bedtools sort -i - > $bed
echo $bed

bedtools closest -d -a $bed -b <(bedtools sort -i $regions) > $overlap

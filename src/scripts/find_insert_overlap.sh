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

while getopts ':h:r:b:o:n:s:i:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    i) insert="$OPTARG"
       ;;
    r) regions="$OPTARG"
       ;;
    b) insert_bed="$OPTARG"
       ;;
    o) overlap="$OPTARG"
       ;;
    n) label="$OPTARG"
       ;;
    s) score="$OPTARG"
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

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


$DIR/insert_to_bed.awk -vL=$label -vS=$score $insert > $insert_bed


bedtools closest -d -a $insert_bed -b <(bedtools sort -i $regions) > $overlap

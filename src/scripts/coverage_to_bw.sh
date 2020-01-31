#!/bin/bash
usage=("$(basename "$0") [-h] -a <...> -b <...> -c <...> -o <...> -m [0-2]

script to run bedtools coverage to get per base coverage on sites of interest.

where:
    -h  show this help text
    -a  bed file with regions of interest
    -b  bam file
    -c  chrom sizes
    -m  mate pair of interest (0 for both)
    -o  output file
    ")

threads=1

while getopts ':h:a:b:o:c:m:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    a) bed="$OPTARG"
       ;;
    b) bam="$OPTARG"
       ;;
    o) out="$OPTARG"
       ;;
    c) cs="$OPTARG"
       ;;
    m) mate="$OPTARG"
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

tmp=$(tempfile)



if [ $mate == 0 ]
then
    echo 'using all reads'
    sambamba view --format bam $bam -o $tmp.bam
elif [ $mate == 1 ]
then
    echo 'using only first mate of read pair'
    sambamba view --format bam -F "first_of_pair" $bam -o $tmp.bam
elif [ $mate == 2 ]
then
    echo 'using only second mate of read pair'
    sambamba view --format bam -F "second_of_pair" $bam -o $tmp.bam
fi

bedtools coverage -d -a <(awk -vOFS='\t' '{print $1, $2, $3}' $bed) \
                  -b $tmp.bam | \
    awk -vOFS='\t' '{print $1, $2+$4-1, $2+$4, $5}' | \
    sort -k1,1 -k2,2n | uniq > $tmp


bedGraphToBigWig $tmp $cs $out

rm $tmp.bam
rm $tmp

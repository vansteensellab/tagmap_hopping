#!/bin/bash
usage="$(basename "$0") [-h] -c <...> -b [...] -c <...> -d <...> -g <INT> -t <...>

script to find closeby bed regions and find max mapq in region

where:
    -h show this help text
    -c chain file
    -i file with insertsions
    -b 1 or 2 bed files for fwd and reverse reads
    -l left over file
    -j out file for insertions
    -r regions
    -n number of sides used to locate home
"



while getopts ':h:c:i:b:l:j:r:n:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    c) chain="$OPTARG"
       ;;
    i) insert="$OPTARG"
       ;;
    b) bed=("$OPTARG")
       ;;
    l) left="$OPTARG"
       ;;
    j) insert_out="$OPTARG"
       ;;
    r) region="$OPTARG"
       ;;
    n) n="$OPTARG"
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


if [ $n == 2 ]
then
    tail -n+2 $insert |  \
        awk -vOFS='\t' '{print $1, $2, $3, "insert_"NR, $18, $19}' > $insert_out
elif [ $n == 1 ]
then
    tail -n+2 $insert |  \
        awk -vOFS='\t' '{print $1, $2, $3 + 1, "insert_"NR, $7, $8}' > $insert_out
fi


cat ${bed[*]} | \
    bedtools intersect -v -a /dev/stdin -b $insert_out | \
    awk -vOFS='\t' '{print $1, $2, $3, "grey_"NR, ".", "."}' | \
    cat $insert_out /dev/stdin | \
    liftOver /dev/stdin $chain $region $left

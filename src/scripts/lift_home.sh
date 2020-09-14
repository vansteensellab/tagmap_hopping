#!/bin/bash
usage="$(basename "$0") [-h] -c <...> -b [...] -c <...> -d <...> -g <INT> -t <...>

script to find closeby bed regions and find max mapq in region

where:
    -h show this help text
    -c chain file
    -d left over file
    -i file with insertsions
    -b 1 or 2 bed files for fwd and reverse reads
    -l regions with in silico coordinate system
    -j out file for insertions
    -r regions with reference genome coordinate system
    -n name column in insertion file
    -s score column in insertion file
"

score="p_adj"
label="name"

while getopts ':h:c:i:b:d:l:j:r:n:s:' option; do
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
    d) left="$OPTARG"
       ;;
    j) insert_out="$OPTARG"
       ;;
    r) region_ref="$OPTARG"
       ;;
    l) region_lift="$OPTARG"
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


awk -vOFS='\t' -vL=$label -vS=$score '
    BEGIN{
        arr["seqnames"] = 1
        arr["start"] = 2
        arr["end"] = 3
        arr[L] = -1
        arr[S] = -1
        arr["strand"] = -1
    } {
    if (NR==1){
        for (i = 1; i <= NF; i++){
            if ($i in arr){
                arr[$i] = i
            }
        }
    } else {
        name = arr[L]==-1?"insert_"NR-1:$arr[L]
        score = arr[S]==-1?".":$arr[S]
        strand = arr["strand"]==-1?".":$arr["strand"]
        print $arr["seqnames"], $arr["start"], $arr["end"], name, score, strand
    }
}' $insert > $insert_out


if [ -z "$bed" ]
then
    cp $insert_out $region_ref
else
    cat ${bed[*]} | \
        bedtools intersect -v -a /dev/stdin -b $insert_out | \
        awk -vOFS='\t' '{print $1, $2, $3, "grey_"NR, ".", "."}' | \
        cat $insert_out /dev/stdin > $region_ref
fi

liftOver -minMatch=0 $region_ref $chain $region_lift $left

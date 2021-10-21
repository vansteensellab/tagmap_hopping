#!/bin/bash
usage="$(basename "$0") [-h] -a <...> -b <...> -c <...> -d <...> -g <INT> -t <...>

script to find closeby bed regions and find max mapq in region

where:
    -h  show this help text
    -a 1 or 2 bed files
    -b 1 or 2 bam files
    -c output for regions that didn't pass the threshold
    -d minimum read depth on each side
    -m minimum for highest mapq on each side
    -g max distance between 2 regions
    -t temp file
    -o output file"




while getopts ':h:a:b:c:d:m:g:o:t:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    a) bed=($OPTARG)
       ;;
    b) echo $OPTARG
       bam=($OPTARG)
       ;;
    c) cut=$OPTARG
       ;;
    d) depth="$OPTARG"
       ;;
    m) mapq="$OPTARG"
       ;;
    g) gap="$OPTARG"
       ;;
    t) temp="$OPTARG"
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



touch $temp

if [ "${bam[1]}" == '' ]
then
    awk -vOFS='\t' '{print $0, NR}' $bed | \
        bedtools intersect -bed -wb -abam $bam -b - > $temp
else
    cat ${bed[0]} ${bed[1]} | sort -k1,1 -k2,2n | bedtools merge -d $gap -i - |
        tee >(bedtools intersect -bed -wb -abam ${bam[0]} -b - >> $temp) | \
              bedtools intersect -bed -wb -abam ${bam[1]} -b - >> $temp
    #
    # cat <(bedtools intersect -a <(awk -vG=$gap '{printf "%s\t%i\t%i\n", $1, $2 - G/2, $3 + G/2}'  ${bed[0]}) \
    #                          -b <(awk -vG=$gap '{printf "%s\t%i\t%i\n", $1, $2 - G/2, $3 + G/2}'  ${bed[1]}) \
    #                          -wb -wa -loj) \
    #     <(bedtools intersect -a <(awk -vG=$gap '{printf "%s\t%i\t%i\n", $1, $2 - G/2, $3 + G/2}'  ${bed[1]}) \
    #                          -b <(awk -vG=$gap '{printf "%s\t%i\t%i\n", $1, $2 - G/2, $3 + G/2}'  ${bed[0]}) \
    #                          -wb -wa -loj | awk -vOFS='\t' '{print $4, $5, $6, $1, $2, $3}') | \
    # tee >(awk -vOFS='\t' '{
    #         if ($2==-1){
    #             seq=$4
    #             start=$5
    #             end=$6
    #         } else if ($5==-1){
    #             seq=$1
    #             start=$2
    #             end=$3
    #         } else {
    #             seq=$1
    #             start=$2<$5?$2:$5
    #             end=$6<$3?$3:$6
    #         }
    #         print seq, start, end, NR}' | \
    #     bedtools intersect -bed -wb -abam ${bam[0]} -b - >> $temp) | \
    #     awk -vOFS='\t' '{
    #             if ($2==-1){
    #                 seq=$4
    #                 start=$5
    #                 end=$6
    #             } else if ($5==-1){
    #                 seq=$1
    #                 start=$2
    #                 end=$3
    #             } else {
    #                 seq=$1
    #                 start=$2<$5?$2:$5
    #                 end=$6<$3?$3:$6
    #             }
    #             print seq, start, end, NR}' | \
    #     bedtools intersect -bed -wb -abam ${bam[1]} -b - >> $temp
fi


touch $cut


awk -F '\t' -vOFS='\t' -vD=$depth -vM=$mapq -vC=$cut '{
       pos = $13"\t"$14"\t"$15;
       if (!(pos in mapq)){
           mapq[pos] = $5
           count[pos] = 1
       } else if ($5>mapq[pos]){
           mapq[pos] = $5
           count[pos] += 1
       } else {
           count[pos] += 1
       }
    }END{
        for (pos in mapq){
            if (mapq[pos] >= M && count[pos] >= D){
                print pos, mapq[pos]
            }  else{
                print pos, mapq[pos] >> C
            }
        }}' $temp > $out


#   tee >(awk '{print $1, $2, $3, NR}' | bedtools intersect -abam $a_bam -b - ) \
#       >(awk '{print $4, $5, $6, NR}' | bedtools intersect -abam $b_bam -b - )

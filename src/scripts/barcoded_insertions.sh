

# awk '{if ($2 > 5){print $0}}' RSTP2_5.txt | sort -k2 -k3 -k4n |
awk '{if ($2 > 1){print $0}}' $1 | sort -k2 -k3 -k4n |
    awk '
    BEGIN{
        OFS="\t"
    }
    {if (NR==1){
        count=$1
        bc=$2
        tm=$3
        loc[$4] = $1
        mapq=$5 * $1
        sign=$6
    } else if (bc!=$2 || tm!=$4){
        max=0
        for (l in loc){
            if (loc[l] > max){
                location=l
                max = loc[l]
            }
        }
        print count, bc, tm, location, mapq/count, sign

        count=$1
        bc=$2
        tm=$3
        delete loc
        loc[$4] = $1
        mapq=$5 * $1
        sign=$6
    } else {
        count += $1
        mapq += ($5 * $1)
        loc[$4] += $1
    }}
    END {
        for (l in loc){
            if (loc[l] > max){
                location=l
                max = loc[l]
            }
        }
        print count, bc, tm, location, mapq, sign
    }' |
    sort -k2,2 -k1,1nr | sed 's/_/ /g' | sed 's/::/ /g' |
    awk '
    BEGIN{
        OFS="\t"
        ORS=""
        print "#barcode\tchr\tori\tpos\treads\tfreq_1\tfreq_2\tunmapped\n"
    }
    {
        if (NR==1){
            bc = $2;
            second=0;
            chr = $3
            unmapped=0
            if (chr=="*"){
                pos = "*"
                ori = "*"
                unmapped= $1
            } else {
                pos = $4 + $6
                ori = $8
            }
            print bc, chr, ori, pos, $1
        } else if (bc != $2){
            printf "\t%i\t%i\n", second, unmapped;
            bc = $2;
            second=0;
            unmapped=0;
            chr = $3
            if (chr=="*"){
                pos = "*"
                ori = "*"
                unmapped= $1
            } else {
                pos = $4 + $6
                ori = $8
            }
            print bc, chr, ori, pos, $1
        } else if (bc == $2 && second==0 && $3!="*"){
            second = $1
        } else if (bc == $2 && $3=="*"){
            unmapped = $1
        }
    } END {
        printf "\t%i\t%i\n", second, unmapped;
    }'

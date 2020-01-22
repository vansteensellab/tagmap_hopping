#!/usr/bin/awk -f
BEGIN{
    OFS="\t"
}
{
    if ($1 ~ /^[^#]/){
        print $1, $2 - spacer, $3 + spacer, $1"_"$2
    }
}

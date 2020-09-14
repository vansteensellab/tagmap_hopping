#!/usr/bin/awk -f

BEGIN {
    OFS="\t"
    arr["seqnames"] = 1
    arr["start"] = 2
    arr["end"] = 3
    arr["start_gap"] = 4
    arr["end_gap"] = 5
    arr[L] = -1
    arr[S] = -1
    arr["strand"] = -1
}
{
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

        start = $arr["start_gap"]=="."?$arr["start"]:$arr["start_gap"]
        end = $arr["end_gap"]=="."?$arr["end"]:$arr["end_gap"]

        print $arr["seqnames"], start, end, name, score, strand
    }
}

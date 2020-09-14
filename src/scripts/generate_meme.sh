#!/usr/bin/env bash
# PWM

printf $1 | awk -vOFS='\t' \
  -vPWMFile=<(printf $1) '
  BEGIN{
    c=0;
    print "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25"
    while ((getline<PWMFile)>0) {
       if($1~/^>/){
         gsub(">","",$1);
         motName=$1
      } else {
        hSizeMot[motName]++;
      }}
      close(PWMFile)
  }{
    if($0~/^>/){
      c++;
      gsub(">","",$1);
      tf=substr($1, 1, 10)
      print "\nMOTIF\t"substr($1, 1, 10)"\nletter-probability matrix: alength= 4 w= "hSizeMot[$1] " nsites= 20 E= 1"
    }
      else{
        if(NF>1){
          sum=$2+$3+$4+$5;
          printf("%.3f\t%.3f\t%.3f\t%.3f\n",$2/sum,$3/sum,$4/sum,$5/sum)
    } }
  }'

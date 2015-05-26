# Accept admix.prob.gz and admix.info files, and convert to plink dosage file and matrixeqtl input snp.txt
 
#!/bin/bash

prob2plink.pl -prob admix.prob.gz -info admix.info -o admix 

awk 'length($2)==1 && length($3)==1 && $7>0.3 && $5>=0.05 {print $1}' admix.info | LANG=en_EN sort > snps.keep 
LANG=en_EN sort admix.plink_dat > plink.sorted 
LANG=en_EN join snps.keep plink.sorted > dosage.keep 
awk '{printf "%s",$1; for(i=1;i<=133;i++){printf "\t%s",($(2*i+2)*2+$(2*i+3))};printf "\n"}' dosage.keep  > snp.txt 


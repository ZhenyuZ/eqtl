#!/bin/bash
mkdir -p ~/SCRATCH/peer/##d##/
cd ~/SCRATCH/peer/##d##/
s3cmd get --skip-existing s3://bioinformatics_scratch/zhenyu/peer/##d##/##f##
s3cmd get --skip-existing s3://bioinformatics_scratch/zhenyu/peer/##d##/evec.txt
s3cmd get --skip-existing s3://bioinformatics_scratch/zhenyu/SNP/##d##.SNP.txt
xvfb-run R CMD BATCH --no-save --no-restore '--args ##d## ##f##' ~/tools/testpeer.r 
s3cmd put ##f##.perm s3://bioinformatics_scratch/zhenyu/peer/##d##/
s3cmd put ##f##.egene s3://bioinformatics_scratch/zhenyu/peer/##d##/

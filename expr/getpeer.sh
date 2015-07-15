#!/bin/bash
# ##d##: disease name
# ##n##: number of PEER factors
mkdir -p ~/SCRATCH/peer/##d##/
cd ~/SCRATCH/peer/##d##/
s3cmd get --skip-existing s3://bioinformatics_scratch/zhenyu/peer/##d##/expr.residuals.rda
s3cmd get --skip-existing s3://bioinformatics_scratch/zhenyu/peer/##d##/clin.txt
R CMD BATCH --no-save --no-restore '--args ##d## ##n##' ~/tools/getpeer.r 
s3cmd put GE-##n##.txt s3://bioinformatics_scratch/zhenyu/peer/##d##/

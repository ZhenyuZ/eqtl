# Mach-admix imputation of chr.dat and chr.ped with 1KG EUR phase and parameters
# Only one chromosome, need to provide chromosome number i from env
# Converted to plink dosage file, and output of matrixeqtl input snp.txt
# Default: states 200, MAF 0.05, Number of samples 133

#!/bin/bash

vcfFile=$(echo "/glusterfs/users/ZHANGZ18/tools/mach/meta/CEU/chr"$i".phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.gz")
recFile=$(echo "/glusterfs/users/ZHANGZ18/tools/mach/meta/admix-paramEst/parameterEstimate/chr"$i".rec")
erateFile=$(echo "/glusterfs/users/ZHANGZ18/tools/mach/meta/admix-paramEst/parameterEstimate/chr"$i".erate")

time mach-admix --geno --dosage --probs --quality --runMode ImputeOnly \
 --crossoverMap $recFile --errorMap $erateFile -d chr.dat -p chr.ped -h \
 $vcfFile --autoflip --vcfReference --states 200 --rounds 30 --prefix admix | tee admix.log &


plink

# convert PED and FAM generated from birdseed2ped to BIM/FAM/BED
plink --ped tcga.raw.ped --map tcga.raw.map --chr 1-22 --make-bed --out tcga.auto

# estimate missing rate to determine cutoffs of next step
plink --threads 32 --bfile tcga.auto --missing --out tcga.auto

		# remove dupplicated snps
		library(data.table)
		bim = fread("tcga.auto.bim", h=F)
		miss = fread("tcga.auto.lmiss", h=T)
		# find duplicate SNPs and manually check
		bim$pos = paste(bim$V1, bim$V4, sep="-")
		dup = bim[pos %in% bim$pos[which(duplicated(bim$pos))]]
		# check missing rate and pick and higher one to remove
		miss[SNP %in% dup$V2]
		dup.snps = c("SNP_A-8387337", "SNP_A-2144818")
		write.table(dup.snps, "duplicated.snps", col.names=F, row.names=F, quote=F)

# SNP missing 0.05 and MAF 0.05 filter
plink --threads 32 --bfile tcga.auto --geno 0.05 --maf 0.05 --exclude duplicated.snps --make-bed --out tcga.auto.geno005.maf005

# Patient missing 0.05 filter
plink --threads 32 --bfile tcga.auto.geno005.maf005 --mind 0.05 --make-bed --out tcga.auto.geno005.maf005.mind005

# SNP missing 0.02 and MAF 0.05 filter
plink --threads 32 --bfile tcga.auto.geno005.maf005.mind005 --geno 0.02 --maf 0.05 --make-bed --out tcga.auto.geno002.maf005.mind005

# patient missing 0.02 filter
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind005 --mind 0.02 --make-bed --out tcga.auto.geno002.maf005.mind002


# generate prune list and prune, using data with full patient
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002 --hwe 1e-12 --make-bed --out tcga.auto.geno002.maf005.mind002.hwe1e12

plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.hwe1e12 --indep-pairwise 50 5 0.2 --out tcga.auto.geno002.maf005.mind002.hwe1e12.prune
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.hwe1e12 --extract tcga.auto.geno002.maf005.mind002.hwe1e12.prune.prune.in --make-bed --out tcga.auto.geno002.maf005.mind002.hwe1e12.prune
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.hwe1e12.prune --genome --min 0.10 --out tcga.auto.geno002.maf005.mind002.hwe1e12.prune


		# pre-processing TCGA and 1KG data before merging
		library(data.table)
		library(dplyr)
		setwd("/mnt/SCRATCH/eqtl/snp/03_population")
		kg = fread("1kg.v5a.20130502.bim", h=F)
		tcga = fread("tcga.auto.geno002.maf005.mind002.bim", h=F)

		tcga = tcga %>%
				mutate( pos = paste(V1, V4), 
						A1 = paste(V1, V4, V5, V6),
						A2 = paste(V1, V4, V6, V5)
				)

		kg = kg %>%
				mutate( pos = paste(V1, V4), 
						A1 = paste(V1, V4, V5, V6),
						A2 = paste(V1, V4, V6, V5)
				)

		common.A1 = intersect(tcga$A1, kg$A1)
		# there are only 5% flipped ref/alt, so ignored 
		# common.A2 = intersect(tcga$A2, kg$A1)

		tcga = data.table(tcga)[A1 %in% common.A1]
		kg = data.table(kg)[A1 %in% common.A1]

		write.table(tcga$V2, "tcga.common.snps", col.names=F, row.names=F, quote=F)
		write.table(kg$V2, "1kg.common.snps", col.names=F, row.names=F, quote=F)


# filter TCGA and 1KG files
plink --bfile tcga.auto.geno002.maf005.mind002 --extract tcga.common.snps --make-bed --out tcga
plink --bfile 1kg.v5a.20130502 --extract 1kg.common.snps --make-bed --out 1kg

# merge
plink --threads 32 --bfile 1kg --bmerge tcga.bed 1kg.bim tcga.fam --make-bed --out tcga.1kg.merged

# exclude LR LD
plink --bfile tcga.1kg.merged --exclude range LR-LD.range --make-bed --out tcga.1kg.merged.nolrld

# prune and prepare for smartpca
plink --bfile tcga.1kg.merged.nolrld -indep-pairwise 50 5 0.2 --make-bed --out tcga.1kg.smartpca
plink --bfile tcga.1kg.smartpca --recode --tab --out tcga.1kg.smartpca
plink --bfile tcga.1kg.smartpca --recode transpose --out tcga.1kg.smartpca

# after PCA and assign population, do hwe separately
plink --threads 32 --bfile ../03_population/tcga.auto.geno002.maf005.mind002 --keep eur.id --make-bed --out tcga.auto.geno002.maf005.mind002.eur
plink --threads 32 --bfile ../03_population/tcga.auto.geno002.maf005.mind002 --keep sas.id --make-bed --out tcga.auto.geno002.maf005.mind002.sas
plink --threads 32 --bfile ../03_population/tcga.auto.geno002.maf005.mind002 --keep afr.id --make-bed --out tcga.auto.geno002.maf005.mind002.afr
plink --threads 32 --bfile ../03_population/tcga.auto.geno002.maf005.mind002 --keep amr.id --make-bed --out tcga.auto.geno002.maf005.mind002.amr
plink --threads 32 --bfile ../03_population/tcga.auto.geno002.maf005.mind002 --keep eas.id --make-bed --out tcga.auto.geno002.maf005.mind002.eas

# make population specific dataset
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.eur --geno 0.02 --maf 0.05 --hwe 1e-6 --make-bed --out eur.auto.geno002.maf005.mind002
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.sas --geno 0.02 --maf 0.05 --hwe 1e-6 --make-bed --out sas.auto.geno002.maf005.mind002
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.eas --geno 0.02 --maf 0.05 --hwe 1e-6 --make-bed --out eas.auto.geno002.maf005.mind002
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.afr --geno 0.02 --maf 0.05 --hwe 1e-6 --make-bed --out afr.auto.geno002.maf005.mind002
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.amr --geno 0.02 --maf 0.05 --hwe 1e-6 --make-bed --out amr.auto.geno002.maf005.mind002

# get hwe stats
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.eur --hardy --out tcga.auto.geno002.maf005.mind002.eur
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.sas --hardy --out tcga.auto.geno002.maf005.mind002.sas
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.eas --hardy --out tcga.auto.geno002.maf005.mind002.eas
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.afr --hardy --out tcga.auto.geno002.maf005.mind002.afr
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.amr --hardy --out tcga.auto.geno002.maf005.mind002.amr

	# extract snps that failed hwe 1e-6 by population
	library(data.tabe)
	snp = NULL
	for(pop in c("eur", "sas", "eas", "amr", "afr")) {
		hwe = fread(paste0("tcga.auto.geno002.maf005.mind002.", pop, ".hwe"))
		snp = c(snp, hwe[P<1e-6]$SNP)
	}
	snp = unique(snp)
	write.table(snp, "snp.fail.hwe1e6.pop.txt", col.names=F, row.names=F, sep="\t", quote=F)

plink --threads 32 --bfile ../03_population/tcga.auto.geno002.maf005.mind002 --exclude snp.fail.hwe1e6.pop.txt --make-bed --out tcga.auto.geno002.maf005.mind002.hwe1e6

# IBD/IBS
plink --threads 32 --bfile ../04_hwe/tcga.auto.geno002.maf005.mind002.hwe1e6 -indep-pairwise 50 5 0.2 --make-bed --out tcga.auto.geno002.maf005.mind002.hwe1e6.prune
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.hwe1e6.prune --extract tcga.auto.geno002.maf005.mind002.hwe1e6.prune.prune.in --exclude range LR-LD.range --make-bed --out tcga.auto.geno002.maf005.mind002.hwe1e6.prune.prune
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.hwe1e6.prune.prune --genome --min 0.05 --out tcga.auto.geno002.maf005.mind002.hwe1e6.prune.prune


# generate VCFs for imputation
plink --threads 32 --bfile tcga.auto.geno002.maf005.mind002.hwe1e6 --recode vcf-iid --out tcga
bgzip -c tcga.vcf > tcga.vcf.gz
tabix -p vcf tcga.vcf.gz
tabix -h tcga.vcf.gz 22 | bgzip -c  > 22.vcf.gz

for i in `seq 1 22`; do  tabix -h tcga.vcf.gz $i | bgzip -c  > $i.vcf.gz; done


# filter by rsq >= 0.8 and maf >= 0.1 using awk
for i in `seq 7 14`; do maf=0.01; rsq=0.8; zcat chr$i.dose.vcf.gz | awk -v maf=$maf -v rsq=$rsq 'BEGIN {FS="\t|;MAF=|;R2="}; {if(substr($1,1,1)=="#") print; else if($9>maf && ($10>rsq || $7=="PASS;GENOTYPED_ONLY")) print}' | bgzip -c > ./mafrsq/$i.vcf.gz & done

# filter by hwe using vcftools
parallel -j 8 vcftools --gzvcf {1}.vcf.gz --keep ../../../03_population/TCGA.{2}.samples.txt --hardy --out {1}.{2}.pvalues ::: {1..5} ::: EUR AFR AMR EAS SAS

# combine by population hwe list and 
# remove hwe snps
parallel -j 8  vcftools --gzvcf {1}.vcf.gz --exclude-positions chr{1}.pop.hwe1e-3.to.exclude.txt --recode --recode-INFO-all --out chr{1}.filtered ::: {1..22}
parallel -j 8  vcftools --gzvcf {1}.vcf.gz --exclude-positions chr{1}.pop.hwe1e-3.to.exclude.txt --recode --recode-INFO-all --stdout | bgzip -c > chr{1}.filtered.vcf.gz ::: {1..22}
for i in `seq 1 8`; do vcftools --gzvcf $i.vcf.gz --exclude-positions chr$i.pop.hwe1e-3.to.exclude.txt --recode --recode-INFO-all --stdout | bgzip -c > chr$i.filtered.vcf.gz & done 
[1] 4179
[2] 4181
[3] 4183
[4] 4185
[5] 4187
[6] 4189
[7] 4191
[8] 4193
[9] 4195
[10] 4197
[11] 4199
[12] 4201
[13] 4203
[14] 4205


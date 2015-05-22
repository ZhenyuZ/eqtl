# Scripts to Process TCGA SNP6 birdseed genotypes downloaded as 
# Merged files from Firehose

# load packages
options(stringsAsFactors=F)
library(data.table)
source("~/github/eqtl/module.annotation.r")

disease = "KIRP"
file = paste0(disease, ".snp__genome_wide_snp_6__broad_mit_edu__Level_2__birdseed_genotype__birdseed.data.txt")

# read firehose Merged SNP6 Birdseed file
data = fread(file, skip=1)
aliquots = unlist(read.delim(file, h=F, nrows=1, stringsAsFactors=F))

nAliquot = (ncol(data)-1)/2
geno = matrix(nrow=nrow(data), ncol=0)

for (i in 1: nAliquot) {
  aliquot = aliquots[2*i]
  
  # Get Blood Normal only
  if(substr(aliquot, 14, 15) == "10") {
    call = data[, 2*i, with=F]
    confidence = data[, 2*i+1, with=F]
  
    # remove call with confidence > 0.1
    call[which(confidence > 0.1)] =  NA
    geno = cbind(geno, as.matrix(call))	
	colnames(geno)[ncol(geno)] = aliquot
  }
}  

# replace rownames with dbSNP names
conv <- read.table("/glusterfs/netapp/homes1/ZHANGZ18/meta/SNP_A2dbsnp.txt", as.is=T)
rownames(geno) <- conv$V2[match(data$"Composite Element REF", conv$V1)]

# remove SNPs without "rs" names
geno <- geno[which(grepl("^rs", rownames(geno))), ]

# load TCGA annotation table
load("~/eqtl/annotation.tbl.rda") 

# filter by annotation table 
aliquots <- colnames(geno)
aliquots <- FilterByAnnotation(aliquots, annot.tbl)
geno <- geno[, which(colnames(geno) %in% aliquots)]

# save the processed dbSNP matrix
save(geno, file=paste0(disease, ".geno.rda"))
  
# save aliquot names
write.table(colnames(geno), paste0(disease, ".geno.aliquots"), col.names=F, row.names=F, quote=F, sep="\t")


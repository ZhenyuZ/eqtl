# SplitDosage.r filter dosage file by plink-filtered snps, and split 
# dosage files by disease types

options(stringsAsFactors=F)
library(data.table)
diseases = unlist(read.table("~/SCRATCH/plink2/diseases.txt"))

# read meta data
info = fread("~/SCRATCH/plink2/tcga.snp.info", h=T)
tcga.fam = fread("~/SCRATCH/plink2/tcga.geno.fam", h=F)
ids = tcga.fam$V2

# get patient.lst, which contains patient barcodes of each disease
# get snp.lst, which contains QC'ed snp list for each disease genotypes
# also write patient by disease as header
n = length(diseases)
patient.lst = array(list(), n)
snp.lst = array(list(), n)
for(i in 1: n) {
	disease = diseases[i]
	# read plink bim file to get snps
  bim = fread(paste0("~/SCRATCH/plink2/", disease, ".bim"), h=F)
  snp.lst[[i]] = bim$V2
  # read plink fam file to get patient ids
  fam = fread(paste0("~/SCRATCH/plink2/", disease, ".fam"), h=F)
  patient = fam$V2
  patient.lst[[i]] = patient
  # write header
  outfile = paste0("~/SCRATCH/SNP/", disease, ".header")
  header = paste(c("id", patient), collapse="\t")
  write.table(header, outfile, col.names=F, row.names=F, quote=F, sep="\t")
}
save(info, ids, diseases, patient.lst, snp.lst, file = "~/SCRATCH/SNP/meta.rda")


##################################################################
# second part: can run separatedly from obove, need load meta.rda
##################################################################
options(stringsAsFactors=F)
library(data.table)
load("~/SCRATCH/SNP/meta.rda")

# make dosage matrix file by chromosome, change colnames and rownames
# as SNP id and patient id, respectively
for (i in 1:22) {
  # read dosage
  dosage.file = paste0("~/SCRATCH/impute2/dosage.chr/chr", i, ".dosage")
  dosage = fread(dosage.file, h=F)

  # filter by existence in info 
  snp = dosage$V1
  m = match(snp, info$impute2.snp)
  snp = info$snp[m]
  w = which(!is.na(snp))
  snp = snp[w]

  # make dosage as matrix
  dosage = as.matrix(dosage[w, 2:ncol(dosage), with=F])
  rownames(dosage) = snp
  colnames(dosage) = ids

  # save dosage by chromosome r data file
  outfile = paste0("~/SCRATCH/impute2/trimed.dosage/chr", i, ".dosage.rda")
  save(dosage, file=outfile)
}


##################################################################
# third part: can run separatedly from obove, need load meta.rda
##################################################################
options(stringsAsFactors=F)
library(data.table)
load("~/SCRATCH/SNP/meta.rda")

# split dosage by disease type
for(chr in 1:22) {
  # load dosage file
  cat("loading Chromosome", chr, "... ... \n")
  dosage.file = paste0("~/SCRATCH/impute2/trimed.dosage/chr", chr, ".dosage.rda")
  load(dosage.file)

  for (i in 1:length(diseases)) {
  	disease = diseases[i]
  	cat("  processing disease type", disease, "... \n")
  	patients = patient.lst[[i]]
  	snps = snp.lst[[1]]

    # trim by id
    m = match(patients, colnames(dosage))
    my.dosage = dosage[, m]
    
    # trim by snp
    w = which(rownames(dosage) %in% snps)
    my.dosage = my.dosage[w, ]

    # filter again by MAF 0.05
    af = apply(my.dosage, 1, mean)/2
    w = which(af >= 0.05 & af <= 0.95)
    my.dosage = my.dosage[w, ]

    # output
    cat("    writing output to hard disk\n")
    outfile = paste0("~/SCRATCH/SNP/", disease, ".chr", chr, ".SNP.txt")
    write.table(my.dosage, outfile, col.names=F, row.names=T, quote=F, sep="\t")
  }
}



    

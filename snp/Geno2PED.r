# Geno2PED.r Read tcga.snp.annotation.txt for SNP annotation
# and convert geno.rda data into plink PED format

options(stringsAsFactors=F)

# read diseases
diseases = unlist(read.delim("diseases.txt", h=F))

# read SNP annotation file
annot = read.delim("./meta/tcga.snp.annotation.txt", h=T)

# read gender and vitalstatus information from clinical data
clin = data.frame(read.delim("./aliquot/clin.summary.txt", h=T))
 
# function to change 0 1 2 to ATGC 
Geno2ATGC <- function(geno) {
  # loop to change 0 1 2 coding to ATGC coding
  for(i in 1: nrow(geno)) {
    g = geno[i, ]
    g[which(is.na(g))] = "0 0"
    g[which(g=="0")] = annot$code0[i]
    g[which(g=="1")] = annot$code1[i]
    g[which(g=="2")] = annot$code2[i]
    geno[i,] = g
  }
  return(geno)
}


for(disease in diseases) {  
  # Get Genotype data
  geno.file = paste0("./snp/", disease, ".geno.rda")
  load(geno.file)

  # Convert to ATGC
  geno = Geno2ATGC(geno)
  
  # Add fid, iid, pid, mid, sex and phenotype
  fid = colnames(geno)
  iid = substr(fid, 1, 12)
  pid = rep(disease, length(fid))
  mid = rep("0", length(fid))
  m = match(iid, clin$patient)
  sex = clin$gender[m]
  phenotype = clin$dead[m]
  
  tped = Reduce(rbind, list(fid,iid,pid,mid,sex,phenotype, geno))
  
  ped.file = paste0("./plink/", disease, ".ped")
  write.table(t(tped), ped.file, col.names=F, row.names=F, quote=F, sep=" ")
}

  

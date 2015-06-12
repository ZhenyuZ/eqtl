# seeQTL.summary.r input seeQTL dataset, calculate beta of 
# Churc1 expression ~ SNP (minor allele = 1, major allele =0)

setwd("/Users/zhenyu/Downloads/genotype_expression_data")

# load data list
suffix = unlist(read.table("suffix.txt", h=F, stringsAsFactors=F))

# load snp summary
load("summary.rda")
snps = rownames(summary)

# init output
output = matrix(NA, length(snps), 2*length(suffix))
rownames(output) = snps

for(i in 1:length(suffix)) {
  s = suffix[i]

  # read edata
  edata = read.delim(paste0("matrix_edata", s), h=T, row.names=1, stringsAsFactors=F)
  if(! "91612" %in% rownames(edata)) break
  edata = unlist(edata["91612", ])
  
  # read mdata  
  mdata = read.delim(paste0("matrix_mdata", s), h=T, row.names=1, stringsAsFactors=F)
  m = match(snps, rownames(mdata))
  if(sum(is.na(m)) == 11) break
  mdata = mdata[m,]
  
  # change to minor alleles
  row.of.major.allele = which(apply(mdata, 1, mean, na.rm=T) > 1.2)
  mdata[row.of.major.allele, ] = 2 -   mdata[row.of.major.allele, ]
  output[, 2*i-1] = apply(mdata, 1, mean, na.rm=T)/2
 
  
  # linear regression edata ~ mdata
  mdata = mdata[which(!is.na(m)), ]
  beta =  apply(mdata, 1, function(x) lm(edata ~ unlist(x))$coefficients[2])

  # write to output
  output[which(!is.na(m)), 2*i] = beta
}

# Add colnames
for (i in 1: length(suffix)) {
  colnames(output)[2*i-1] = "maf"
  colnames(output)[2*i] = suffix[i]
}

write.table(output, "seeQTL.direction.txt", col.names=T, row.names=T, sep="\t", quote=F)




# ProcessGeneLevelMethylation.r take firebrowse meth.by_mean.data.txt file,
# filter by aliquot id and order by RNA-SeqV2 gene names, and output 
# met.by.gene.txt files for each disease type

options(stringsAsFactors=F)
library(data.table)

# read disease types
diseases = unlist(read.table("./diseases.txt", h=F))

expr.sample.file = "./met2/unc.edu.006eb1bc-b058-49b9-8c2d-3c50f2427d95.1214207.rsem.genes.results"

# read rnaseq-V2 sample file to extract gene names
expr = fread(expr.sample.file, h=T)
id = expr$gene_id
gene = sapply(id, function(x) strsplit(x, "\\|")[[1]][1])
entrez = sapply(id, function(x) strsplit(x, "\\|")[[1]][2])

for (disease in diseases) {
  # read met by mean file from firebrowse
  met.mean.file = paste0("./met2/", disease, ".meth.by_mean.data.txt")
  met = fread(met.mean.file, h=T)
  
  # read common patient
  common.patient.file = paste0("./plink/geno/", disease, ".patient")
  common.patient = unlist(read.table(common.patient.file, h=F))
  
  # filter by aliquot
  met.patient = substr(colnames(met), 1, 12)
  setnames(met, colnames(met), met.patient)
  m = match(common.patient, met.patient)
  if(sum(is.na(m)) > 0) print("WARNING! missing data")
  met.gene = met[, 1, with=F]
  met = met[, m, with=F]
  
  # re-arrange gene order
  m = match(gene, unlist(met.gene))
  output = as.matrix(met[m, ])
  rownames(output) = id
  mode(output) = "numeric"
    
  # output 
  outFile = paste0("./met2/", disease, ".met.by.gene.txt")
  write.table(output, outFile, col.names=T, row.names=T, sep="\t", quote=F)
  outFile2 = paste0("./met2/", disease, ".met.by.gene.rda")
  met = output
  save(met, file=outFile2)
}

  
  
    

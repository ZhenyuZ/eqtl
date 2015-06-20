# CalculateGeneLevelMethylation.r takes matrix methylation matrix, methylation 
# table calculated by MapMet2Gene.r, and output gene level methylation of all 
# disease types


# Function to impute NA in a numeric matrix.  If a whole row is NA, fill with 
# mean value of the matrix; if a row is partial NA, fill with row mean  
ImputeMatrix <- function(data) {
  if(mode(data) != "numeric") {
    print("Warning! ImputeMatrix input is not a numeric matrix.")
  }
  # get data mean and row means
  data.mean = mean(data, na.rm=T)  
  row.mean = apply(data, 1, function(x) mean(x, na.rm=T))
  
  # fix whole row of NA
  data[which(is.nan(row.mean)), ] = data.mean

  # fix partial row of NA
  w = which(is.na(data))
  data[w] = row.mean[w %% nrow(data)]
  
  # return imputed data
  return(data)
}   


options(stringsAsFactors=F)
library(data.table)

diseases = unlist(read.table("diseases.txt", h=F))
met27.table = data.frame(read.table("./meta/met.table.27.txt", h=T))
met450.table = data.frame(read.table("./meta/met.table.450.txt", h=T))

summary.file = "./aliquot/euro.aliquot.summary.txt"
summary = fread(summary.file, h=T)

for (disease in diseases) {
  print(disease)
  if (disease == "OV") {
    met.table = met27.table
  } else {
    met.table = met450.table
  }
    
  patient.file = paste0("./plink/geno/", disease, ".patient")
  met.file = paste0("./met/", disease, ".raw.met.txt")
  log.file = paste0(disease, ".log")
  out.file = paste0("./met/", disease, ".processed.met.txt")
  
  # read met data
  print("  reading methylation by probe data")
  met = fread(met.file, h=T)
  
  # filter by common patient
  print("  trimming methylation data by patient")
  patient = unlist(read.table(patient.file, h=F))
  m = match(patient, summary$patient)
  if(sum(is.na(m)) > 0) {
    cat("missing patient in summary file", file=log.file, append=T)
  }
  aliquot = summary$Met[m]
  m = match(aliquot, colnames(met))
  if(sum(is.na(m)) > 0) {
    cat("missing patient in met file", file=log.file, append=T)
  }  
  met = met[, m, with=F]  
  met = as.matrix(met)  
  colnames(met) = patient
  n = length(patient)

  # impute 
  print("  imputation of missing methylation value ")
  met = ImputeMatrix(met)

  # calculate gene level methylation
  out = matrix(0, 20531, n)
  colnames(out) = patient
  for (i in 1: nrow(met.table)) {
    if(i %% 10000 == 0) {
      print(paste("    ", i, "out of", nrow(met.table), "done"))
    }  
    met.row = met.table$met.row[i]
    expr.row = met.table$expr.row[i]
    factor = met.table$count[i]
    out[expr.row, ] = out[expr.row, ] + met[met.row, ]/factor
  }
  rownames(out) = 1:20531
  
  # output 
  print("  writing result to disk")
  write.table(out, out.file, col.names=T, row.names=T, sep="\t", quote=F)
}  
  
    
 


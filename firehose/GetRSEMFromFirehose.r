# Input disease and firehose normalized rsem expression data
# Output all expression, expression of primary tumor and normal
# tissue, and corresponding aliquots. 
GetRSEMFromFirehose <- function(disease, file) {
  # load packages
  options(stringsAsFactors=F)
  library(data.table)
  source("~/github/eqtl/module.annotation.r")
  
  # init 
  disease = toupper(disease)

  # read firehose rsem.gene.normalized
  expr = fread(file, h=T)
  data = expr[2: nrow(expr), 2: ncol(expr), with=F]
  data = as.matrix(data)
  mode(data) = "numeric"
  rownames(data) = unlist(expr[2:nrow(expr), 1, with=F])

  # save data
  save(data, file=paste0(disease, ".rsem.rda"))

  # load TCGA annotation table
  load("~/eqtl/meta/annotation.tbl.rda") 

  # filter by annotation table 
  aliquots <- colnames(data)
  aliquots <- FilterByAnnotation(aliquots, annot.tbl)
  data <- data[, which(colnames(data) %in% aliquots)]

  # save solid tumor data 
  aliquots <- colnames(data)
  tp = aliquots[substr(aliquots, 14, 15) == "01"]
  expr = log2(data[, which(colnames(data) %in% tp)]+1)
  save(expr, file=paste0(disease, ".tp.expr.rda"))
  write.table(tp, paste0(disease, ".tp.aliquots"), col.names=F, row.names=F, quote=F)

  # save solid normal data 
  nt = aliquots[substr(aliquots, 14, 15) == "11"]
  expr = log2(data[, which(colnames(data) %in% nt)]+1)
  save(expr, file=paste0(disease, ".nt.expr.rda"))
  write.table(nt, paste0(disease, ".nt.aliquots"), col.names=F, row.names=F, quote=F)

}



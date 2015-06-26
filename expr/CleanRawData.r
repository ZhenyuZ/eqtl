# CleanRawData.r take met, cnv, and rnaseqv2 data, and make 
# numeric matrix with the same gene list and patient names

options(stringsAsFactors=F)
library(data.table)

diseases = unlist(read.table("./diseases.txt"))

for (disease in diseases) {

# need indent 
summary.file = "./aliquot/euro.aliquot.summary.txt"
expr.file = paste0("./rnaseqv2/", disease, ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt")
met.file = paste0("./met2/", disease, ".met.by.gene.rda")
cnv.file = paste0("./cnv/", disease, ".CN.txt")
out.file = paste0("./expr/", disease, "/expr.met.cnv.raw.rda")
patient.out.file = paste0("./expr/", disease, "/patient")
aliquot.out.file = paste0("./expr/", disease, "/expr.aliquot")

# read summary
summary = data.frame(read.table(summary.file, h=T, colClasses="character", sep="\t", comment.char=""))
patient = summary$patient[which(summary$disease==disease)]
expr.aliquot = summary$Expr[which(summary$disease==disease)]

# read first line of rnaseqv2 dataset
header = readLines(expr.file, n=1)
header = strsplit(header, "\t")[[1]]
header = header[2: (length(header)-1)]

# read the whole rnaseqv2 dataset
expr = read.table(expr.file, skip=2, row.names=1, colClasses=c("character", rep("numeric", length(header))), sep="\t")

# filter expr by common patient
expr.aliquot = expr.aliquot[expr.aliquot %in% header]
patient = substr(expr.aliquot, 1, 12)
m = match(expr.aliquot, header)
if(sum(is.na(m)) > 0) print("missing expr aliquot")
expr = as.matrix(expr[, m])
colnames(expr) = expr.aliquot
gene = rownames(expr)

# read met data
load(met.file)
m = match(patient, colnames(met))
if(sum(is.na(m)) > 0) print("wrong methylation order")
met = met[, m]

# read cnv data
cnv = read.table(cnv.file, h=T, sep="\t", check.names=F, colClasses="numeric")
m = match(patient, colnames(cnv))
if(sum(is.na(m)) > 0) print("wrong copy number order")
cnv = as.matrix(cnv)
cnv = cnv[, m]
rownames(cnv) = gene

# output
write.table(patient, patient.out.file, col.names=F, row.names=F, quote=F)
write.table(expr.aliquot, aliquot.out.file, col.names=F, row.names=F, quote=F)
save(expr, cnv, met, expr.aliquot, patient, file=out.file)

}


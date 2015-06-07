# GetMatchedAliquots is a script to find matched aliquots from:
# germline birdseed genotype
# somatic RnaseqV2 expression
# somatic CNA
# somatic methylation
# The scripts read all diseases from "../diseases.txt"
# aliquots informations from
# "./snp/[disease].geno.aliquots"
# "./rnaseqv2/[disease].rnaseqv2.aliquots"
# "./met/[disease].met.aliquots"
# "./cnv/[disease].cnv.aliquots"
# Then it will find common portions among somatic aliquots, and common
# patients among somatic and germline. 
# The output will be written to "../aliquot/[disease].matched.aliquots"


options(stringsAsFactors=F)


# read disease list
diseases = unlist(read.table("diseases.txt", h=F))
s = data.frame(matrix(0, length(diseases), 6))
colnames(s) = c("disease", "num.SNP", "num.Expr", "num.Met", "num.CNA", "num.somatic.common")


for (i in 1: length(diseases)) {
  disease = diseases[i]
  
  # process somatic rnaseqv2 aliquots
  expr.file = paste0("./rnaseqv2/", disease, ".rnaseqv2.aliquots")
  expr.aliquots = unlist(read.table(expr.file, h=F))
  expr.aliquots = expr.aliquots[which(substr(expr.aliquots, 14,15)=="01")]
  expr.portions = substr(expr.aliquots, 1, 19)
  
  # process somatic methylation450k aliquots
  met.file = paste0("./met/", disease, ".met.aliquots")
  met.aliquots = unlist(read.table(met.file, h=F))
  met.aliquots = met.aliquots[which(substr(met.aliquots, 14,15)=="01")]
  met.portions = substr(met.aliquots, 1, 19)

  # process somatic CNV aliquots
  cnv.file = paste0("./cnv/", disease, ".cnv.aliquots")
  cnv.aliquots = unlist(read.table(cnv.file, h=F))
  cnv.aliquots = cnv.aliquots[which(substr(cnv.aliquots, 14,15)=="01")]
  cnv.portions = substr(cnv.aliquots, 1, 19)
  
  # common tumor portion and patient
  common.tumor.portions = Reduce(intersect, list(expr.portions, met.portions, cnv.portions))
  tumor.patients = substr(common.tumor.portions, 1, 12)
  
  # process germline genotype aliquots
  geno.file = paste0("./snp/", disease, ".geno.aliquots")
  geno.aliquots = unlist(read.table(geno.file, h=F))
  geno.aliquots = geno.aliquots[which(substr(geno.aliquots, 14,15)=="10")]
  normal.patients = substr(geno.aliquots, 1, 12)  
  
  # common patients and aliquots
  common.patients = intersect(normal.patients, tumor.patients)  
  tumor.portions = common.tumor.portions[match(common.patients, tumor.patients)]
  
  # Aggregate patients and aliquots
  output = data.frame(matrix("", length(common.patients), 10))
  colnames(output) = c("Patient", "SNP", "Expr", "Met", "CNA")
  output$patients = common.patients
  output$SNP = geno.aliquots[match(common.patients, normal.patients)]
  output$Expr = expr.aliquots[match(tumor.portions, expr.portions)] 
  output$Met = met.aliquots[match(tumor.portions, met.portions)] 
  output$CNA = cnv.aliquots[match(tumor.portions, cnv.portions)] 
  
  s$disease[i] = disease
  s$num.SNP[i] = length(geno.aliquots)
  s$num.Expr[i] = length(expr.aliquots)
  s$num.Met[i] = length(met.aliquots)
  s$num.CNA[i] = length(cnv.aliquots)
  s$num.somatic.common[i] = length(tumor.patients)
  
  # output 
  out.file = paste0("./aliquot/", disease, ".matched.aliquots")
  write.table(output, out.file, col.names=T, row.names=F, sep="\t", quote=F)
}
  
# output summary
out.file = paste0("./aliquot/number.summary.txt")
write.table(s, out.file, col.names=T, row.names=F, sep="\t", quote=F)
  

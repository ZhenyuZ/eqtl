# Extract TP aliquots from Firehose SNP6 CNV Merged data file
# write output as [disease].cnv.aliquots

# load packages
options(stringsAsFactors=F)
library(data.table)
source("~/github/eqtl/module.annotation.r")

# init
diseases = read.delim("../diseases.txt", h=F)
diseases = unlist(diseases)

# load annotation table
load("~/eqtl/annotation.tbl.rda") 

# Do for each disease
for (disease in diseases) {
  file = paste0(disease, ".snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
  data = fread(file)
  aliquots = unique(data$Sample)
  aliquots = aliquots[which(substr(aliquots,14,15)=="01")]
  aliquots <- FilterByAnnotation(aliquots, annot.tbl)
  write.table(aliquots, paste0(disease, ".cnv.aliquots"), col.names=F, row.names=F, quote=F, sep="\t")
}


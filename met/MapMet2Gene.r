# MapMet2Gene.r take sample level 3 methylation 27k or 450k file, and sample 
# RNA-SeqV2 gene quantification file, calculate a met.table file, with four 
# fields: 
#   gene: probe affected gene name
#   met.row: probe row number in methylation file
#   expr.row: corresponding gene row number in expression file
#   count: number of appearance of this gene

options(stringsAsFactors=F)

# read sample met level3 file
met = data.frame(read.table("jhu-usc.edu_LUSC.HumanMethylation27.4.lvl-3.TCGA-07-0227-20A-01D-1096-05.txt", h=T, stringsAsFactors=F, sep="\t", skip=1))

# make of table of methylation affected gene, and the row index of the probe
met.table = data.frame(matrix(NA, 2*nrow(met), 2))
colnames(met.table) = c("gene", "met.row")
genes = met$Gene_Symbol
index = 1
for (i in 1: nrow(met)) {
  if(i %% 100 == 0) print(i)
  gs = genes[i]
  if(nchar(gs) == 0 | is.na(gs)) next
  gs = strsplit(gs, ";")[[1]]
  for(g in gs) {
    met.table[index, ] = c(g, i)
    index = index + 1
  }
}
met.table = met.table[1: (index-1), ]

# read sample rnaseqv2 level3 file
expr = data.frame(read.table("unc.edu.006eb1bc-b058-49b9-8c2d-3c50f2427d95.1214207.rsem.genes.results", h=T, stringsAsFactors=F))

# extract gene name
expr$gene = sapply(expr$gene_id, function(x) strsplit(x, "\\|")[[1]][1])
expr$entrez = sapply(expr$gene_id, function(x) strsplit(x, "\\|")[[1]][2])

# match gene position
genes = met$Gene_Symbol
m = match(met.table$gene, expr$gene)
met.table$expr.row = m 
met.table = met.table[-which(is.na(m)), ]

# count gene appearance
count = table(met.table$gene)
count = count[match(met.table$gene, names(count))]
met.table$count = count

# output
write.table(met.table, "met.table.txt", col.names=T, row.names=F, quote=F, sep="\t")
save(met.table, file="met.table.27.rda")


    


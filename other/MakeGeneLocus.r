# Input TCGA gaf2 description file, and output genelocs.txt

# init
options(stringsAsFactors=F)

# read example example.rsem.genes.results 
example = read.delim("~/eqtl/meta/example.rsem.genes.results", h=T)

# read gaf and extract gene only
gaf = read.delim("TCGA.hg19.June2011.gaf", h=T)
gaf = gaf[gaf$FeatureType == "gene", ]

# find and filter genes with multiple gene range definition in gaf2
w = which(grepl(";", gaf$GeneLocus))
gaf = gaf[-w, ]

# keep only genes that exist in the rsem.genes.results
comm = intersect(example$gene_id, gaf$Gene)
gaf = gaf[match(comm, gaf$Gene), ]

# output init
output = data.frame(matrix("", nrow(gaf), 9))
colnames(output) = c("geneid", "chr", "s1", "s2", "symbol", "entrez", "strand", "gene.length", "transcript.length")

# add gene names
output$geneid = gaf$Gene
output$symbol = sapply(output$geneid, function(x) strsplit(x, "\\|")[[1]][1])
output$entrez = sapply(output$geneid, function(x) strsplit(x, "\\|")[[1]][2])

# add strand 
output$strand = gaf$strand

# add gene location
locus = gaf$GeneLocus
gaf$chr = sapply(locus, function(x) strsplit(x, ":")[[1]][1])
gaf$strand = sapply(locus, function(x) strsplit(x, ":")[[1]][3])
locus = sapply(locus, function(x) strsplit(x, ":")[[1]][2])
output$s1 = sapply(locus, function(x) strsplit(x, "-")[[1]][1])
output$s2 = sapply(locus, function(x) strsplit(x, "-")[[1]][2])

# calculate gene length
output$gene.length = -sapply(locus, function(x) eval(parse(text=x))) + 1
# sanity check summary(output$gene.length)

# calculate transcript length
output$transcript.length = sapply(gaf$FeatureCoordinates, function(x) {y=strsplit(x, "-")[[1]]; y[length(y)]})

write.table(output, "~/eqtl/meta/genelocs.txt", col.names=T, row.names=F, quote=F, sep="\t")


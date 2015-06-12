# ExtractAffyAnnot is a script to read SNP6 array annotation file
# GenomeWideSNP_6.na[xx].annot.csv file and extract the most 
# important information to create a smaller file na35.CEU.txt

options(stringsAsFactors=F)
out.file = "./meta/na35.CEU.txt"

# read affy annotation file
annot.file = "./meta/GenomeWideSNP_6.na35.annot.csv"
annot = read.csv(annot.file, h=T, comment.char="#",colClasses="character")

# Extract ID information
output = data.frame(annot$dbSNP.RS.ID)
colnames(output) = "dbSNP"
output$probeset = annot$Probe.Set.ID

# Extract position information
output$chr = annot$Chromosome
output$pos = annot$Physical.Position
output$strand = annot$Strand
output$distance = sapply(annot$Genetic.Map, function(x) strsplit(x, " //")[[1]][1])

# Extract B Allele information
output$Allele.A = annot$Allele.A
output$Allele.B = annot$Allele.B
output$CEU.Freq.B = sapply(annot$Allele.Frequencies, function(x) strsplit(x, " //")[[1]][1])
output$Strand.Versus.dbSNP = annot$Strand.Versus.dbSNP

# Extract Minor Allele information
output$CEU.Minor.Allele = sapply(annot$Minor.Allele, function(x) strsplit(x, " //")[[1]][1])
output$CEU.Minor.Allele.Frequency = sapply(annot$Minor.Allele.Frequency, function(x) strsplit(x, " //")[[1]][1])

# write output
write.table(output, out.file, col.names=T, row.names=F, quote=F, sep="\t")

  
  

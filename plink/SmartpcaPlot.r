# SmartpcaPlot.r take smartpca output evec file, and plot it

# define input
evec.file = "tcga.hapmap3.prune.evec"
pop.file = "relationships_w_pops_121708.txt"

# read smartpca evec output
evec = data.frame(read.table(evec.file,  h=F , stringsAsFactors=F, row.names=1, comment.char="#"))

# clean evec 
evec = evec[, -ncol(evec)]
colnames(evec) = paste0("PC", 1:ncol(evec)) 
iid = sapply(rownames(evec), function(x) strsplit(x, ":")[[1]][1])
rownames(evec) = iid
 
# read 1000G pop info
pop.info = data.frame(read.table(pop.file, h=T, stringsAsFactors=F)) 
 
# generate population data
pop = rep("", nrow(evec))
m = match(rownames(evec), pop.info$IID)
pop = pop.info$population[m]
pop[which(substr(rownames(evec), 1, 4) == "TCGA")] = "TCGA"
evec$population = pop

# define color table
color.table = data.frame(matrix("", 12, 2))
colnames(color.table) = c("population", "color")
color.table$population = c("ASW","CHB","CHD","GIH","JPT","LWK","MEX","MKK","YRI","CEU","TSI","TCGA")
color.table$color = c("red","lightgreen","green","orange","darkgreen","darkred","purple","pink","darkred","blue","lightblue","#BEBEBE80")

# assign color
color = array("", length(pop))
for (i in 1: nrow(color.table)) {
  color[which(pop == color.table$population[i])] = color.table$color[i]
}
evec$color = color

png("tcga.hapmap.smartpca.PC1-2.png", width=960, height=720)
with(evec, plot(PC1, PC2, col=color, pch=19, main="Genotype PCA of combined TCGA and Hapmap data", cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5))
legend("bottomleft", color.table$population, col = color.table$color, pch = 19, cex=1.5)
dev.off()

png("tcga.hapmap.smartpca.PC3-4.png", width=960, height=720)
with(evec, plot(PC3, PC4, col=color, pch=19, main="Genotype PCA of combined TCGA and Hapmap data", cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5))
legend("bottomleft", color.table$population, col = color.table$color, pch = 19, cex=1.5)
dev.off()

# Based on the plots above, determine which TCGA patients are European 
# European descendant PC2 > 0.012 - 2 * PC1 ; PC4 > -0.04 - 9 * PC3
evec$euro = with(evec, (PC2 > 0.012 - 2 * PC1) & (PC4 > -0.04 - 9 * PC3))
select = with(evec, euro & population=="TCGA")
tcga.evec = evec[select,]
write.table(tcga.evec, "../tcga.evec.txt", col.names=T, row.names=T, quote=F, sep="\t")
write.table(evec, "../tcga.hapmap3.txt", col.names=T, row.names=T, quote=F, sep="\t")


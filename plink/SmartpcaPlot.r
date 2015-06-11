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

with(evec, plot(PC1, PC2, col=color, pch=19))

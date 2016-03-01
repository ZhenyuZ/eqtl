options(stringsAsFactors=F)
library(data.table)

# read disease list
diseases = read.table("../all.disease", h=F)$V1

# init
data = array(list(), length(diseases))
names(data) = diseases

# Bin all positions in (-range, +range) with windows size, and return number in each bin
GetCount = function(pos, window = 10000, range = 1000000) {
	# int
	half.bin = ceiling(range / window)
	bins = -half.bin : half.bin
	count = rep(0, length(bins))
	names(count) = bins

	# bin
	index = floor((pos + window / 2) / window)

	# count
	tbl = table(index)
	m = match(names(tbl), bins)
	count[m] = tbl
	return(count)
}

# Bin distances of each distance data
for(i in 1: length(diseases)) {
	# load distance data
	print(disease)
	disease = diseases[i]
	load(paste0(disease, ".gene.snp.distance.rda"))

	# summarize disease data
	egene.count = GetCount(egene.pos)
	index = as.numeric(names(egene.count))
	disease.data = data.frame(index)
	disease.data$egene = egene.count / sum(egene.count)
	non.egene.count = GetCount(non.egene.pos)
	disease.data$non.egene.prop = non.egene.count / sum(non.egene.count)
	all.pair.count = GetCount(all.pair.pos)
	disease.data$all = all.pair.count / sum(all.pair.count)

	data[[i]] = disease.data
}
save(data, file="snp.distance.bin.rda")

# Plot seperate graphs
png(file = "eq.distance.png", width=2160, height=5120) 
par(mfrow=c(7, 3))
par(mar=c(6.1, 5.1, 3.1, 1.1), mgp=c(1,1,1))
par(cex.main=5, cex.lab=3, cex.axis=2.5, lwd=4)

for(i in 1: length(data)) {	
	mydata = data[[i]]
	disease = names(data)[i]
	print(disease)

	# Plot eGene Best SNPs
	plot(index, mydata$egene, pch="", axes=F, ylim=c(0, 0.2), xlim=c(-100, 100), xlab ="Distance to TSS (kb)", ylab="Fraction of eQTLs", main=disease)
	lines(index, mydata$egene, lwd=3, col="blue")
	axis(1, pos=0, at=c(-100, -50, 0, 50, 100), labels=c("-1000", "-500", "0", "500", "1000"), tck=0.01)
	axis(2, pos=min(index), at=c(0, 0.05, 0.1, 0.15, 0.2), tck=0.01)

	# Plot non-eGene Best SNPs
	par(new=T)
	plot(index, mydata$non.egene, pch="", axes=F, ylim=c(0, 0.2), xlim=c(-100, 100), xlab="", ylab="", main="")
	lines(index, mydata$non.egene, col="red")

	# Plot All SNPs
	par(new=T)
	plot(index, mydata$all, pch="", axes=F, ylim=c(0, 0.2), xlim=c(-100, 100), xlab="", ylab="", main="")
	lines(index, mydata$all, col="black")	
}
dev.off()


# Get cross-disease summary data
egene = rep(0, length(index))
names(egene) = index
non.egene = egene
all = egene
for(i in 1: 18) {	
	egene = egene + data[[i]]$egene
	non.egene = non.egene + data[[i]]$non.egene
	all = all + data[[i]]$all
}
egene = egene / 18
non.egene = non.egene / 18
all = all / 18
mydata = data.frame(cbind(index, egene, non.egene, all))
save(mydata, file="summary.snp.distance.bin.rda")


# Plot summary distance plot
png(file = "eq.summary.distance.png", width=960, height=960) 
par(mfrow=c(1, 1))
par(mgp=c(1,0.5,1))
par(cex.lab=1.2, lwd=3)

plot(mydata$index, mydata$egene, pch="", axes=F, ylim=c(0, 0.16), xlim=c(-100, 100), xlab ="Distance to TSS (kb)", ylab="Fraction of eQTLs")
lines(mydata$index, mydata$egene, col="blue")
axis(1, pos=0, at=seq(-100, 100, 20), labels=seq(-1000, 1000, 200), tck=0.01)
axis(2, pos=min(mydata$index), at=c(0, 0.05, 0.1, 0.15, 0.2), tck=0.01)

# Plot non-eGene Best SNPs
par(new=T)
plot(mydata$index, mydata$non.egene, pch="", axes=F, ylim=c(0, 0.2), xlim=c(-100, 100), xlab="", ylab="", main="")
lines(mydata$index, mydata$non.egene, col="red")

# Plot All SNPs
par(new=T)
plot(mydata$index, mydata$all, pch="", axes=F, ylim=c(0, 0.2), xlim=c(-100, 100), xlab="", ylab="", main="")
lines(mydata$index, mydata$all, col="black")	

# add legend
c1 = "Significant eQTLs, FDR<10%"
c2 = "Non-ignificant eQTLs, FDR<10%"
c3 = "All SNP-gene pairs tested (Null)"
legend("topright", c(c1, c2, c3), lwd = rep(4,3), col = c("blue", "red", "black"), cex=1.2)

dev.off()











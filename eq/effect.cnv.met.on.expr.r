# Read Data
load("expr.met.cnv.processed.rda")

# Filter gene with mean expression < 5
mean.expr = apply(expr, 1, mean)
w = which(mean.expr >= 5)
expr = expr[w, ]
cnv = cnv[w, ]
met = met[w, ]

# CNV Raw
par(mfrow=c(2,2), mgp=c(2.7,1.5,1))
plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1,"Copy Number\nAlteration", cex=3)
plot(density(cnv[100, ]), xlab="Segment Mean", col="blue", lwd=3, main="One gene of all samples")
plot(density(cnv[, 100]), xlab="Segment Mean", col="blue", lwd=3, main="All gene of one samples")
plot(density(cnv), xlab="Segment Mean", col="blue", lwd=3, main="All genes of all samples")

# MET Raw
par(mfrow=c(2,2), mgp=c(2.7,1.5,1))
plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1,"Methylations", cex=3)
plot(density(met[100, ]), xlab="Segment Mean", col="blue", lwd=3, main="One gene of all samples")
plot(density(met[, 100]), xlab="Segment Mean", col="blue", lwd=3, main="All gene of one samples")
plot(density(met), xlab="Segment Mean", col="blue", lwd=3, main="All genes of all samples")

# Init output
s = data.frame(matrix(NA, length(w), 8))
colnames(s) = c("cnv.rsq", "cnv.p", "cnv.beta", "met.rsq", "met.p", "met.beta", "both.rsq", "both.p")

# Collect Regression Stats
for (i in 1: length(w)) {
	if(i %% 100 == 0) {
		cat(i, "out of", length(w), "genes done!\n")
	}
	my.expr = expr[i, ]
	my.cnv = cnv[i, ]
	my.met = met[i, ]

	var.cnv = var(my.cnv)
	var.met = var(my.met)

	if (var.cnv > 0) {	
		l = summary(lm(my.expr ~ my.cnv ))	
		s$cnv.rsq[i] = l$r.squared
		s$cnv.p[i] = l$coefficients[2,4]
		s$cnv.beta[i] = l$coefficients[2,1]
	}

	if (var.met > 0) {
		l = summary(lm(my.expr ~ my.met ))	
		s$met.rsq[i] = l$r.squared
		s$met.p[i] = l$coefficients[2,4]
		s$met.beta[i] = l$coefficients[2,1]		
	}

	if (var.cnv + var.met > 0) {
		l = summary(lm(my.expr ~ my.cnv + my.met))
		s$both.rsq[i] = l$r.squared
		f = l$fstatistic
		s$both.p[i] = 1 - pf(f[1], f[2], f[3])
	}
}
save(s, file="summary.rda")

# CNV Summary
par(mfrow=c(2,2), mgp=c(2.7,1.5,1))
rsq = s$cnv.rsq
w = which(!is.na(rsq))
rsq = rsq[w]
p = s$cnv.p[w]
beta = s$cnv.beta[w]
plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1,"Copy Number\nAlteration", cex=3)
percent = round(sum(p<0.05) / length(p) * 100)
hist(p, breaks=20, xlim=c(0,1), xlab="p-value", col="blue", lwd=3, main=paste0(percent, "% of genes significant regulated"))
percent = round(mean(rsq) * 100)
plot(density(rsq), xlim=c(0,1), xlab="R-squared", col="blue", lwd=3, main=paste0(percent, "% of gene expression explained"))
percent = round(sum(beta>=0) / length(beta) * 100)
hist(beta, xlab="Beta", col="blue", lwd=3, main=paste0(percent, "% with positive correlation"))


# MET Summary
par(mfrow=c(2,2), mgp=c(2.7,1.5,1))
rsq = s$met.rsq
w = which(!is.na(rsq))
rsq = rsq[w]
p = s$met.p[w]
beta = s$met.beta[w]
plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1,"Methylation", cex=3)
percent = round(sum(p<0.05) / length(p) * 100)
hist(p, breaks=20, xlim=c(0,1), xlab="p-value", col="blue", lwd=3, main=paste0(percent, "% of genes significant regulated"))
percent = round(mean(rsq) * 100)
plot(density(rsq), xlim=c(0,1), xlab="R-squared", col="blue", lwd=3, main=paste0(percent, "% of gene expression explained"))
percent = round(sum(beta>=0) / length(beta) * 100)
hist(beta, xlab="Beta", col="blue", lwd=3, main=paste0(percent, "% with positive correlation"))

# Both Summary
par(mfrow=c(1,2), mgp=c(2.7,1.5,1))
rsq = s$both.rsq
w = which(!is.na(rsq))
rsq = rsq[w]
p = s$both.p[w]
percent = round(sum(p<0.05) / length(p) * 100)
hist(p, breaks=20, xlim=c(0,1), xlab="p-value", col="blue", lwd=3, main=paste0(percent, "% of genes significant \nregulated by either \nCNV or Met"))
percent = round(mean(rsq) * 100)
plot(density(rsq), xlim=c(0,1), xlab="R-squared", col="blue", lwd=3, main=paste0(percent, "% of average gene \nexpression explained by \nboth CNV and Met"))












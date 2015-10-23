options(stringsAsFactors=F)
library(qvalue)

diseases = read.table("disease", h=F, stringsAsFactors=F)$V1
block.size = 1000

# Collect minimun p-values of permutation
for(disease in diseases){
print(disease)
for (g in 1:109) {
	for (r in 1:10) {
		filename = paste0(disease, "/result/", "G", g, ".R", r)
		print(filename)
		result = as.matrix(read.table(filename, h=T, row.names=1))
		print(dim(result))
		if (r == 1) {
			repeat.block = result
		} else {
			result = result[, 2:(1+block.size)]
			repeat.block = cbind(repeat.block, result)
		}
	}
	if (g == 1) {
		perm = repeat.block
	} else {
		perm = rbind(perm, repeat.block)
	}
}
print("processing...")
# Calculate imperical p-values, and q-values
colnames(perm) = c("real", paste0("p", 1:(ncol(perm)-1)))
emp.p = apply(perm, 1, function(x) rank(x)[[1]]/ncol(perm))
q.obj = qvalue(emp.p)
q = q.obj$qvalues
result = data.frame(cbind(rownames(perm), perm[,1], emp.p, q))
rownames(result) = NULL
colnames(result) = c("gene", "min.p", "emp.p", "fdr")

# Save result
save(perm, file=paste0("/scratch/zzhang4/", disease, ".perm.data.rda"))
save(q.obj, file = paste0("/scratch/zzhang4/", disease, ".qobj.rda"))
write.table(result, paste0("/scratch/zzhang4/", disease, ".fdr.txt"), col.names=T, row.names=F, sep="\t", quote=F)
}




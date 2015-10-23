args = (commandArgs(TRUE))
disease = args[1]
options(stringsAsFactors=F)
library(qvalue)

block.size = 1000

# Collect minimun p-values of permutation
for (g in 1:109) {
	for (r in 1:10) {
		filename = paste0("G", g, ".R", r)
		result = as.matrix(read.table(filename, h=T, row.names=1))
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

while read -r disease
do
	cd ./$disease/result
	R CMD BATCH --no-save --no-restore "--args disease=$disease" /scratch/zzhang4/collect.perm.r $disease.collect.perm.out &
	cd /scratch/zzhang4/
done < "disease3"



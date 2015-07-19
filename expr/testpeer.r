args = commandArgs(T)
disease = args[1]
expression_file_name = args[2]

library(MatrixEQTL)
library(qvalue)

# init file path
SNP_file_name = paste0(disease, ".SNP.txt");
snps_location_file_name = "~/peer/snpsloc.txt";
gene_location_file_name = "~/peer/geneloc.txt";
covariates_file_name = paste0("evec.txt");

useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-200;
pvOutputThreshold_tra = 0;

# Error covariance matrix
errorCovariance = numeric();

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## load expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

# Output file name
output_file_name_cis = "out.cis";
output_file_name_tra = "out.tra";

# init output
permutations <- 100
nID = ncol(gene)
nGene <- nrow(gene)
output <- matrix(nrow=nGene, ncol=permutations)

for(i in 0:permutations) {
	cat("Permutation", i, "\n")
	rm(.Random.seed)
	newseed <- as.integer(runif(1)*1000000)
	set.seed(newseed)
    
    ord = sample(nID, replace=F)
    if(i==0) ord = 1:nID
    gene = gene$ColumnSubsample(ord)

	me = Matrix_eQTL_main(
	snps=snps,
	gene=gene,
	cvrt=cvrt,
	output_file_name = output_file_name_tra,
	pvOutputThreshold = pvOutputThreshold_tra,
	useModel = useModel,
	errorCovariance = errorCovariance,
	verbose = TRUE,
	output_file_name.cis = output_file_name_cis,
	pvOutputThreshold.cis = pvOutputThreshold_cis,
	snpspos = snpspos,
	genepos = genepos,
	cisDist = cisDist,
	noFDRsaveMemory = T,
	min.pv.by.genesnp = T)

    if(i==0) {
    	real = me$cis$min.pv.gene
    } else {
    	output[,i] <- me$cis$min.pv.gene
    }

}

# output 
colnames(output) = paste0("p", 1:permutations)
output = cbind(real, output)
outfile = paste0(expression_file_name, ".perm")
write.table(output, outfile, col.names=T, row.names=T, sep="\t", quote=F)

# imperical p-value
rank = apply(output, 1, rank)
pval = (rank[1,] - 1)/permutations
qval = qvalue(pval)

# get eGene numbers
q = qval$qvalues
output = c(expression_file_name, sum(q < 0.01), sum(q < 0.05), sum(q < 0.1),sum(q < 0.2))
output = paste(output, collapse=" ")
outfile = paste0(expression_file_name, ".egene")
write.table(output, outfile, col.names=F, row.names=F, quote=F)



# TestPeer.r get disease names, detect all peer residual expression files in 
# ~/SCRATCH/eqtl/expr/peer/[disease]/, run matrixeqtl, and output the following
# expression file name, max FDR, pair comparison at FDR 0.01, pair comparison at FDR 0.1,
# number of eGene at FDR 0.01, and number of eGene at FDR 0.1
# execute as R CMD BATCH --no-save --no-restore '--args LUAD' TestPeer.r & 

args = commandArgs(T)
disease = args[1]

# init
options(stringsAsFactors=F)
library(MatrixEQTL)
setwd(paste0("~/SCRATCH/eqtl/expr/peer/", disease))

# download SNP data
cmd = paste0("s3cmd get --skip-existing s3://bioinformatics_scratch/zhenyu/SNP/", disease, ".SNP.txt")
system(cmd)

# init file path
SNP_file_name = paste0(disease, ".SNP.txt");
snps_location_file_name = "~/SCRATCH/eqtl/meta/snpsloc.txt";
gene_location_file_name = "~/SCRATCH/eqtl/meta/geneloc.txt";
covariates_file_name = paste0("~/SCRATCH/eqtl/expr/", disease, "/evec.txt");

expr.files = list.files(pattern="\\w+.+GE.PEER-\\d+.txt$")
expr.files = expr.files[which(grepl(disease, expr.files))]

useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-3;
pvOutputThreshold_tra = 0;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

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

for (i in 1: length(expr.files)) {
  ## Load gene expression data
  expression_file_name = expr.files[i]; 

  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);


  me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');

  # output
  eq = me$cis$eqtls
  w1 = which(eq$FDR < 0.01)
  w2 = which(eq$FDR < 0.1)
  n1 = length(unique(eq$gene[w1]))
  n2 = length(unique(eq$gene[w2]))
  output = c(expression_file_name, max(eq$FDR), length(w1), length(w2), n1, n2)
  output = paste(output, collapse=" ")
  outfile = paste0(expression_file_name, ".eq.summary")
  write.table(output, outfile, row.names=F, col.names=F, sep="\t", quote=F)
  # pdf("plot.pdf")
  # plot(me)
  # dev.off()
}



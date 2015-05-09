# Take array of SNP6 birdseed geneotype urls, the names you want to name each sample, and TCGA credential, and return matrix of genotypes (with quality cutoff 0.1)
GetGenotype <- function(urls, names, username, password) {
  test = read.delim(textConnection(getURL(urls[1], username=username, password=password)), skip=1)
  geno = matrix(NA, nrow(test), length(urls))
  rownames(geno) = test$Composite.Element.REF
  colnames(geno) = names

  failedIndex = NULL
  for(i in 1: length(names)) {
    data <- read.delim(textConnection(getURL(urls[i], username=username, password=password)), skip=1)
    call <- data$Call
    call[which(data$Confidence > 0.1)] <- NA
	if(length(call) == nrow(geno)) {
		geno[, i]  = call 
	} else {
		failedIndex = c(failedIndex, i)
	}
  }
  if(length(failedIndex) > 0) {
	cat("Expression data of some aliquots not retrieved.  Please try again! \n")
  }
  return(geno)
}

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
# Ask for TCGA username and password, and return a data frame with $username as username and $password as password
GetTCGACredential <- function() {
  # load package
  require(RCurl)

  valid <- F
  try <- 0

  # Try three times until succeed
  while (!valid & try < 3) {
    # Ask for username and password
    username <- readline("What is your TCGA username?")
    password <- readline("What is your TCGA password?")

    # Test username and password
    secureEntry <- paste0("https://", username, ":", password, "@tcga-data-secure.nci.nih.gov/tcgafiles/tcga4yeo/tumor/")
    valid <- grepl("Index of", getURL(secureEntry))
    try <- try + 1
  }
  
  # Feedback and return
  if(valid) {
    cat("Valid username and password\n")
  } else {
    cat("Username/password invalid: you may not access TCGA protected data\n")
    username <- ""
    password <- ""
  }   
  cred <- cbind.data.frame(username, password, stringsAsFactors=F)
  return(cred)
}

# Read TCGA Table from url or file, return data frame
# Need TCGA username and password to access protected dataset  	
GetTCGATable <- function(url, username="", password="") {
  # load package
  require(RCurl)
  
  # read table file
  if(grepl("^http", url)) {
    s <- try(getURL(url, username=username, password=password), silent = TRUE);
    if (class(s) == "try-error") 
      return(NULL);    
    tbl <- read.delim(textConnection(s), quote="\"", as.is=T)	
  } else {
    tbl = read.delim(url, quote="\"", as.is=T)
  }
  
  # remove the first one or two lines starting with "CDE_ID" or "bcr_"
  if(grepl("bcr_", tbl[1,1]) | grepl("CDE_ID", tbl[1,1]))
    tbl	<- tbl[-1,]
  if(grepl("bcr_", tbl[1,1]) | grepl("CDE_ID", tbl[1,1]))
    tbl	<- tbl[-1,]
    
  # return the data frame  
  return(tbl);
}

# read SNP6 sdrf table, output processed table with urls
ProcessSNP6Sdrf <- function(sdrf, disease) {
  # init entry point
  public.link <- paste0("https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/", disease, "/cgcc/broad.mit.edu/genome_wide_snp_6/snp/")
  protected.link <- paste0("https://tcga-data-secure.nci.nih.gov/tcgafiles/tcga4yeo/tumor/", disease, "/cgcc/broad.mit.edu/genome_wide_snp_6/snp/")
  
  # extract columns of CEL, Birdseed and Hg19.nocnv
  uuid <- sdrf$Extract.Name
  aliquot <- sdrf$Comment..TCGA.Barcode.
  file.cel <- sdrf$Array.Data.File
  url.cel <- sdrf$Comment..TCGA.Archive.Name.
  url.cel <- paste0(protected.link, url.cel, "/", file.cel)
  file.birdseed <- sdrf$Derived.Array.Data.Matrix.File.1 
  url.birdseed <- sdrf$Comment..TCGA.Archive.Name..2
  url.birdseed <- paste0(protected.link, url.birdseed, "/", file.birdseed)
  file.seg <- sdrf$Derived.Array.Data.File.3 
  url.seg <- sdrf$Comment..TCGA.Archive.Name..9
  url.seg <- paste0(public.link, url.seg, "/", file.seg)
  
  # output 
  output = cbind.data.frame(uuid, aliquot, file.cel, url.cel, file.birdseed, url.birdseed, file.seg, url.seg, stringsAsFactors=F)
}


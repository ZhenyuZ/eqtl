# Read TCGA Table from url or file, return data frame  	
GetTCGATable <- function(url) {
  # load package
  require(RCurl)
  
  # read table file
  if(grepl("^http", url)) {
    s <- try(getURL(url), silent = TRUE);
    if (class(s) == "try-error") 
      return(NULL);    
    tbl <- read.delim(textConnection(s), quote="\"", check.names=F, as.is=T)	
  } else {
    tbl = read.delim(url, quote="\"", check.names=F, as.is=T)
  }
  
  # remove the first one or two lines starting with "CDE_ID" or "bcr_"
  if(grepl("bcr_", tbl[1,1]) | grepl("CDE_ID", tbl[1,1]))
    tbl	<- tbl[-1,]
  if(grepl("bcr_", tbl[1,1]) | grepl("CDE_ID", tbl[1,1]))
    tbl	<- tbl[-1,]
    
  # return the data frame  
  return(tbl);
}


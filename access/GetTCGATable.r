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


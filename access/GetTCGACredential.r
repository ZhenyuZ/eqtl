# Ask for TCGA username and password, and return secure site url
GetTCGACredential <- function() {
  # load package
  require(RCurl)

  valid <- F
  try <- 0

  # Try three times until succeed	
  while (!valid & try < 3) {
  #	Ask for username and password
    user <- readline("What is your TCGA username?")
    passwd <- readline("What is your TCGA password?")	
    
    # Test username and password
    secureEntry <- paste0("https://", user, ":", passwd, "@tcga-data-secure.nci.nih.gov/tcgafiles/tcga4yeo/tumor/") 
    valid <- grepl("Index of", getURL(secureEntry))
    try <- try + 1
  }
  
  if(!valid) {
    secureEntry = NA
  }
  return(secureEntry)
}



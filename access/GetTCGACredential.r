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
    print("Your TCGA username and password are validated!")
    cred <- cbind.data.frame(username, password, stringsAsFactors=F)
  } else {
    print("Invalid username/password combination!")
    cred <- NULL
  }   
  return(secureEntry)
}


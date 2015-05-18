# Get TCGA Annotation JSON list of the input "disease"
# if no disease name designated, get annotation for all types
GetAnnotationJSON <- function(disease="") {
  # load packages	
  require(RCurl)
  require(RJSONIO)
  
  # get annotation file and parce to json of dccAnnotation
  # if no disease designated, get all annotation information
  if (disease == "") {
    url <- "https://tcga-data.nci.nih.gov/annotations/resources/searchannotations/json?"
    cat(paste("querying TCGA annotation database of all disease types ... ...\n"))
  } else {  
    url <- paste0("https://tcga-data.nci.nih.gov/annotations/resources/searchannotations/json?disease=", disease)
    cat(paste("querying TCGA annotation database of disease", disease, "... ...\n"))
  }  
  annot.file <- getURL(url)
  annot.json <- fromJSON(annot.file)$dccAnnotation
}


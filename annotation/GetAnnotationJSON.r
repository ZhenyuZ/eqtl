# Get TCGA Annotation JSON list of the input "disease"
GetAnnotationJSON <- function(disease) {
  # load packages	
  require(RCurl)
  require(RJSONIO)
  
  # get annotation file and parce to json of dccAnnotation
  url <- paste0("https://tcga-data.nci.nih.gov/annotations/resources/searchannotations/json?disease=", disease)
  cat(paste("querying TCGA annotation database of disease", disease, "... ...\n"))
  annot.file <- getURL(url)
  annot.json <- fromJSON(annot.file)$dccAnnotation
}



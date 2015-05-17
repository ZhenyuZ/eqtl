
# Remove TCGA barcode that exists in the annotation database
# Input: TCGA barcode
#        Annotation tbl
#        array of categoryId to keep.  To remove all annotation barcode, use "none". Default is c("6", "204"), corresponding to prior and synchronous malignancy
#        recursive flag. Default is TRUE, to remove not only the exact barcode, but all derived barcode
# Output: Filtered TCGA barcode
FilterByAnnotation <- function(barcode, annot.tbl, keep=c("6", "204"), recursive=T) {  
  # remove annotation with categoryId in character array "keep"
  barcodeToRemove <- unique(annot.tbl$item[! annot.tbl$categoryId %in% keep])
  
  if (recursive) {
    # recursively find annotation barcode, as well as any derived barcode, from the query barcode
    regexToRemove <- paste0("(", paste(barcodeToRemove, collapse="|"), ")")
    indexToRemove <- which(grepl(regexToRemove, barcode))
  } else { 
    # only find annotation barcode of exact mach
    indexToRemove <- which(barcode %in% barcodeToRemove)
  } 
   
  # remove the barcode detected 
  if(length(indexToRemove) > 0) {
    barcode = barcode[-w]
  }
  return(barcode)
}

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


# Input Annotation JSON list, output Annotation Table with 3 columns "id", "item" and "categoryId".  If a outFile name is given, the function will save the table into a RData file.
GetAnnotationTable <- function(annot.json, outFile=NA){
  # parse annotation json
  id <- as.character(sapply(annot.json, function(x) x$id))
  item <- as.character(sapply(annot.json, function(x) x$items[[1]]$item))
  categoryId <- as.character(sapply(annot.json, function(x) x$annotationCategory$categoryId))
  annot.tbl <- cbind.data.frame(id, item, categoryId, stringsAsFactors=F)
  
  # save and return annotation table
  if(!is.na(outFile)) {
    save(annot.tbl, file=outFile)
  }
  return(annot.tbl)
}
# Input vector of TCGA barcodes and an Annotation Table, output Annotation List with barcode as names, and combined categoryId as values
QueryBarcodeAnnotation <- function(barcode, annot.tbl){
	annot.lst <- sapply(barcode, function(y) paste(sort(unique(annot.tbl$categoryId[which(sapply(annot.tbl$item, function(x) grepl(x, y)))])), collapse="|") )
	names(annot.lst) <- barcode
	return(annot.lst)
}


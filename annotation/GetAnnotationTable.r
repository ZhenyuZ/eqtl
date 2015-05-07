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

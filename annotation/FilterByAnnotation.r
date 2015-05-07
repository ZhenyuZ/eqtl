
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


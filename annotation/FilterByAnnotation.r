# Remove TCGA barcode that exists in the annotation database
# Input: TCGA barcode
#        Annotation tbl
#        array of categoryId to keep.  To remove all annotation barcode, use "none". Default is c("6", "204"), corresponding to prior and synchronous malignancy
#        recursive flag. Default is TRUE, to remove not only the exact barcode, but all derived barcode
# Output: Filtered TCGA barcode
FilterByAnnotation <- function(barcodes, annot.tbl, keep=c("6", "204"), recursive=T) {
  # remove annotation with categoryId in character array "keep"
  barcodeToRemove <- unique(annot.tbl$item[! annot.tbl$categoryId %in% keep])

  if (recursive) {
    # recursively find annotation barcode, as well as any derived barcode, from the query barcode
    isFound <- sapply(barcodes, function(x) 
                       (substr(x, 1, 12) %in% barcodeToRemove) | 
                       (substr(x, 1, 15) %in% barcodeToRemove) | 
                       (substr(x, 1, 16) %in% barcodeToRemove) | 
                       (substr(x, 1, 19) %in% barcodeToRemove) | 
                       (substr(x, 1, 20) %in% barcodeToRemove) | 
                       (substr(x, 1, 23) %in% barcodeToRemove) | 
                       (substr(x, 1, 24) %in% barcodeToRemove) | 
                       (substr(x, 1, 27) %in% barcodeToRemove) | 
                       (substr(x, 1, 28) %in% barcodeToRemove) )
    barcodes <- barcodes[which(!isFound)]                  
  } else {
    # only find annotation barcode of exact mach
    barcodes <- barcodes[which(!barcodes %in% barcodeToRemove)]
  }
  return(barcodes)
}

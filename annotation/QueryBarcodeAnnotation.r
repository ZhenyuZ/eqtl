# Input vector of TCGA barcodes and an Annotation Table, output Annotation List with barcode as names, and combined categoryId as values
QueryBarcodeAnnotation <- function(barcode, annot.tbl){
	annot.lst <- sapply(barcode, function(y) paste(sort(unique(annot.tbl$categoryId[which(sapply(annot.tbl$item, function(x) grepl(x, y)))])), collapse="|") )
	names(annot.lst) <- barcode
	return(annot.lst)
}


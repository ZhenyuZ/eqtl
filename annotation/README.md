# TCGA Annotation Database Query
Functions to query and filter TCGA annotation database

## GetAnnotationJSON.r
Query Annotation Database with the disease type, such as "LUAD" or "OV", return json format of the database values <br>
`GetAnnotationJSON(disease)`

## GetAnnotationTable.r
Convert json format into character data frame of three columns "id", "item" and "categoryId".  
If outFile RData file is indicated, the result will also be saved <br> 
`GetAnnotationTable <- function(annot.json, outFile=NA)`

## QueryBarcodeAnnotation.r
Input vector of TCGA barcodes and an Annotation Table, output Annotation List with barcode as names, 
and combined categoryId (separated by "|") as values <br>
`QueryBarcodeAnnotation(barcode, annot.tbl)`

## FilterByAnnotation.r
Remove TCGA barcode that exists in the annotation database <br>
Input:  <br>
        * TCGA barcode <br>
        * Annotation tbl <br>
        * array of categoryId to keep.  To remove all annotation barcode, use "none". Default is c("6", "204"), 
          corresponding to prior and synchronous malignancy <br>
        * recursive flag. Default is TRUE, to remove not only the exact barcode, but all derived barcode <br>
Output: Filtered TCGA barcode <br>
`FilterByAnnotation(function(barcode, annot.tbl, keep=c("6", "204"), recursive=T))`

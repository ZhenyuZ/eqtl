# TCGA Assember
This directory contains key files of the [TCGA Assembler](http://www.compgenome.org/TCGA-Assembler/) <br>
[TCGA Assembler](http://www.compgenome.org/TCGA-Assembler/) is developed by Yitan Zhu, Yuan Ji, et al.  <br>
Please cite: <br>
`Y. Zhu, P. Qiu, Y. Ji*. TCGA-Assembler: Open-Source Software for Retrieving and Processing TCGA Data. Nature Methods. 11:599-600, 2014. | doi:10.1038/nmeth.2956`

## Module_C
Module_C is a nasty modification of Module_A and Module_B, and the intent is to get TCGA sample Barcode without download the real data <br>
Module_C contains the following functions <br>
### GetSampleID
Get Sample IDs of the designed tissueCode from SDRF.  Not working currently due to recent SDRF changes. <br>
`GetSampleID(sdrf, tissueCode = "")`

### GetGenotypeDataSample
Get Sample IDs of the designed tissueCode from SDRF file.  Not working currently due to recent SDRF changes. <br>
`GetGenotypeDataSample(sdrf, tissueCode="10") `

### GetCNADataSample
Get SNP6 CNV Sample IDs. Not tested after recent SDRF changes. <br>
`GetCNADataSample(traverseResultFile = travFile, saveFolderName = type, cancerType = type, assayPlatform = "genome_wide_snp_6", inputPatientIDs = NULL, tissueCode = "")`

### GetMethylationDataSample 
Get Methylation Array Sample IDs. Not tested after recent SDRF changes. <br>
`GetMethylationDataSample(traverseResultFile = travFile, saveFolderName = type, cancerType = type, assayPlatform, inputPatientIDs = NULL, tissueCode = "")`

### GetRNASeqDataSample
Get RNA-Seq Sample IDs. Not tested after recent SDRF changes. <br>
`GetRNASeqDataSample(traverseResultFile = travFile, saveFolderName = type, cancerType = type, assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results", inputPatientIDs = NULL, tissueCode = "")`

### filterID
Return TCGA Sample IDs that of a list of Patient IDs (4 character code) in "filter" <br>
`filterID (sample, filter)`



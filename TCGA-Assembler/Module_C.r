
# Convert sdrf variables to array of sample ids, and extract only matched tissue
GetSampleID <- function(sdrf, tissueCode = "")
{
  samples = NULL
  ids = sdrf[,1]
  for (i in 1: length(ids))
  {
  	id = strsplit(ids[i],"-")[[1]]
	tissue = substr(id[4],1,2) 
	if(tissueCode == "" || tissueCode == tissue) 
	{
		samples = append(samples, paste(id[1],id[2],id[3],tissue, sep="-"))
	}
  }
  return(samples)
}  
   
GetGenotypeDataSample <- function(sdrf, tissueCode="10") 
{
  snp6.sdrf = read.table(sdrf, h=T, sep="\t")
  snp6.aliquot = snp6.sdrf$Comment..TCGA.Barcode.
  snp6.NB.aliquot = snp6.aliquot[which(sapply(snp6.aliquot, function(x) substr(x, 14,15)) == tissueCode)]
  snp6.sample = unique(substr(snp6.NB.aliquot, 1, 15))
  return(snp6.sample)
}   
   
GetCNADataSample <- function(traverseResultFile = travFile, saveFolderName = type, cancerType = type, assayPlatform = "genome_wide_snp_6", inputPatientIDs = NULL, tissueCode = "")
{
  options(warn=-1);
  
  writeLines("**********************************************************************************");
  writeLines("");
  writeLines(paste("Getting information about copy number data of ", cancerType, " patients.", sep = ""));
  
  if (!is.null(inputPatientIDs))
  {
    inputPatientIDs = toupper(gsub(pattern = "\\.", replacement = "-", x = inputPatientIDs));
  }
  
  # load directory traverse result and create folder to store output files
  writeLines("Load information of TCGA data files.");
  # dir.create(path = saveFolderName, recursive = TRUE);
  SpecificID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/broad\\.mit\\.edu/", assayPlatform, "/", sep = "")), x = upper_file_url, ignore.case = FALSE);
  IDLevel_3InFile_Url = SpecificID[grep(pattern = toupper("Level_3"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
  SpecificName = paste(cancerType, "__broad.mit.edu__", assayPlatform, sep = "");
  
  # download and process Sample and Data Relationship Format (SDRF) file
  ind = SpecificID[grepEnd(pattern = toupper("\\.sdrf\\.txt"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
  if (length(ind) == 0)
  {
    writeLines("Program existed due to missing SDRF file."); 
    return();
  }
  if (length(ind) > 1)
  {
    URL = GetNewestURL(AllURL = file_url[ind]);
  }else{
    URL = file_url[ind];
  }
  downloadResult = urlReadTable(url = URL);
  if (downloadResult$errorFlag != 0)
  {
    writeLines("Error in downloading SDRF file.");
    return();
  }

  sdrf = toupper(downloadResult$data);
  level_3_filename_column = which(sdrf[1, ] == toupper("Derived Array Data File"));
  DataLevelColID = which(sdrf[1, ] == toupper("Comment [TCGA Data Level]"));
  DataLevelColID = DataLevelColID[(length(DataLevelColID)-length(level_3_filename_column)+1):length(DataLevelColID)];
  BarcodeColID = which(sdrf[1, ] == toupper("Comment [TCGA Barcode]"));
  colnames(sdrf) = sdrf[1, ];
  sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
  Level3_ID = 1:dim(sdrf)[1];
  for (DataLevelColIDIndex in 1:length(DataLevelColID))
  {
    Level3_ID = sort(intersect(Level3_ID, union(which(sdrf[, DataLevelColID[DataLevelColIDIndex]] == "LEVEL_3"),
                    which(sdrf[, DataLevelColID[DataLevelColIDIndex]] == "LEVEL 3"))), decreasing = FALSE);
  }
  if (length(Level3_ID) == 0)
  {
    writeLines("Error: there are no Level 3 data");
    return();
  }
  sdrf = unique(sdrf[Level3_ID, sort(c(BarcodeColID, level_3_filename_column, DataLevelColID), decreasing = FALSE), drop = FALSE]);  
  
  # If specific patient TCGA barcodes are inputted, only download data of the specified samples.
  if (!is.null(inputPatientIDs))
  {
    indInputPatientID = c();
    for (i in 1:length(inputPatientIDs))
    {
      indInputPatientID = c(indInputPatientID, grepBeginning(pattern = toupper(inputPatientIDs[i]), x = sdrf[, 1], ignore.case = FALSE));
    }
    if (length(indInputPatientID) == 0)
    {
      writeLines("No Level 3 data for the inputted TCGA barcodes.");
      return();      
    }else{
      sdrf = sdrf[indInputPatientID, , drop = FALSE];
    }
  }
  samples = GetSampleID(sdrf, tissueCode = tissueCode)
  return(samples)
}


GetMethylationDataSample <- function(traverseResultFile = travFile, saveFolderName = type, cancerType = type, assayPlatform, inputPatientIDs = NULL, tissueCode = "")
{
  options(warn=-1);
  
  writeLines("**********************************************************************************");
  writeLines("");
  writeLines(paste("Getting information about DNA methylation data of ", cancerType, " patients.", sep = ""));
  
  # Check whether specific TCGA patient IDs are inputted. 
  if (!is.null(inputPatientIDs))
  {
    inputPatientIDs = toupper(gsub(pattern = "\\.", replacement = "-", x = inputPatientIDs));
  }
  

  # load directory traverse result and create folder to store output files
  writeLines("Load information of TCGA data files.");
  # dir.create(path = saveFolderName, recursive = TRUE);
  SpecificID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/jhu-usc\\.edu/", assayPlatform, "/", sep = "")), x = upper_file_url, ignore.case = FALSE);
  MeLevel3ID = SpecificID[grep(pattern = toupper("Level_3"), x = upper_file_url[SpecificID], ignore.case = FALSE)];

  # search for the Sample and Data Relationship Format (SDRF) file of the specified platform and cancer type
  ind = SpecificID[grepEnd(pattern = toupper("\\.sdrf\\.txt"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
  if (length(ind) == 0)
  {
    writeLines("Program existed due to missing SDRF file."); 
    return();
  }
  if (length(ind) > 1)
  {
    URL = GetNewestURL(AllURL = file_url[ind]);
  }else{
    URL = file_url[ind];
  }
  downloadResult = urlReadTable(url = URL);
  if (downloadResult$errorFlag != 0)
  {
    writeLines("Error in downloading SDRF file.");
    return();
  }
  sdrf = toupper(downloadResult$data);
  
  # Process SDRF file, identify the columns of level 3 data file name and TCGA sample barcodes.
  level_3_filename_column = max(grep(pattern = "Data Matrix File", x = sdrf[1, ], ignore.case = TRUE));
  DataLevelColID = max(grep(pattern = "TCGA Data Level", x = sdrf[1, ], ignore.case = TRUE));
  TCGABarcodeID = min(grep(pattern = "TCGA Barcode", x = sdrf[1, ], ignore.case = TRUE));
  colnames(sdrf) = sdrf[1, ];
  sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
  sdrf = sdrf[!duplicated(sdrf[, level_3_filename_column]), c(TCGABarcodeID, DataLevelColID, level_3_filename_column), drop = FALSE];
  SDRFID = sort(union(which(sdrf[, 2] == "LEVEL_3"), which(sdrf[, 2] == "LEVEL 3")), decreasing = FALSE);
  if (length(SDRFID) == 0)
  {
    writeLines("Error: there are no Level 3 data");
    return();
  }
  sdrf = sdrf[SDRFID, , drop = FALSE];  
  
  # If specific patient TCGA barcodes are inputted, only download the specified samples.
  if (!is.null(inputPatientIDs))
  {
    indInputPatientID = c();
    for (i in 1:length(inputPatientIDs))
    {
      indInputPatientID = c(indInputPatientID, grepBeginning(pattern = inputPatientIDs[i], x = sdrf[, 1], ignore.case = FALSE));
    }
    if (length(indInputPatientID) == 0)
    {
      writeLines("No Level 3 data for the inputted TCGA barcodes.");
      return();      
    }else{
      sdrf = sdrf[indInputPatientID, , drop = FALSE];
    }
  }
  samples = GetSampleID(sdrf, tissueCode = tissueCode)
  return(samples)
}
  


GetRNASeqDataSample <-function (traverseResultFile = travFile, saveFolderName = type, cancerType = type, assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results", inputPatientIDs = NULL, tissueCode = "")
{
  options(warn=-1);
  
  writeLines("**********************************************************************************");
  writeLines("");
  writeLines(paste("Getting information about RNA-seq data of ", cancerType, " patients.", sep = ""));
  
  # Check whether specific TCGA patient IDs are inputted. 
  if (!is.null(inputPatientIDs))
  {
    inputPatientIDs = toupper(gsub(pattern = "\\.", replacement = "-", x = inputPatientIDs));
  }

  # load directory traverse result and create folder to store output files
  writeLines("Load information of TCGA data files.");
  # dir.create(path = saveFolderName, recursive = TRUE);
  if (assayPlatform == "RNASeqV1")
  {
    platform = c("illuminaga_rnaseq", "illuminahiseq_rnaseq"); 
    Institution = c("unc.edu", "bcgsc.ca");    
  }
  if (assayPlatform == "RNASeqV2")
  {
    platform = c("illuminaga_rnaseqv2", "illuminahiseq_rnaseqv2"); 
    Institution = c("unc.edu");    
  }
  
  # For returning downloaded data
  downloadedData = vector("list", 0);     
  dataIndex = 0;     
  
  # download RNASeqV2 data
  if (assayPlatform == "RNASeqV2")
  {
    for (IDin in 1:length(Institution))
    {
      for (IDpl in 1:length(platform))
      {
        SpecificName = paste(cancerType, "__", Institution[IDin], "__", platform[IDpl], sep = "");
        SpecificID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/", Institution[IDin], "/", platform[IDpl], "/", sep = "")), x = upper_file_url, ignore.case = FALSE);
        RNALevel3ID = SpecificID[grep(pattern = toupper("Level_3"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        ind = SpecificID[grepEnd(pattern = toupper("sdrf\\.txt"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        if (length(ind) == 0)
        {
          next;
        }
        if (length(ind) > 1)
        {
          URL = GetNewestURL(AllURL = file_url[ind]);
        }else{
          URL = file_url[ind];
        }
        downloadResult = urlReadTable(url = URL);
        if (downloadResult$errorFlag != 0)
        {
          next;
        }
        sdrf = toupper(downloadResult$data);
        
        level_3_filename_column = max(grep(pattern = "Derived Data File", x = sdrf[1, ], ignore.case = TRUE));
        DataLevelColID = max(grep(pattern = "Comment \\[TCGA Data Level]", x = sdrf[1, ], ignore.case = TRUE));
        ExtractNameColID = grep(pattern = "Comment \\[TCGA Barcode]", x = sdrf[1, ], ignore.case = TRUE);
        RefGenomeColID = grep(pattern = "Comment \\[Genome reference]", x = sdrf[1, ], ignore.case = TRUE);
        if (length(ExtractNameColID) == 0)
        {
          ExtractNameColID = min(grep(pattern = "Extract Name", x = sdrf[1, ], ignore.case = TRUE));          
        }
        colnames(sdrf) = sdrf[1, ];
        sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
        Level3_ID = sort(union(which(sdrf[, DataLevelColID] == "LEVEL_3"), which(sdrf[, DataLevelColID] == "LEVEL 3")), decreasing = FALSE);
        if (length(Level3_ID) == 0)
        {
          next;
        }
        sdrf = sdrf[Level3_ID, c(ExtractNameColID, level_3_filename_column, RefGenomeColID), drop = FALSE];
        sdrf = sdrf[!duplicated(sdrf[, 2]), , drop = FALSE];
        
        # Only keep the file information for the data types that should be downloaded.
        keepID = c();
        for (keep_i in 1:length(dataType))
        {
          keepID = c(keepID, grep(pattern = dataType[keep_i], x = sdrf[, 2], ignore.case = TRUE))
        }
        sdrf = sdrf[sort(unique(keepID), decreasing = FALSE), , drop = FALSE];
        
        # If specific patient TCGA barcodes are inputted, only download the specified samples.
        if (!is.null(inputPatientIDs))
        {
          indInputPatientID = c();
          for (i in 1:length(inputPatientIDs))
          {
            indInputPatientID = c(indInputPatientID, grepBeginning(pattern = toupper(inputPatientIDs[i]), x = sdrf[, 1], ignore.case = FALSE));
          }
          if (length(indInputPatientID) == 0)
          {
            next;
          }else{
            sdrf = sdrf[indInputPatientID, , drop = FALSE];
          }
        }
      }  
    }
  } 
  samples = GetSampleID(sdrf, tissueCode = tissueCode)
  return(samples)
}  	  

# filter TCGA barcode by list of patient IDs in filter
filterID <- function (sample, filter)
{
	if(length(sample) > 0 )
	{
		ids = unlist(strsplit(sample,"-"))[4*(1:length(sample))-1]
		sample = sample[which(ids %in% filter)]
	}
	return(sample)
}


  

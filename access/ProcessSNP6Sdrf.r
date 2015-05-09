# read SNP6 sdrf table, output processed table with urls
ProcessSNP6Sdrf <- function(sdrf, disease) {
  # init entry point
  public.link <- paste0("https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/", disease, "/cgcc/broad.mit.edu/genome_wide_snp_6/snp/")
  protected.link <- paste0("https://tcga-data-secure.nci.nih.gov/tcgafiles/tcga4yeo/tumor/", disease, "/cgcc/broad.mit.edu/genome_wide_snp_6/snp/")
  
  # extract columns of CEL, Birdseed and Hg19.nocnv
  uuid <- sdrf$Extract.Name
  aliquot <- sdrf$Comment..TCGA.Barcode.
  file.cel <- sdrf$Array.Data.File
  url.cel <- sdrf$Comment..TCGA.Archive.Name.
  url.cel <- paste0(protected.link, url.cel, "/", file.cel)
  file.birdseed <- sdrf$Derived.Array.Data.Matrix.File.1 
  url.birdseed <- sdrf$Comment..TCGA.Archive.Name..2
  url.birdseed <- paste0(protected.link, url.birdseed, "/", file.birdseed)
  file.seg <- sdrf$Derived.Array.Data.File.3 
  url.seg <- sdrf$Comment..TCGA.Archive.Name..9
  url.seg <- paste0(public.link, url.seg, "/", file.seg)
  
  # output 
  output = cbind.data.frame(uuid, aliquot, file.cel, url.cel, file.birdseed, url.birdseed, file.seg, url.seg, stringsAsFactors=F)
}


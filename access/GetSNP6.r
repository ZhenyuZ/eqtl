# Scripts to download, extract and harmonized TCGA SNP6 birdseed genotypes
# need to run line by line with modifications

# load packages
options(stringsAsFactors=F)
source("/Users/zhenyu/github/eqtl/module.access.r")

# provide sdrf that contain file and aliquot information
sdrf.link <- "https://tcga-data-secure.nci.nih.gov/tcgafiles/tcga4yeo/tumor/luad/cgcc/broad.mit.edu/genome_wide_snp_6/snp/broad.mit.edu_LUAD.Genome_Wide_SNP_6.mage-tab.1.2012.0/broad.mit.edu_LUAD.Genome_Wide_SNP_6.sdrf.txt"

# get username and password for TCGA protected data
cred <- GetTCGACredential()

# read sdrf as table
sdrf <- GetTCGATable(sdrf.link, cred$username, cred$password)

# extract information from sdrf
file.info <- ProcessSNP6Sdrf(sdrf, "luad")

# download all birdseed data provided by the urls, and combine to a matrix file
geno <- GetGenotype(file.info$url.birdseed, file.info$file.birdseed, cred$username, cred$password)

# save the birdseed matrix
save(geno, file="luad.birdseed.rda")

# replace colnames with TCGA aliquot barcode, and rownames with dbSNP names
colnames(geno) = file.info$aliquot
conv <- read.table("/glusterfs/netapp/homes1/ZHANGZ18/meta/SNP_A2dbsnp.txt", as.is=T)
rownames(geno) <- conv$V2[match(rownames(geno), conv$V1)]

# remove SNPs without "rs" names
geno <- geno[which(grepl("^rs", rownames(geno))), ]

# save the processed dbSNP matrix
save(geno, file="luad.geno.rda")


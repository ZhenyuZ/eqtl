# SplitDosage.r filter dosage file by plink-filtered snps, and split 
# dosage files by disease types

options(stringsAsFactors=F)
library(data.table)
diseases = unlist(read.table("~/SCRATCH/plink2/diseases.txt"))

# read meta data
info = fread("~/SCRATCH/plink2/tcga.snp.info", h=T)
tcga.fam = fread("~/SCRATCH/plink2/tcga.geno.fam", h=F)
ids = tcga.fam$V2

for (i in 1:22) {
  # read dosage
  dosage.file = paste0("~/SCRATCH/impute2/dosage.chr/chr", i, ".dosage")
  dosage = fread(dosage.file, h=F)

  # filter by existence in info 
  snp = dosage$V1
  m = match(snp, info$impute2.snp)
  snp = info$snp[m]
  w = which(!is.na(snp))
  snp = snp[w]

  # make dosage as matrix
  dosage = as.matrix(dosage[w, 2:ncol(dosage), with=F])
  rownames(dosage) = snp
  colnames(dosage) = ids

  # save dosage by chromosome r data file
  outfile = paste0("~/SCRATCH/impute2/trimed.dosage/chr", i, ".dosage.rda")
  save(dosage, file=outfile)

  # read disease bim and fam, and split dosage by disease type
  for (disease in diseases) {
    disease = "UCEC"
    
    # read disease bim and fam
    bim = fread(paste0("~/SCRATCH/plink2/", disease, ".bim"), h=F)
    setnames(bim, c("chr", "snp", "dist", "pos", "ref", "alt"))
    fam = fread(paste0("~/SCRATCH/plink2/", disease, ".fam"), h=F)
    setnames(fam, c("aliquot", "id", "V3", "V4", "V5", "V6"))

    # trim by id
    m = match(fam$id, colnames(dosage))
    my.dosage = dosage[, m]
    
    # trim by snp
    w = which(rownames(dosage) %in% bim$snp)
    my.dosage = my.dosage[w, ]
    
    # output
    outfile = paste0("~/SCRATCH/SNP/", disease, ".header")
    header = paste(c("id", colnames(my.dosage)), collapse="\t")
    write.table(header, outfile, col.names=F, row.names=F, quote=F, sep="\t")
    outfile = paste0("~/SCRATCH/SNP/", disease, ".chr", i, ".SNP.txt")
    write.table(dosage, outfile, col.names=F, row.names=T, quote=F, sep="\t")
  }
}



    

m = match(bim$snp, info$impute2.snp)

while read -r line; do awk '{printf $2; for(i=6; i<NF; i+=3) {printf "\t"$(i+0)*0+$(i+1)*1+$(i+2)*2}; printf "\n"}' $line/$line.impute2 > dosage/$line.dosage; done < "chunks" 


      
  
  while read -r line
  do 
      echo $line
  done < $file
  
  
while read -r line; 
do


# change directory
dir = paste0("/home/ubuntu/SCRATCH/eqtl/expr/", disease)
setwd(dir)
  
# load expr, met, cnv, expr.aliquot and patient data
load("expr.residuals.rda")
  
# read clinical info
clin.file = paste0("~/SCRATCH/eqtl/clin/nationwidechildrens.org_clinical_patient_", tolower(disease), ".txt") 
clin = read.table(clin.file, h=T, sep="\t", quote="")
  
# manually pick potential clinical covariates
keep = c("bcr_patient_barcode", "gender", "ajcc_pathologic_tumor_stage", "ajcc_clinical_tumor_stage", "age_at_initial_pathologic_diagnosis", "clinical_stage", "histologic_subtype", "age_at_diagnosis", "menopause_status", "histological_type", "tumor_grade", "tobacco_smoking_history_indicator")


# Stop here, check clinical variable names
print(colnames(clin))
optional = c("tumor_focality", "braf_mutation_result", "ras_mutation_result")
  
# only get clinical variable in keep and optional, and check their values   
# and also match with the patients in expr 
clin2 = clin[, which(colnames(clin) %in% c(keep, optional)) ]
clin2 = clin2[-c(1,2), ]
m = match(substr(colnames(expr), 1, 12), clin2$bcr_patient_barcode)
clin2 = clin2[m, ]


# display clin2
  for(i in 1: ncol(clin2)) { 
    if (colnames(clin2)[i] != "bcr_patient_barcode"){
      cat(paste(i, disease, colnames(clin2)[i]), "\n")
      print(table(clin2[,i]))
      cat("\n\n")
    }  
  }
  
# further remove variables by values  
clin3=clin2[, c(1,2,4)]

# display clin3
  for(i in 1: ncol(clin3)) { 
    if (colnames(clin3)[i] != "bcr_patient_barcode"){
      cat(paste(i, disease, colnames(clin3)[i]), "\n")
      print(table(clin3[,i]))
      cat("\n\n")
    }  
  }
  
# code each variable  
i=1
  
i=i+1
temp = clin3[,i]
table(temp)

# code tumor_grade
temp[which(temp=="GB")] = 0
temp[which(temp=="G1")] = 1
temp[which(temp=="G2")] = 2
temp[which(temp=="G3")] = 3
temp[which(temp=="G4")] = 4
temp = imp(temp)
clin3[,i] = temp

# code testis cancer serum marker
temp[which(temp=="S0")] = 0
temp[which(temp=="S1")] = 1
temp[which(temp=="S2")] = 2
temp[which(temp=="S3")] = 3
temp[which(temp=="S4")] = 4
temp = imp(temp)
clin3[,i] = temp


# code for sarcoma subtype
temp[which(temp == "Conventional leiomyosarcoma")] = 2
temp[which(temp == "Poorly differentiated or pleomorphic or epithelioid leiomyosarcoma")] = 1
temp[which(temp == "Well-differentiated leiomyosarcoma (resembling leiomyoma)")] = 3
temp = imp(temp)
clin3[,i] = temp

# code level variables
temp[which(temp=="Elevated")] = 1
temp[which(temp=="Normal")] = 0
temp[which(temp=="Low")] = -1
temp = imp(temp)
clin3[,i] = temp

# code Yes/No variables
temp[which(grepl("YES", temp))] = 1
temp[which(grepl("NO", temp))] = 0
temp = imp(temp)
clin3[,i] = temp

# code gender    
temp[which(grepl("FEMALE", temp))] = 2
temp[which(grepl("MALE", temp))] = 1
clin3[,i] = temp

# code manopause
temp[which(grepl("^Pre", temp))] = -1
temp[which(grepl("^Post", temp))] = 1
temp[which(grepl("^Inde", temp))] = 0
temp[which(grepl("^Peri", temp))] = 0
clin3[,i] = temp

# code positive/negative variable
temp[which(grepl("Negative", temp))] = -1
temp[which(grepl("Positive", temp))] = 1
temp[which(grepl("Indeterminate", temp))] = 0
temp = imp(temp)
clin3[,i] = temp

# code smoking indicator
temp[which(grepl("Lifelong Non-smoker", temp))] = 0
temp[which(grepl("Current smoker", temp))] = 3
temp[which(grepl("Current reformed smoker for > 15 years", temp))] = 1
temp[which(grepl("Current reformed smoker for < or = 15 years", temp))] = 2
w = which(temp %in% c(1,2))
m = mean(as.numeric(temp[w]))
temp[which(grepl("Current Reformed Smoker, Duration Not Specified", temp))] = m
temp = imp(temp)
clin3[,i] = temp

# code others
temp[which(grepl("Non-Papillary", temp))] = 0
temp[which(grepl("Papillary", temp))] = 1
temp = imp(temp)
colnames(clin3)[i] = "Papillary"
clin3[,i] = temp
  
# transform tumor_grade 
getGrade = function(grade) {
  grade = gsub("Stage ", "", grade)
  out = grade
  out[which(grepl("^I", grade))] = 1
  out[which(grepl("^II", grade))] = 2
  out[which(grepl("^III", grade))] = 3
  out[which(grepl("^IV", grade))] = 4
  out[which(grepl("^V", grade))] = 5
  
  w = which(! out %in% c(1,2,3,4,5))
  print(table(grade[-w])) 
  return(out)
}      

# imputation missing values by average
imp = function(data) {
  data = as.matrix(data)
  mode(data) = "numeric"
  for(i in 1: ncol(data)) {
    temp = data[,i]
    w = which(is.na(temp))
    temp[w] = mean(temp, na.rm=T)
    data[, i] = temp
  }
  return(data)
}    

# qqplot to examine the effect of each clinical variable on expression
clin = clin3
for(i in 2: ncol(clin)) {
  cov = as.numeric(clin[,i])
  name = colnames(clin)[i]
  p = apply(expr, 1, function(x) cor.test(x, cov, method="pearson")$p.value)
  out.file=paste0("qqplot.", i, ".png")
  png(out.file, width=960, height=960)
  qq(p, main=name)
  dev.off()
}


      

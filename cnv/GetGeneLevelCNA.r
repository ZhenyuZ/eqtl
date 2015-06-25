# GetGeneLevelCNA.r
# Input TCGA segmentation data (long format of multiple aliquots)
# Calculate Gene Level copy number based on the following steps:
# 1. for each gene, calculated overlapped segments
# 2. regenerated copy number value from segment mean, and calculate 
# average weighted on the range of segment overlaps (assume 2 for
# any region there is no segment overlaps
# output matrix format of gene x patient copy number values
#
# This code is quite slow, mainly b/c it need to loop through all 
# segments for each gene in each patient.  Potentially, it could 
# be speed up by collapse the aggregate copy number calculation 
# without using loop.  

# load library and init
options(stringsAsFactors=F)
library("GenomicRanges")

# read TCGA GAF2 location file to define gene locations
locfile="./meta/geneloc.for.cnv.txt"
geneloc=read.table(locfile, h=T, sep="\t", colClasses="character")
loc = GRanges(seqnames = geneloc[,2],
               ranges = IRanges(start=as.numeric(geneloc[,3]), end=as.numeric(geneloc[,4])),
               strand = "+",
               name = geneloc[,1])
gene.length = width(loc)               
n = length(loc)

# read diseases
diseases = unlist(read.table("diseases.txt"))

# read patient aliquot mapping
euro.file = "./aliquot/euro.aliquot.summary.txt"
euro = data.frame(read.table(euro.file, h=T, comment.char=""))

# Calculate gene level CNV for each disease type
for (disease in diseases){

  # read long format combined segment file for one disease type
  segfile = paste0("./cnv/", disease, ".snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
  all.cnvseg = data.frame(read.table(segfile, h=T, colClasses="character"))

  # filter by patient and aliquot ID
  patient.file = paste0("./plink/geno/", disease, ".patient")
  patients = unlist(read.table(patient.file, h=F))
  aliquots = euro$CNA[match(patients, euro$patient)]

  # sanity check to see if all aliquots in the summary file are in our current dataset
  all.aliquots = unique(all.cnvseg$Sample)
  if(sum(!aliquots %in% all.aliquots)> 0) {
    print ("WARNINGS! missing CNA aliquot")
  }
  
  # init output
  output = matrix(0, n, length(aliquots)) 
  colnames(output) = patients
  
  # Calculate gene level CNV for each patient
  for (patient.index in 1:length(patients)) {
    
    # extract segment information for one patient
    aliquot = aliquots[patient.index]
    cnvseg = all.cnvseg[which(all.cnvseg$Sample==aliquot), ]  
  
    # load into grange object
    seg = GRanges(seqnames = paste("chr",cnvseg[,2],sep=""),
                  ranges = IRanges(start=as.numeric(cnvseg[,3]), end=as.numeric(cnvseg[,4])),
               strand = "+",
                Num_Probes = as.numeric(cnvseg[,5]),
               Segment_Mean = as.numeric(cnvseg[,6]) )
    
    segm = values(seg)[,"Segment_Mean"]

    # detect range overlap between gene locations and segment ranges
    Hits = findOverlaps(loc,seg)
    hit = data.frame(cbind(queryHits(Hits), subjectHits(Hits)))
    colnames(hit) = c("loc.index", "seg.index")
    
    # extract overlap range borders, width, and segment mean
    myloc = loc[hit$loc.index]
    myseg = seg[hit$seg.index]
    hit$loc.start = start(myloc)
    hit$loc.end = end(myloc)
    hit$seg.start = start(myseg)
    hit$seg.end = end(myseg)
    hit$overlap = apply(hit, 1, function(x){y=x[3:6]; o = order(y); y[o[3]] - y[o[2]] + 1})
    hit$segm = segm[hit$seg.index]

    # init cnv matrix to accept output, first column to be aggregated copy 
    # number in the gene location; second column to be aggregated length of 
    # overlapped ranges (which finally should reach the length of that gene
    cnv=data.frame(matrix(0, length(loc), 2))
    colnames(cnv) = c("aggregate.cn", "aggregate.length")

    # loop through the overlap detection data.frame, and aggregate copy number
    # and overlapped length.  
    # by definition: segmean = log2(cn/2)
    # recover copy number: cn = 2^segmean * 2
    for(i in 1: nrow(hit)){
      loc.index = hit$loc.index[i]
      seg.index = hit$seg.index[i]
      cnv$aggregate.cn[loc.index] = cnv$aggregate.cn[loc.index] + (2^hit$segm[i]) * hit$overlap[i] * 2
      cnv$aggregate.length[loc.index] = cnv$aggregate.length[loc.index] + hit$overlap[i]
    }
    
    # aggregate regions that are not covered by any segments, and divided
    # by gene length to get average copy number per gene
    cnv$aggregate.cn = cnv$aggregate.cn + 2 * (gene.length - cnv$aggregate.length)
    output[,patient.index] = cnv$aggregate.cn/gene.length
  }  
  
  # output
  out.file = paste0("./cnv/", disease, ".CN.txt")
  write.table(output, out.file, sep="\t", quote=F, col.names=T, row.names=F)
}

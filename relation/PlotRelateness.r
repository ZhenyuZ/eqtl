# PlotRelateness.r takes disease IBD measurement, filter by preset
# cutoff value (PI_HAT > 0.1875), and plot their relationship

# load library
library(qgraph)

# set PI_HAT cutoff value = 0.1875 (between 2 and 3 degree)
cutoff = 0.1875

# read disease types
diseases = unlist(read.table("../../diseases.txt", h=F, stringsAsFactors=F))

# read plink gnome IBD output and generate relationship plots by disease



for(disease in diseases) { 
  # read pair-wise plink genome IBD output
  ibd.file = paste0(disease, ".IBD.genome")
  ibd = data.frame(read.table(ibd.file, h=T, colClasses=c(rep("character", 6), rep("numeric",8))))
  
  # filter for PI_HAT > preset cutoff value
  w = which(ibd$PI_HAT > cutoff)
  if(length(w) == 0) {
    break
  }
  ibd = ibd[w, ]
  
  # define edges
  edges <- data.frame(
    from = ibd$IID1,
    to = ibd$IID2,
    thickness = ibd$PI_HAT
  )
  
  # plot relationship
  out.file = paste0(disease, ".relationship.png")
  png(out.file, width=960, height=960)
  qgraph(edges,esize=5,gray=TRUE)
  dev.off()
}



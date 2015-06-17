# Cytoband2Chunk.r take cytoband.txt and provide chunk ranges 
options(scipen=999)

# read cytoband info
cyto.file = "cytoBand.txt"
cyto = data.frame(read.delim(cyto.file, h=F, stringsAsFactors=F))
colnames(cyto) = c("chr", "start", "end", "band", "property")


output = matrix(NA, 0, 4)
colnames(output) = c("name", "chr", "start", "end")

# loop to fill output
for (i in 1:22) {
  # extract chr specific data
  chr = paste0("chr", i)
  c = cyto[cyto$chr == chr, ]
  
  # calculate ends of short and long arm
  pband = which(substr(c$band, 1, 1) == "p")
  qband = which(substr(c$band, 1, 1) == "q")
  minp = min(c$start[pband])
  maxp = max(c$end[pband])
  minq = min(c$start[qband])
  maxq = max(c$end[qband])  
  
  # sanity check
  if(maxp != minq) {
    print("Warning! centromere location does not match.")
  }
  
  # Calculate chunks
  chunk = 5000000
  center = maxp
  len = maxq
  
  # init chunk finding
  end = chunk
  s = end
  
  # chunk p-arm
  while (end < center){
    end = end + chunk
    s = append(s, end)
  }
  s = s[-length(s)]
  if (center - s[length(s)] + chunk < 7000000) {
    s[length(s)] = center
  } else {
    s = append(s, center)
  }
  end = center	
  
  # chunk q-arm
  while (end < len){
    end = end + chunk
    s = append(s, end)
  }
  s = s[-length(s)]
  if (len - s[length(s)]  + chunk < 7000000) {
    s[length(s)] = len
  } else {
    s = append(s, len)
  }
  
  # combine name chr start and end info in output
  names = paste0("chr", i, ".", (1:length(s)))
  output.block = cbind(names, i, c(1, (s + 1)[1:(length(s)-1)]), s)
  output = rbind(output, output.block)
}  
        
write.table(output, "chunk.range", col.names=F, row.names=F, sep="\t", quote=F)
  
  
  
  

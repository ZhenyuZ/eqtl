# Run PEER analsys with disease type and number of factor as input
# R CMD BATCH --no-save --no-restore '--args LUAD 60' GetPeerResiduals.r LUAD.60.out & 

args = commandArgs(T)
disease = args[1]
nPeer = as.numeric(args[2])

library(peer)
options(stringsAsFactors=F)

# my.write.table use write.table with both colnames and rownames
# with additional parameters sep="\t" and quote=F
# In addition, it add the field name of rownames at the begining 
# of the file, with default "id", so that all rows (including first)
# will have the same number of fields
my.write.table <- function(data, file, add="id"){
  # create and open the file connection
  f.con <- file(file, open = 'wt')
  # close on exit 
  on.exit(close(f.con))
  # write header
  header = paste(c(add, colnames(data)), collapse="\t")
  writeLines(header, con=f.con)
  # write data
  write.table(data, f.con, col.names=F, row.names=T, sep="\t", quote=F)
}


# read expr
expr.file = "expr.residuals.rda"
load(expr.file)

# read clin
clin.file = "clin.txt"
clin = read.table(clin.file, h=T, row.names=1, check.names=F)

# PEER analsys  
t = t(expr)
covs = t(clin)

model = PEER()
PEER_setPhenoMean(model, as.matrix(t))
PEER_setCovariates(model, as.matrix(covs))
PEER_setNk(model, nPeer)
PEER_getNk(model)
PEER_update(model)

# save PEER residuals
resPeer = t(PEER_getResiduals(model))
rownames(resPeer) = rownames(expr)
colnames(resPeer) = substr(colnames(expr), 1, 12)

outfile = paste0("GE-", nPeer, ".txt")

my.write.table(resPeer, outfile)







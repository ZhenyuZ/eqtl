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

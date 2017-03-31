

############  filter_out_strand_break()  ############################


filter_out_strand_break <- function(mat){
  strand_break_out <- as.character(sapply(as.character(colnames(mat)), function(x) return (paste0(strsplit(x, "_")[[1]][-3], collapse="_"))))
  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat[l,], strand_break_out, sum))
  }
  rownames(mat_filtered) <- rownames(mat)

  return(mat_filtered)
}

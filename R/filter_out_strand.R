

##############   filter_out_strand()  #############################

filter_out_strand <- function(mat){
  strand_out <- as.character(sapply(as.character(colnames(mat)), function(x) return (paste0(strsplit(x, "_")[[1]][-2], collapse="_"))))
  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat[l,], strand_out, sum))
  }
  rownames(mat_filtered) <- rownames(mat)

  return(mat_filtered)
}

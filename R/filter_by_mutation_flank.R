


############   filter_by_mutation_flank()  ###############################


filter_by_mutation_flank <-  function(mat){
  mutation_flank <- as.character(sapply(as.character(colnames(mat)), function(x) return (strsplit(x,"_")[[1]][1])))

  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat[l,], mutation_flank, sum))
  }
  rownames(mat_filtered) <- rownames(mat)

  return(mat_filtered)
}

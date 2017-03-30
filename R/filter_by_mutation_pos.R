


############   filter_by_mutation_pos()  ###############################


filter_by_mutation_pos <-  function(mat, max_pos = 10, flanking_bases=1){
  input_pos = 1:max_pos
  mutations <- as.character(sapply(as.character(colnames(mat)), function(x) return (paste0(strsplit(x, "")[[1]][(flanking_bases + 1): (flanking_bases + 4)], collapse=""))))
  pos <- as.numeric(sapply(as.character(colnames(mat)), function(x) return (strsplit(x,"_")[[1]][4])))
  indices <- which(!is.na(match(pos, input_pos)))

  mutations_pos <- paste0(mutations, "_", pos)
  mutations_pos_reduced <- mutations_pos[indices]
  mat_reduced <- mat[,indices]

  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat_reduced[l,], mutations_pos_reduced,  sum))
  }
  rownames(mat_filtered) <- rownames(mat)

  return(mat_filtered)
}

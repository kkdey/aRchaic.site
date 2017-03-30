
##############   filter_by_mutation_pos_flank()  #############################


filter_by_mutation_flank_pos <-  function(mat, max_pos = 20, flanking_bases=1){

  input_pos <- 1:max_pos
  mutation_flank <- as.character(sapply(as.character(colnames(mat)), function(x) return (strsplit(x,"_")[[1]][1])))
  pos <- as.numeric(sapply(as.character(colnames(mat)), function(x) return (tail(strsplit(x,"_")[[1]],1))))
  pos_indices <- which(!is.na(match(pos, input_pos)))
  mutation_flank_pos <- paste0(mutation_flank, "_", pos)
  mutation_flank_pos_reduced <- mutation_flank_pos[pos_indices]
  mat_reduced <- mat[,pos_indices]

  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat_reduced[l,], mutation_flank_pos_reduced,  sum))
  }
  rownames(mat_filtered) <- rownames(mat)

  return(mat_filtered)

}



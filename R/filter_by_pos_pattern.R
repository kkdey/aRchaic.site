

##############   filter_by_pos_pattern()  #############################

filter_by_pos_pattern <-  function(mat, max_pos = 20, pattern = "C->T", flanking_bases=1){

  input_pos <- 1:max_pos
  mutations <- as.character(sapply(as.character(colnames(mat)), function(x) return (paste0(strsplit(x, "")[[1]][(flanking_bases + 1): (flanking_bases + 4)], collapse=""))))
  CtoTindices <- which(mutations == pattern)

  pos <- as.numeric(sapply(as.character(colnames(mat)), function(x) return (tail(strsplit(x,"_")[[1]],1))))
  pos_indices <- which(!is.na(match(pos, input_pos)))

  matched_indices <- intersect(CtoTindices, pos_indices)

  mutations_pos <- paste0(mutations, "_", pos)
  mutations_pos_reduced <- mutations_pos[matched_indices]
  mat_reduced <- mat[,matched_indices]

  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat_reduced[l,], mutations_pos_reduced, sum))
  }
  rownames(mat_filtered) <- rownames(mat)
  return(mat_filtered)
}

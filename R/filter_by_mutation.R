


##############   filter_by_mutation()  ##############################


filter_by_mutation <-  function(mat, flanking_bases=1){
  mutations <- as.character(sapply(as.character(colnames(mat)), function(x) return (paste0(strsplit(x, "")[[1]][(flanking_bases + 1): (flanking_bases + 4)], collapse=""))))

  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat[l,], mutations, sum))
  }
  rownames(mat_filtered) <- rownames(mat)

  return(mat_filtered)
}

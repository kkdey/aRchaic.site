library(dplyr)
library(Logolas)
library(plyr)
library(grid)
library(gridBase)
library(limma)


aRchaic_pca_beta =  function(mat,
                             labs = NULL,
                             type = c("mutation", "mutation-flank", "mutation-flank-pos",
                                      "specific-mutation-pos", "wo-strand", "wo-strand-break"),
                             pattern = NULL,
                             flanking_bases = 1,
                             max_pos = 20,
                             normalize=TRUE,
                             cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                 "hotpink","burlywood","darkkhaki","yellow","darkgray","deepskyblue",
                                 "brown4","darkorchid","magenta", "azure1","azure4"),
                             pcs_to_plot = c("PC1", "PC2", "PC3"),
                             lay = NULL,
                             plot_width = 10,
                             plot_height = 7,
                             output_dir = NULL,
                             output_name = "pca"){

  if(type == "mutation"){
    mat_reduced <- filter_by_mutation(mat)
  }else if (type == "mutation-flank"){
    mat_reduced <- filter_by_mutation_flank(mat)
  }else if (type == "mutation-flank-pos"){
    mat_reduced <- filter_by_mutation_flank_pos(mat, max_pos, flanking_bases)
  } else if (type == "specific-mutation-pos"){
    if(is.null(pattern)){
      stop("for this option, user needs to provide a mutation type: C->T, C->A C->G, T->A, T->C or T->G")
    }
    mat_reduced <- filter_by_pos_pattern(mat, max_pos, pattern = pattern, flanking_bases)
  }else if (type == "wo-strand"){
    mat_reduced <- filter_out_strand(mat)
  }else if (type == "wo-strand-break"){
    mat_reduced <- filter_out_strand_break(mat)
  }else {
    stop("the type of filtering does not match with the possible options")
  }

  rownames(mat_reduced) <- rownames(mat)
  if(normalize){
    voom_mat_reduced <- t(limma::voom(t(mat_reduced))$E);
    pr <- prcomp(voom_mat_reduced)
  }else{
    pr <- prcomp(mat_reduced)
  }

 if(is.null(labs)){
   labs <- rep(1, dim(voom_mat_reduced)[1])
 }

  pc_data_frame <- cbind.data.frame(pr$x, labs)

  graphList <- vector(mode="list")

  indices <- array(0, length(pcs_to_plot))

  pc_data_frame <- cbind.data.frame(pr$x, labs)

  for(num in 1:length(pcs_to_plot)){
    indices[num] <- match(paste0(pcs_to_plot[num]), colnames(pc_data_frame))
  }

  pc_data_frame_filtered <- pc_data_frame[,indices]

  total_comb <- length(pcs_to_plot)*(length(pcs_to_plot) - 1)/2

  if(is.null(lay)){
    lay <- matrix(0, length(pcs_to_plot) - 1, length(pcs_to_plot) - 1)
    l <- 1
    for(m in 1:nrow(lay)){
      for(n in 1:ncol(lay)){
        if(m > n){
          lay[m,n] <- NA
        }else{
          lay[m,n] = l
          l <- l+1
        }
      }
    }
  }

  lnum <- 1
  for(m in 1:(length(pcs_to_plot)-1)){
    for(n in (m+1):(length(pcs_to_plot))){
      a <- ggplot2::ggplot(pc_data_frame_filtered, ggplot2::aes_string(x = paste0("PC", indices[m]), y = paste0("PC", indices[n]))) +
        geom_point(aes(colour = factor(labs))) + ggplot2::scale_colour_manual(values = cols, guide = ggplot2::guide_legend(title = "Groups")) +
        ggplot2::xlab(paste0(pcs_to_plot[m])) + ggplot2::ylab(paste0(pcs_to_plot[n]))
      graphList[[lnum]] <- a
      lnum = lnum +1
    }
  }

  plot.new()
  grid.newpage()
  a <- do.call("grid.arrange",
               args = list(grobs=graphList,
                           layout_matrix = lay))

  ggplot2::ggsave(paste0(output_dir, output_name, ".png"), a, width = plot_width,
                  height = plot_height)


}

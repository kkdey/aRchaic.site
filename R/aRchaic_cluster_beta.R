library(dplyr)
library(Logolas)
library(plyr)
library(grid)
library(gridBase)
library(limma)

aRchaic_cluster_beta = function(mat,
                                K,
                                type = c("mutation", "mutation-flank", "mutation-pos", "mutation-flank-pos",
                                         "specific-mutation-pos", "wo-strand", "wo-strand-break"),
                                tol=0.01,
                                pattern = NULL,
                                max_pos = 20,
                                labs = NULL,
                                levels = NULL,
                                flanking_bases = 1,
                                gom_method = "independent",
                                topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                               "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                                structure.control = list(),
                                graph.control = list(),
                                topics.control = list(),
                                output_dir = NULL,
                                structure_width = 5,
                                structure_height = 8,
                                inflation = rep(2,1,2),
                                output_width = 1200,
                                output_height = 700){

if(type == "mutation"){
  aRchaic_cluster_beta_mutation(mat=mat, K=K, tol=tol,
                                labs = labs, levels = levels, flanking_bases = flanking_bases,
                                gom_method = gom_method, topic_cols = topic_cols,
                                structure.control = structure.control,
                                logo.control = graph.control, topics.control = topics.control,
                                output_dir = output_dir, structure_width = structure_width,
                                structure_height = structure_height,
                                output_width = output_width, output_height = output_height)
}
else if(type == "mutation-flank"){
  aRchaic_cluster_beta_mutation_flank(mat=mat, K=K, tol=tol,
                                labs = labs, levels = levels, flanking_bases = flanking_bases,
                                gom_method = gom_method, topic_cols = topic_cols,
                                structure.control = structure.control,
                                logo.control = graph.control, topics.control = topics.control,
                                output_dir = output_dir, structure_width = structure_width,
                                structure_height = structure_height, inflation = inflation,
                                output_width = output_width, output_height = output_height)
}
else if(type == "mutation-pos"){
  aRchaic_cluster_beta_mutation_pos(mat=mat, K=K, tol=tol, max_pos = max_pos,
          labs = labs, levels = levels, flanking_bases = flanking_bases,
          gom_method = gom_method, topic_cols = topic_cols, structure.control = structure.control,
          logo.control = graph.control, topics.control = topics.control,
          output_dir = output_dir, structure_width = structure_width,
          structure_height = structure_height,
           output_width = output_width, output_height = output_height)
}
else if(type == "mutation-flank-pos"){
  aRchaic_cluster_beta_mutation_flank_pos(mat=mat, K=K, tol=tol, max_pos = max_pos,
                                    labs = labs, levels = levels, flanking_bases = flanking_bases,
                                    gom_method = gom_method, topic_cols = topic_cols, structure.control = structure.control,
                                    logo.control = graph.control, topics.control = topics.control,
                                    output_dir = output_dir, structure_width = structure_width,
                                    structure_height = structure_height, inflation = inflation,
                                    output_width = output_width, output_height = output_height)
}
else if(type == "specific-mutation-pos"){
  if(is.null(pattern)){
    stop("for this option, user needs to provide a mutation type: C->T, C->A C->G, T->A, T->C or T->G")
  }
  aRchaic_cluster_beta_pos(mat=mat, pattern = pattern, K=K, tol=tol, max_pos = max_pos,
                           labs = labs, levels = levels, flanking_bases = flanking_bases,
                           gom_method = gom_method, topic_cols = topic_cols,
                           structure.control = structure.control,
                           graph.control = graph.control,
                           topics.control = topics.control,
                           output_dir = output_dir, structure_width = structure_width,
                           structure_height = structure_height,
                           output_width = output_width, output_height = output_height)
}
else if(type == "wo-strand"){
  aRchaic_cluster_beta_wo_strand(mat=mat, K=K, tol=tol, max_pos = max_pos,
                                labs = labs, levels = levels, flanking_bases = flanking_bases,
                                gom_method = gom_method, topic_cols = topic_cols, structure.control = structure.control,
                                logo.control = graph.control, topics.control = topics.control,
                                output_dir = output_dir, structure_width = structure_width,
                                tructure_height = structure_height, inflation = inflation,
                                output_width = output_width, output_height = output_height)
}
else if(type == "wo-strand-break"){
    aRchaic_cluster_beta_wo_strand_break(mat=mat, K=K, tol=tol, max_pos = max_pos,
                                        labs = labs, levels = levels, flanking_bases = flanking_bases,
                                        gom_method = gom_method, topic_cols = topic_cols, structure.control = structure.control,
                                        logo.control = graph.control, topics.control = topics.control,
                                        output_dir = output_dir, structure_width = structure_width,
                                        tructure_height = structure_height, inflation = inflation,
                                        output_width = output_width, output_height = output_height)
}else{
  stop("the type of filtering does not match with the possible options: see documentation")
  }
}

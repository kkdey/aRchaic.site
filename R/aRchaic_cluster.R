


######################################  aRchaic cluster   ###########################################################################3

library(dplyr)
library(Logolas)
library(plyr)
library(grid)
library(gridBase)



aRchaic_cluster = function(folders,
                           K,
                           tol=100,
                           labs = NULL,
                           levels = NULL,
                           run_from = c("start", "gom", "plot"),
                           run_index = NULL,
                           breaks = c(-1, seq(1,20,1)),
                           flanking_bases = 1,
                           gom_method = "independent",
                           topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                      "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                      "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                           structure.control = list(),
                           logo.control = list(),
                           topics.control = list(),
                           output_dir = NULL,
                           save_plot = TRUE,
                           structure_width = 5,
                           structure_height = 8){

  structure.control.default <- list(yaxis_label = "aRchaic pops",
                                 order_sample = FALSE,
                                 figure_title = paste0("  StructurePlot: K=", K,""),
                                 axis_tick = list(axis_ticks_length = .1,
                                                  axis_ticks_lwd_y = .1,
                                                  axis_ticks_lwd_x = .1,
                                                  axis_label_size = 7,
                                                  axis_label_face = "bold"),
                                 legend_title_size = 8,
                                 legend_key_size = 0.4,
                                 legend_text_size = 5)

  logo.control.default <- list(sig_names = NULL, ic.scale=TRUE,
                               max_pos = 20, flanking_bases=1,
                               yscale_change = TRUE, xaxis=TRUE,
                               yaxis=TRUE, xlab = " ", xaxis_fontsize=20,
                               xlab_fontsize=10, title_aligner = 11,
                               y_fontsize=20, title_fontsize = 35,
                               mut_width=2, start=0.0001,
                               renyi_alpha = 5, inflation_factor = c(3,1,3),
                               pop_names = paste0("Cluster : ", 1:K),
                               logoport_x = 0.25, logoport_y= 0.50, logoport_width= 0.25, logoport_height= 0.50,
                               lineport_x = 0.9, lineport_y=0.40, lineport_width=0.32, lineport_height=0.28,
                               breaklogoport_x = 0.94, breaklogoport_y = 0.40, breaklogoport_width=0.30, breaklogoport_height=0.45,
                               barport_x = 0.60, barport_y=0.60, barport_width=0.25, barport_height=0.25,
                               output_width = 1200, output_height = 700)

  topics.control.default <- list(bf = FALSE, kill = 2, ord = TRUE, verb = 1, admix = TRUE,
                        nbundles = 1, use_squarem = FALSE, init.adapt = FALSE, type = "full",
                        light = 1, method_admix = 1, sample_init = TRUE, tmax = 10000)

  structure.control <- modifyList(structure.control, structure.control.default)
  logo.control <- modifyList(logo.control, logo.control.default)
  topics.control <- modifyList(topics.control, topics.control.default)



  if(is.null(labs)){
    labs <- c()
    for(i in 1:length(folders)){
      temp <- setdiff(list.files(folders[i], pattern = ".csv"), list.files(folders[i], pattern = ".csv#"))
      labs <- c(labs, rep(tail(strsplit(folders[i], "/")[[1]],1), length(temp)))
    }
  }

  if(is.null(levels)){
    levels <- unique(labs)
  }

  ##########################  Check if the folder names actually exist #######################################

  message("Checking if the folders exist")

  for(i in 1:length(folders)){
    if(!file.exists(folders[i]))
      stop("A folder in the folder list does not exist:  aborting")
  }

  datalist <- vector("list", length(folders))

  #########################  If the user wants to run from scratch  ##########################################

  if(run_from == "start"){
    if(is.null(run_index)){
      run_index <- 1:length(folders)
    }
    if(sum(run_index - 1:length(folders))^2 == 0){
      for(i in 1:length(folders)){
        file.remove(paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda"))
      }
    }else{
      folders1 <- folders[run_index]
      for(i in 1:length(folders1)){
        file.remove(paste0(folders1[i], tail(strsplit(folders1[i], "/")[[1]],1), ".rda"))
      }
    }
  }

  #########################  Run aggregation functions on MutationFeatureFormat #############################################


  for(i in 1:length(folders)){
    if(!file.exists(paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda"))){
      message ("Processing the MutationFeatureFormat files in the directory")
      out <- aggregate_signature_counts(dir = paste0(folders[i]),
                                        pattern = NULL,
                                        breaks = breaks,
                                        flanking_bases = flanking_bases)
      clubbed_data <- club_signature_counts(out, flanking_bases = 1)
      save(clubbed_data, file = paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda"))
      datalist[[i]] <- clubbed_data
    }else{
      datalist[[i]] <- get(load(paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda")))
    }
  }

  ######################  Pooling data from different folders if present  ############################################

  if(run_from == "start"){
    message("Pooling the data from multiple sources")
  }

  sig_names <- colnames(datalist[[1]])
  row_names_pool <- rownames(datalist[[1]])
  if(length(datalist) >= 2){
    for(num in 2:length(datalist)){
      sig_names <- union(sig_names, colnames(datalist[[num]]))
      row_names_pool <- c(row_names_pool, rownames(datalist[[num]]))
    }
  }

  pooled_data <- matrix(0, length(row_names_pool), length(sig_names))
  rownames(pooled_data) <- row_names_pool
  colnames(pooled_data) <- sig_names

  for(num in 1:length(datalist)){
    pooled_data[match(rownames(datalist[[num]]), rownames(pooled_data)), match(colnames(datalist[[num]]), sig_names)] <- datalist[[num]]
  }

  zero_sum_rows <- which(rowSums(pooled_data) == 0)
  if(length(zero_sum_rows) > 0){
    pooled_data <- pooled_data[-zero_sum_rows, ]
    labs <- labs[-zero_sum_rows]
  }


  ########################  Grade of Membership Model  ######################################################

  if(is.null(output_dir)){
    output_dir <- paste0(getwd(), "/")
  }

  if(run_from == "gom" | run_from == "start"){
    if(gom_method == "full"){
        message("Fitting the Grade of Membership Model - full version - due to Matt Taddy")
        suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = pooled_data, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
        save(topic_clus, file = paste0(output_dir, "model.rda"))
    }else if(gom_method == "independent"){
      message("Fitting the Grade of Membership Model - full version - due to Y. Shiraichi and M. Stephens")

      signature_set <- colnames(pooled_data)
      sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:8])))
      new_sig_split <- matrix(0, dim(sig_split)[1], 3);
      new_sig_split[,1] <- sig_split[,1]
      new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,2:5], collapse="")))
      new_sig_split[,3] <- sig_split[,6]

      levels(new_sig_split[,1]) <- c("0", "1", "2", "3", "4")

      pos <- t(sapply(1:length(signature_set), function(x)
      {
        y = strsplit(signature_set[x], "")[[1]]
        return(paste(y[12:length(y)], collapse=""))
      }))



      mat <- matrix(0, dim(new_sig_split)[1], dim(new_sig_split)[2])
      for(k in 1:dim(new_sig_split)[2]){
        temp <- as.factor(new_sig_split[,k])
        mat[,k] <- as.numeric(as.matrix(plyr::mapvalues(temp, from = levels(temp), to = 0:(length(levels(temp))-1))))
      }

      pos <- as.numeric(pos)
      pos <- pos - min(pos)
      pos <- factor(pos, levels = 0:21)

      signatures <- mat;
      signature_pos <- cbind.data.frame(signatures, pos)

      suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = pooled_data, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
      save(topic_clus, file = paste0(output_dir, "model.rda"))
    }
  }else if(run_from == "plot"){
    if(gom_method == "full"){
      if(!file.exists(paste0(output_dir, "model.rda"))){
        message("Fitting the Grade of Membership Model - full version - due to Matt Taddy")
        suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = pooled_data, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
        save(topic_clus, file = paste0(output_dir, "model.rda"))
      }else{
        topic_clus <- get(load(paste0(output_dir, "model.rda")))
      }
    }else if(gom_method == "independent"){
      message("Fitting the Grade of Membership Model - independent version - due to Y. Shiraichi and M. Stephens")
      signature_set <- colnames(pooled_data)
      sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:8])))
      new_sig_split <- matrix(0, dim(sig_split)[1], 3);
      new_sig_split[,1] <- sig_split[,1]
      new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,2:5], collapse="")))
      new_sig_split[,3] <- sig_split[,6]

      levels(new_sig_split[,1]) <- c("0", "1", "2", "3", "4")

      pos <- t(sapply(1:length(signature_set), function(x)
      {
        y = strsplit(signature_set[x], "")[[1]]
        return(paste(y[12:length(y)], collapse=""))
      }))



      mat <- matrix(0, dim(new_sig_split)[1], dim(new_sig_split)[2])
      for(k in 1:dim(new_sig_split)[2]){
        temp <- as.factor(new_sig_split[,k])
        mat[,k] <- as.numeric(as.matrix(plyr::mapvalues(temp, from = levels(temp), to = 0:(length(levels(temp))-1))))
      }

      pos <- as.numeric(pos)
      pos <- pos - min(pos)
      pos <- factor(pos, levels = 0:21)

      signatures <- mat;
      signature_pos <- cbind.data.frame(signatures, pos)

      if(!file.exists(paste0(output_dir, "model.rda"))){
        suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = pooled_data, K=K, tol=tol, model = "independent", signatures = signature_pos), topics.control)))
        save(topic_clus, file = paste0(output_dir, "model.rda"))
      }else{
        topic_clus <- get(load(paste0(output_dir, "model.rda")))
      }
    }

  }else {
    stop("run from must be either from start (which clears everything out and restarts) or from gom (which does clustering and follow up)
         or from plot (which just does the plots)")
  }


  message ("Structure plot and Logo plot representations of clusters")

  omega <- topic_clus$omega
  annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(labs, levels = levels)
  )

  if(save_plot){
    if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}
    plot.new()
    grid.newpage()
    do.call(StructureGGplot, append(list(omega= omega,
                                         annotation = annotation,
                                         palette = topic_cols),
                                                structure.control))
    ggplot2::ggsave(paste0(output_dir, "structure.png"), width = structure_width,
                    height = structure_height)
  }

  if(save_plot){
    if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}
    plot.new()
    do.call(damageLogo_five, append(list(theta_pool = topic_clus$theta,
                                output_dir = output_dir, save_plot = save_plot),
            logo.control))
    graphics.off()
  }else if(!save_plot){
    plot.new()
    if(is.null(output_dir)){ output_dir <- paste0(getwd(), "/")}
    do.call(StructureGGplot, append(list(omega= omega, annotation = annotation, palette = topic_cols), structure.control))
    do.call(damageLogo_five, append(list(theta_pool = topic_clus$theta, output_dir = output_dir, save_plot = save_plot),
            logo.control))
  }

  message ("Finished")
}

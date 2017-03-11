

###################   aRchaic view operator   ################################################33


library(dplyr)
library(Logolas)
library(plyr)
library(grid)
library(gridBase)
library(ggplot2)


aRchaic_view = function(dir,
                        file,
                        breaks = c(-1, seq(1,20,1)),
                        flanking_bases =1,
                        logo.control = list(),
                        title = NULL,
                        output_dir = NULL,
                        save_plot = TRUE){

  if(is.null(title)){
    title <- strsplit(file, ".csv")[[1]][1]
  }
  logo.control.default <- list(sig_names = NULL, ic.scale=TRUE,
                               max_pos = 20, flanking_bases=1,
                               yscale_change = TRUE, xaxis=TRUE,
                               yaxis=TRUE, xlab = " ", xaxis_fontsize=5,
                               xlab_fontsize=10, title_aligner = 18,
                               y_fontsize=10, title_fontsize = 20,
                               mut_width=2, start=0.0001,
                               renyi_alpha = 1, pop_names=title,
                               logoport_x = 0.25, logoport_y= 0.50, logoport_width= 0.28, logoport_height= 0.40,
                               lineport_x = 0.9, lineport_y=0.40, lineport_width=0.32, lineport_height=0.28,
                               breaklogoport_x = 1.00, breaklogoport_y = 0.40, breaklogoport_width=0.40, breaklogoport_height=0.50,
                               barport_x = 0.58, barport_y=0.60, barport_width=0.25, barport_height=0.25,
                               output_width = 1200, output_height = 700)

  logo.control <- modifyList(logo.control, logo.control.default)


  if(file.exists(paste0(dir, tail(strsplit(dir, "/")[[1]],1), ".rda"))){
    message("MutationFeatureFormat file present: skipping the signature aggregation step")
    mff_dat <- get(load(paste0(dir, tail(strsplit(dir, "/")[[1]],1), ".rda")))
    index <- grep(paste0(file), rownames(mff_dat))
    clubbed_counts <- mff_dat[index, ]
    clubbed_counts_norm <- clubbed_counts/ sum(clubbed_counts)

    temp <- as.matrix(clubbed_counts_norm)
    rownames(temp) <- names(clubbed_counts_norm)

  }else{
    message("MutationFeatureFormat file not present: performing the signature aggregation step ")
    pattern =file
    out <- aggregate_signature_counts(dir = dir,
                                      pattern = pattern,
                                      breaks = breaks,
                                      flanking_bases = flanking_bases)
    clubbed_counts <- club_signature_counts(out, flanking_bases = 1)
    clubbed_counts_norm <- clubbed_counts/ sum(clubbed_counts)

    temp <- t(clubbed_counts_norm)
  }

  if(is.null(output_dir)){
    output_dir <- paste0(getwd(), "/")
  }


  do.call(damageLogo_five, append(list(theta_pool = temp,
                                       output_dir = output_dir,
                                       save_plot = save_plot),
                                  logo.control))
  }

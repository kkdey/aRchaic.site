
StructureGGplot <- function(omega, annotation = NULL,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            figure_title = "",
                            yaxis_label = "Tissue type",
                            order_sample = TRUE,
                            sample_order_decreasing = TRUE,
                            split_line = list(split_lwd = 1,
                                              split_col = "white"),
                            plot_labels = TRUE,
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 3,
                                             axis_label_face = "bold"),
                            legend_title_size = 8,
                            legend_key_size = 0.4,
                            legend_text_size = 5,
                            output_dir = NULL,
                            output_width = 400,
                            output_height = 700) {


  # check if the number of colors is same as or more than the number of clusters
  if (dim(omega)[2] > length(palette)) {
    stop("Color choices is smaller than the number of clusters!")
  }

  # check if rownames of omega are unique
  if(length(unique(rownames(omega))) != NROW(omega)) {
    stop("omega rownames are not unique!")
  }

  # check the annotation data.frame
  if (is.null(annotation)) null_annotation <- TRUE
  if (!is.null(annotation)) null_annotation <- FALSE

  if (null_annotation) {
    annotation <- data.frame(
      sample_id = paste("X", c(1:NROW(omega))),
      tissue_label = rep("NA", NROW(omega)) )
  } else if (!null_annotation) {
    if (!is.data.frame(annotation))
      stop("annotation must be a data.frame")
    if (!all.equal(colnames(annotation), c("sample_id", "tissue_label")) ) {
      stop("annotation data.frame column names must be sample_id and tissue_label")
    }
    if ( length(unique(annotation$sample_id)) != NROW(omega)) {
      stop("sample_id is not unique")
    }
  }

  df_ord <- do.call(rbind,
                    lapply(1:nlevels(annotation$tissue_label), function(ii) {
                      temp_label <- levels(annotation$tissue_label)[ii]
                      temp_df <- omega[which(annotation$tissue_label == temp_label), , drop=FALSE]

                      is_single_sample <- (nrow(temp_df) == 1)
                      # find the dominant cluster in each sample
                      if ( is_single_sample ) {
                        each_sample_order <- which.max(temp_df)
                      } else {
                        each_sample_order <- apply(temp_df, 1, which.max)
                      }

                      # find the dominant cluster across samples
                      sample_order <- as.numeric(attr(table(each_sample_order), "name")[1])

                      if (order_sample == TRUE & !is_single_sample) {
                        # reorder the matrix
                        temp_df_ord <- temp_df[order(temp_df[ , sample_order],
                                                     decreasing = sample_order_decreasing), ]
                      } else {
                        temp_df_ord <- temp_df
                      }
                      temp_df_ord
                    }) )

  df_mlt <- reshape2::melt(t(df_ord))
  df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
                                             "Var2" = "document"))
  df_mlt$document <- factor(df_mlt$document)
  df_mlt$topic <- factor(df_mlt$topic)

  # set blank background
  ggplot2::theme_set(ggplot2::theme_bw(base_size = 12)) +
    ggplot2::theme_update( panel.grid.minor.x = ggplot2::element_blank(),
                           panel.grid.minor.y = ggplot2::element_blank(),
                           panel.grid.major.x = ggplot2::element_blank(),
                           panel.grid.major.y = ggplot2::element_blank() )

  # inflat nubmers to avoid rounding errors
  value_ifl <- 10000

  # number of ticks for the weight axis, including 0 and 1
  ticks_number <- 6

  # set axis tick positions
  tissue_count <- table(droplevels(annotation$tissue_label))
  tissue_count_cumsum <- cumsum(table(droplevels(annotation$tissue_label)))
  tissue_names <- levels(droplevels(annotation$tissue_label))

  # if more than 2 levels in the phenotype of interest
  if (length(tissue_names) > 1) {

    tissue_breaks <- sapply(1:length(tissue_count), function(i) {
      if (i == 1) {
        if (tissue_count[i] == 1) bk <- 1
        if (tissue_count[i] > 1)  bk <- (tissue_count_cumsum[i] - 0)/2
        return(bk)
      }
      if (i > 1) {
        if (tissue_count[i] == 1) bk_interval <- 1
        if (tissue_count[i] > 1 ) {
          bk_interval <- (tissue_count_cumsum[i] - tissue_count_cumsum[i-1])/2 }
        bk <- tissue_count_cumsum[i-1] + bk_interval
        return(bk)
      }
    })
    names(tissue_breaks) <- tissue_names

    # make ggplot
    a <- ggplot2::ggplot(df_mlt,
                         ggplot2::aes(x = df_mlt$document,
                                      y = df_mlt$value*10000,
                                      fill = factor(df_mlt$topic)) ) +
      ggplot2::xlab(yaxis_label) + ggplot2::ylab("") +
      ggplot2::scale_fill_manual(values = palette) +
      ggplot2::theme(legend.position = "right",
                     legend.key.size = ggplot2::unit(legend_key_size, "cm"),
                     legend.text = ggplot2::element_text(size = legend_text_size),
                     ##<-- TBD: center legend title
                     #              legend.title = element_text(hjust = 1),
                     axis.text = ggplot2::element_text(size = axis_tick$axis_label_size,
                                                       face = axis_tick$axis_label_face),
                     axis.ticks.y = ggplot2::element_line(size = axis_tick$axis_ticks_lwd_y),
                     axis.ticks.x = ggplot2::element_line(size = axis_tick$axis_ticks_lwd_x),
                     axis.ticks.length = ggplot2::unit(axis_tick$axis_ticks_length, "cm"),
                     title = ggplot2::element_text(size = legend_title_size) ) +
      ggplot2::ggtitle(figure_title) +
      ggplot2::scale_y_continuous( breaks = seq(0, value_ifl, length.out = ticks_number),
                                   labels = seq(0, 1, 1/(ticks_number -1 ) ) ) +
      # Add tissue axis labels
      # ggplot2::scale_x_discrete(breaks = as.character(as.numeric(levels(df_mlt$document)[round(tissue_breaks)])),
      #                           labels = names(tissue_breaks)) +
      ggplot2::scale_x_discrete(breaks = as.character((levels(df_mlt$document)[round(tissue_breaks)])),
                                labels = names(tissue_breaks)) +
      # Add legend title
      ggplot2::labs(fill = "Clusters") +
      ggplot2::coord_flip()


    # width = 1: increase bar width and in turn remove space
    # between bars
    b <- a + ggplot2::geom_bar(stat = "identity",
                               position = "stack",
                               width = 1)
    # sample labels option
    if (plot_labels == TRUE) {
      b
    } else {
      b <- b + theme(axis.text.y = element_blank())
    }

    # remove plot border
    b <- b + cowplot::panel_border(remove = TRUE)

    # Add demarcation
    b <- b + ggplot2::geom_vline(
      xintercept = cumsum(table(droplevels(annotation$tissue_label)))[
        -length(table(droplevels(annotation$tissue_label)))] + .5,
      col = split_line$split_col,
      size = split_line$split_lwd)
    b
  } else if (null_annotation) {
    # make ggplot
    a <- ggplot2::ggplot(df_mlt,
                         ggplot2::aes(x = df_mlt$document,
                                      y = df_mlt$value*10000,
                                      fill = factor(df_mlt$topic)) ) +
      ggplot2::xlab(yaxis_label) + ggplot2::ylab("") +
      ggplot2::scale_fill_manual(values = palette) +
      ggplot2::theme(legend.position = "right",
                     legend.key.size = ggplot2::unit(legend_key_size, "cm"),
                     legend.text = ggplot2::element_text(size = legend_text_size),
                     ##<-- TBD: center legend title
                     #              legend.title = element_text(hjust = 1),
                     axis.text = ggplot2::element_text(size = axis_tick$axis_label_size,
                                                       face = axis_tick$axis_label_face),
                     axis.ticks.y = ggplot2::element_line(size = axis_tick$axis_ticks_lwd_y),
                     axis.ticks.length = ggplot2::unit(axis_tick$axis_ticks_length, "cm"),
                     title = ggplot2::element_text(size = legend_title_size) ) +
      ggplot2::ggtitle(figure_title) +
      ggplot2::scale_y_continuous( breaks = seq(0, value_ifl, length.out = ticks_number),
                                   labels = seq(0, 1, 1/(ticks_number -1 ) ) ) +
      ggplot2::scale_x_discrete(breaks = NULL) +
      # Add legend title
      ggplot2::labs(fill = "Clusters") +
      ggplot2::coord_flip()

    # width = 1: increase bar width and in turn remove space
    # between bars
    b <- a + ggplot2::geom_bar(stat = "identity",
                               position = "stack",
                               width = 1)
    # sample labels option
    if (plot_labels == TRUE) {
      b <- b + cowplot::panel_border(remove = TRUE)
     # png(paste0(output_dir, "structure.png"), width = output_width, height = output_height)
     #  plot.new()
      return(b)
     # graphics.off()
    } else {
      b <- b + theme(axis.text.y = element_blank())
      return(b)
     # png(paste0(output_dir, "structure.png"), width = output_width, height = output_height)
     # plot.new()
     # graphics.off()
    }

    # remove plot border
    #b <- b + cowplot::panel_border(remove = TRUE)
    #b

    #   if(!save_structure){
    #     print(b)
    #   }else{
    #     filename = paste0(output_dir, "structure.png")
    #     png(filename, width = output_width, height = output_height)
    #     print(b)
    #     dev.off()
    #   }
    #
  }

  # if (!plot_labels) {
  #     b
  # } else {
  #     b <- cowplot::ggdraw(cowplot::switch_axis_position((b), axis = "y"))
  #     b
  # }
}

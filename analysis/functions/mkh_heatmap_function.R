# mkh heatmap plotting function
# gene is the WT gene sequence
# order is the order of the variants on the y-axis
# chunk_size is the number of AA to plot per row
# scores_df is the dataframe containing the scores
# score_column is the column in scores_df containing the scores
# variant_names is an optional vector of names to use for the variants
# output_file is the file to save the plot to

# example variant order and label vectors below for missense only and full variant set with indels

order_missense <- c("A", "G", "M", "V", "L", "I", "T", "S", "C", "Q", "N", "Y", "W", "F", "E", "D", "H", "K", "R", "P")
order_full <- c("A", "G", "M", "V", "L", "I", "T", "S", "C", "Q", "N", "Y", "W", "F", "E", "D", "H", "K", "R", "P", "D_1", "D_2", "D_3", "I_1", "I_2", "I_3")

label_missense <- c("A", "G", "M", "V", "L", "I", "T", "S", "C", "Q", "N", "Y", "W", "F", "E", "D", "H", "K", "R", "P")
label_full <- c("A", "G", "M", "V", "L", "I", "T", "S", "C", "Q", "N", "Y", "W", "F", "E", "D", "H", "K", "R", "P", "Del x1", "Del x2", "Del x3", "Ins x1(G)", "Ins x2(GS)", "Ins x3(GSG)")

create_heatmap <- function(gene, order, chunk_size, scores_df, score_column, variant_names = NULL, output_file = NULL) {
  
  # Check if variant_names were provided; if not, use the original variant names
  if (is.null(variant_names)) {
    variant_names <- order
  }
  
  scores_df <- scores_df %>%
    filter(variants %in% order)
  
  gene_chunks <- str_split(gene, '')[[1]]
  num_chunks <- ceiling(length(gene_chunks) / chunk_size)
  
  min_score <- min(scores_df[, score_column], na.rm = TRUE)
  max_score <- max(scores_df[, score_column], na.rm = TRUE)
  
  heatmap_plots <- list()
  
  for (i in 1:num_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, length(gene_chunks))
    
    gene_chunk <- gene_chunks[start_idx:end_idx]
    scores_chunk <- scores_df[scores_df$pos %in% c(start_idx:(start_idx + length(gene_chunk) - 1)), ]
    
    # Replace the variant names with the desired names
    scores_chunk$variants <- factor(scores_chunk$variants, levels = order, labels = variant_names)
    
    
    heatmap_plot <- ggplot(data = scores_chunk, aes(x = pos, y = factor(variants, level = variant_names), fill = !!sym(score_column))) +
      geom_tile(aes(color = "grey"), linewidth = 0.2, position = "identity") +
      scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, l1 = 0.2, l3 = 0.2, p1 = 0.9, p3 = .4, p4 = 0.7, rev = FALSE, na.value = 'grey', limits = c(min_score, max_score)) +
      theme(
        panel.background = element_rect(fill = "grey", size = 0.1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
        panel.grid.minor = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(family = "sans", size = 9, angle = 0, hjust = 0.5, vjust = 1, margin = margin(t = 1)),
        axis.text = element_text(family = "sans", size = 9, color = "black"),
        axis.text.x.top = element_text(family = "sans",  size = 9, angle = 0, hjust = 0.5, margin = margin(r = 0)),
        axis.text.y = element_text(family = "sans",  size = 9,  margin = margin(r = 0)) ) +
      scale_x_continuous(breaks = seq(0, nchar(gene), by = 5), expand = c(0, 0),
                         sec.axis = sec_axis(trans = ~., name = "Sequence", breaks = seq(start_idx, end_idx), labels = gene_chunk, guide = derive())) +
      coord_fixed(ratio = 1) +
      scale_color_manual(values = c(NA, 'grey'), guide = FALSE) +
      labs(y = "Mutation", x = "Position") +
      geom_tile(data = subset(scores_chunk, !is.wt), aes(color = "grey"), size = 0.2, position = "identity", show.legend = FALSE) +
      geom_tile(data = subset(scores_chunk, is.wt), aes(color = "grey"), size = 0.2, position = "identity", show.legend = FALSE)
    
    heatmap_plots[[i]] <- heatmap_plot
  }
  heatmap_combined <- ggarrange(plotlist = heatmap_plots, nrow = num_chunks, ncol = 1)
  ggsave(output_file, heatmap_combined, height = 10, width = 20, dpi = 450)
}
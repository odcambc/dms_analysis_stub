library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(tidyverse)
library(colorspace)

# Function for plotting and saving heatmaps.
# label is the column for plotting
# sequence is a string of the WT sequence
# output_file is the location to save
# the p#_in parameters set the palette parameters:
# p1: Power parameter for chroma coordinates in first sequential palette
# p2: Power parameter for luminance coordinates in first sequential palette
# p3: Power parameter for chroma coordinates in second sequential palette
# p4: Power parameter for luminance coordinates in second sequential palette
# high: upper limit for plotting scores. Defaults to 99th percentile of scores.
# low: lower limit for plotting scores. Defaults to 1st percentile of scores
# print: if TRUE, it also prints a plot rather than just saving it. Note that
#        the dimensions might be off.

## TODO: more gracefully handle missing values. If they're NA, they're colored.
## If they are just not included, they're empty.

## TODO: the WT squares are outlined in green, but they get overdrawn.

default_order <- c("G", "A", "I", "L", "V",
                    "C", "M",
                    "S", "T", 
                    "K", "R", "H",
                    "D", "E", "N", "Q",
                    "F", "W", "Y",
                    "P",
                    "X",
                    "D_1", "D_2", "D_3",
                    "I_1", "I_2", "I_3")

default_variant_names <- c("G", "A", "I", "L", "V",
                    "C", "M",
                    "S", "T", 
                    "K", "R", "H",
                    "D", "E", "N", "Q",
                    "F", "W", "Y",
                    "P",
                    "*",
                    "Del x1", "Del x2", "Del x3",
                    "Ins x1 (G)", "Ins x2 (GS)", "Ins x3 (GSG)")

print_heatmap <- function(df, label, sequence, output_file = NULL, low = NA, 
                          high = NA, order_in = default_order, 
                          names = default_variant_names, p1_in = 0.9, 
                          p2_in = NA, p3_in = 0.4, p4_in = NA, print = FALSE,
                          invert_scale = FALSE, scale = FALSE, 
                          row_length = 100) {
  
  if (scale == TRUE) {
    scale_position = "right"
  } else {
    scale_position = "none"
  }
  
  legend_position = "none"
  
  label <- enquo(label)
  
  sequence_len = max(df$pos, na.rm = TRUE)
  
  strips = ceiling(sequence_len/row_length)
  
  #min_value = min(dplyr::select(df, !! label), na.rm = TRUE)
  #max_value = max(dplyr::select(df, !! label), na.rm = TRUE)
  
  # By default, use the 1% and 99% quantiles for low, high limits
  
  if (is.na(low)) {
    low = quantile(df %>% dplyr::select(!! label), probs = c(0.01), na.rm = TRUE)[[1]]
  }
  if (is.na(high)) {
    high = quantile(df %>% dplyr::select(!! label), probs = c(0.99), na.rm = TRUE)[[1]]
  }
  
  row_plots <- list()
  

  for (row in 0:(strips - 1)) {
    
    lower_limit = 1 + (row * row_length)
    upper_limit = min(sequence_len, row_length + (row * row_length))
    
    strip_sequence = str_split(substr(sequence, lower_limit, upper_limit), '')[[1]]
    
    # if this is the last iteration, allow the scale to be added
    
    if (row == (strips - 1)) {
      legend_position = "right"
    }
    
    # TODO: add size aes to tile outlines
    
    new_row = ggplot(data = df %>% filter(pos %in% c(lower_limit:upper_limit)) %>%
                       mutate(score = pmin(pmax(!! label, low, na.rm = FALSE), high, na.rm = FALSE)), 
                     aes(x = pos,
                         y = factor(variants, level = order_in, labels = names),
                         fill = score, height = 0.9, width = 0.9, color = 'black'), size = 3) +
      geom_tile(fill = NA, height = 0.9, width = 0.9, color = 'black', size = 0.4) +
      geom_tile(aes(size = as.factor(is.wt), color = as.factor(is.wt))) +
      scale_fill_continuous_divergingx(palette = 'RdBu',
                                       mid = 0,
                                       l1 = 0.2,
                                       l3 = 0.2,
                                       p1 = p1_in,
                                       p2 = p2_in,
                                       p3 = p3_in,
                                       p4 = p4_in,
                                       rev=invert_scale,
                                       limits = c(low, high),
                                       na.value = 'lightyellow') + 
      scale_x_continuous(breaks = seq(lower_limit, upper_limit, by = 5),
                         expand = c(0,0),
                         sec.axis = sec_axis(
                           trans = ~.,
                           name = "Sequence",
                           breaks = seq(lower_limit, upper_limit),
                           labels = strip_sequence,
                           guide = derive()
                         )) +
      coord_fixed(ratio = 1) +
      theme(
        panel.background = element_rect(fill = "lightyellow"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 4),
        axis.text = element_text(size = 4),
        axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
        axis.ticks = element_blank(),
        legend.position=legend_position
      ) +
      guides(color = FALSE,
             size = FALSE) + 
      scale_size_manual(values = c("TRUE" = 0.2, "FALSE" = 0.1)) +
      scale_color_manual(values = c("TRUE" = "lightgreen", "FALSE" = "black")) +
      labs(y = "Mutation", x = "Position")
    
    row_plots <- append(row_plots, list(new_row))
    
  }
  
  heatmap <- plot_grid(plotlist = as.list(row_plots), nrow = strips, ncol = 1)
  
  if (!is.null(output_file)) {
    ggsave(output_file, height = 2.5*strips, width = 8.5, heatmap)
  }
  
  if (print) {
    print(heatmap)
  }
  
}
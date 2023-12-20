library(MASS)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(tidyverse)
library(colorspace)

# Function for calculating rough cutoffs for significant mutation effect scores
# based on the distribution of synonymous scores. This assumes that the
# synonymous distribution is approximately gaussian, which might not hold.
# df is the dataframe of scores
# label is the column of scores in the dataframe

# This returns a pair of numbers: the upper and lower cutoffs, which are 
# determined by the mean +/- 2 standard deviations of the fit to the synonymous
# distribution. It also returns the mean and standard deviation of the fit.

calc_syn_cutoffs <- function(df, label)
{
  label <- enquo(label)
  
  fit <- fitdistr(
    df %>% filter(!is.na(!! label) & mutation_type == "S") %>% dplyr::select(!! label) %>% unlist(),
    densfun = "normal"
  )
  
  cutoff_1 <- fit$estimate[["mean"]] + 2 * fit$estimate[["sd"]]
  cutoff_2 <- fit$estimate[["mean"]] - 2 * fit$estimate[["sd"]]
  mean <- fit$estimate[["mean"]]
  sd <- fit$estimate[["sd"]]
  
  return( c(cutoff_1, cutoff_2, mean, sd) )
}


# This function plots a scatter plot of two scores, with cutoffs for significant
# mutation effects calculated from the synonymous distribution. It also plots
# the distribution of synonymous scores.

# df: dataframe of scores
# label_x: column of scores for x axis
# label_y: column of scores for y axis
# output_file: file to save plot to
# xlab: label for x axis
# ylab: label for y axis

plot_score_scatter <- function(df, label_x, label_y, output_file = NULL,
                               xlab = NULL, ylab = NULL, fixed = FALSE,
                               density = FALSE)
{
  
  label_x = enquo(label_x)
  label_y = enquo(label_y)
  
  # first, calculate upper and lower cutoffs for label_x and label_y
  x_syn_fit <- calc_syn_cutoffs(df, !!label_x)
  y_syn_fit <- calc_syn_cutoffs(df, !!label_y)
  
  # plot the scatter plot
  
  score_plot <- ggplot(df %>% filter(mutation_type != "X")) +
    geom_point(aes(y = !!label_y, x = !!label_x), alpha = 0.2) +
    geom_vline(xintercept = x_syn_fit[1], linetype = 2) +
    geom_vline(xintercept = x_syn_fit[2], linetype = 2) +
    geom_hline(yintercept = y_syn_fit[1], linetype = 2) +
    geom_hline(yintercept = y_syn_fit[2], linetype = 2) +
    theme_classic()
  
  if (density) {
    score_plot <- score_plot + geom_density_2d(data = df %>% filter(mutation_type %in% c("M", "S", "I", "N")),
                                               aes(y = !!label_y,
                                                   x = !!label_x, 
                                                   color = mutation_type),
                                               na.rm = TRUE,
                                               alpha = 0.8,
                                               adjust = 1,
                                               show.legend = FALSE,
                                               contour_var = "ndensity"
    ) + geom_density_2d(data = df %>% filter(mutation_type =="D"),
                                                 aes(y = !!label_y,
                                                     x = !!label_x, 
                                                     color = mutation_type),
                                              na.rm = TRUE,
                                              alpha = 0.8,
                                              adjust = 1,
                                              show.legend = FALSE,
                                              contour_var = "density",
                                              bins = 4
      ) + scale_color_manual(values = c("M" = "#009E73",
                                    "S" = "#CC79A7",
                                    "I" = "#56B4E9",
                                    "D" = "#D55E00",
                                    "N" = "grey"))
  }
  
  if (fixed) {
    score_plot <- score_plot + coord_fixed(ratio = 1, expand = FALSE)
  }
  
  if (!is.null(xlab)) {
    score_plot <- score_plot + xlab(xlab)
  }
  
  if (!is.null(ylab)) {
    score_plot <- score_plot + ylab(ylab)
  }
  
  # plot the first density plot
  
  dens1 <- ggplot(df %>% filter(mutation_type != "X"),
                  aes(x = !!label_y, fill = mutation_type)) + 
    geom_density(alpha = 0.4) + 
    geom_function(fun = dnorm, args = list(mean = y_syn_fit[3],
                                           sd = y_syn_fit[4]),
                  color = "red", linetype = 2) +
    geom_vline(xintercept = y_syn_fit[1], linetype = 2) +
    geom_vline(xintercept = y_syn_fit[2], linetype = 2) +
    scale_fill_manual(
      values = c("M" = "#009E73",
                 "S" = "#CC79A7",
                 "I" = "#56B4E9",
                 "D" = "#D55E00",
                 "N" = "grey")) +
    theme_void() + 
    coord_flip() +
    theme(legend.position = "none")
  
  # plot the second density plot
  
  dens2 <- ggplot(df %>% filter(mutation_type != "X"),
                  aes(x = !!label_x, fill = mutation_type)) + 
    geom_density(alpha = 0.4) + 
    geom_function(fun = dnorm, args = list(mean = x_syn_fit[3],
                                           sd = x_syn_fit[4]),
                  color = "red", linetype = 2) +
    geom_vline(xintercept = x_syn_fit[1], linetype = 2) +
    geom_vline(xintercept = x_syn_fit[2], linetype = 2) +
    scale_fill_manual(
      values = c("M" = "#009E73",
                 "S" = "#CC79A7",
                 "I" = "#56B4E9",
                 "D" = "#D55E00",
                 "N" = "grey")) +
    theme_void() + 
    theme(legend.position = "none")
  
  # combine the plots
  
  score_scatter <- dens2 + plot_spacer() + score_plot + dens1 + 
    plot_layout(ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
  
  # save the plot if desired
  if (!is.null(output_file)) {
    ggsave(output_file, width = 5, height = 5, score_scatter)
  }
  
  return(score_scatter)
}

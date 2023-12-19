library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(tidyverse)
library(colorspace)



#Color scheme:
  #D55E00 - Vermillion - deletions
  #56B4E9 - Sky blue - insertions
  #009E73 - bluish green - substitutions
  #CC79A7 - reddish purple - synonymous
 
# A color-blind palette with grey or black:

# The palette with grey:
cbp1 <- c("#D55E00", "#56B4E9", "#009E73", "#CC79A7",
          "#E69F00", "#F0E442", "#0072B2", "#999999")

# Levels for displaying mutations

default_order <- c("A", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L",
           "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
           "D_1", "D_2", "D_3",
           "I_1", "I_2", "I_3")

default_variant_names <- c("A", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L",
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
                   "Del x1", "Del x2", "Del x3",
                   "Ins x1 (G)", "Ins x2 (GS)", "Ins x3 (GSG)") 

plot_counts_df <- function(variant_df, experiment) {
  
  # Plot 1: coverage
  plot <- variant_df %>% filter(mutation_type != "S") %>% 
    ggplot(aes(x=pos, y=count, color=length)) + 
    geom_point() +
    facet_grid(rows = vars(mutation_type), scales = "free") +
    labs(x = "Position") + 
    theme_classic()
  
  print(plot)
  
  # Plot 2: counts averaged across positions, by type
  plot <- variant_df %>% filter(mutation_type != "S") %>% 
    ggplot(aes(x=factor(mutation, level = default_order),
               y=count)) + 
    geom_boxplot(outlier.colour=NA) +
    xlab("Variant") + 
    ylab("Counts") +
    labs(x = "Mutation") + 
    ggtitle(paste(title, " observed counts")) +
    theme_classic() 
  print(plot)
  
  # Plot 3: barchart per position
  plot <- variant_df %>% filter(mutation_type != "S") %>% 
    group_by(pos) %>% 
    summarise(count = mean(count),
              chunk = chunk,
              pos=pos,
              mutation_type=mutation_type,
              mutation=mutation,
              length=length) %>%
    ggplot(aes(x=pos, y=count, fill=mutation_type)) +
    geom_bar(stat="identity") +
    ggtitle(paste(title, " counts")) +
    theme_classic()
  
  print(plot)
  
  
  # Plot 4: heatmap
  plot <- variant_df %>% filter(mutation_type != "S") %>% 
    ggplot(aes(pos, y=factor(mutation, level = default_order), fill=count)) + 
    geom_tile() +
    scale_fill_gradient(high="blue", low="red", limits=c(1, 1400)) +
    xlab("Position") +
    ylab("Variant") +
    ggtitle(paste(title, " observed counts")) +
    theme_classic()
  
  print(plot)
  
  # Plot 5: chunk normalized position dot plot
  plot <- variant_df %>% filter(mutation_type != "S") %>% 
    ggplot(aes(chunk_pos, count)) +
    geom_point() +
    ggtitle(paste(title, " counts within chunk position")) +
    theme_classic()
  
  print(plot)
  
  # Plot 6: chunk normalized position bar plot
  plot <- variant_df %>% filter(mutation_type != "S") %>% 
    ggplot(aes(x=factor(chunk_pos), y=count, color=mutation_type)) +
    geom_boxplot() +
    ggtitle(paste(title, " counts within chunk position")) +
    theme_classic()
  
  print(plot)
  
  
  # Plot 7: Lorenz plot
  plot <- variant_df %>% 
    ggplot(aes(x = count, color = mutation_type)) + 
    stat_lorenz() +
    coord_fixed() + 
    geom_abline(linetype = "dashed") +
    hrbrthemes::scale_x_percent() +
    hrbrthemes::scale_y_percent() +
    hrbrthemes::theme_ipsum_rc() +
    theme(legend.title = element_blank()) +
    labs(x = "Cumulative Percentage of variants",
         y = "Cumulative Percentage of reads",
         title = paste(title, " library diversity")) +
    annotate_ineq(variant_df$count)
  print(plot)
  
}

process_df <- function(df, chunk, chunks) {
  
  # Creates an empty data frame to hold variants
  # TODO: deal with the case where the first row has some empty columns
  # example: OLS_2
  #variantCounts_colnames <- c("counts", "coverage", "mean_length", "length_NT", "NT", "length_codon", "codon", "AA", "mutations")
  
  variants_df <- as.data.frame(setNames(replicate(9,numeric(0), simplify = F),  c("count", "pos", "chunk_pos", "type", "name", "codon", "chunk", "mutation", "len")))
  
  for  (index in 1:nrow(df)) {
    # Loop through the variants discovered by ASM. Filtering strategy:
    #   -Reject any variants with multiple substitutions 
    #   -Reject any frameshifting mutations
    #   (designed mutations are singles, interpreting multi-codon indels as singles)
    #   -Accept G, GS, and GSG insertions
    #   -Accept 1,2,3x deletions
    #   -Further process "insdel" calls
    # For each one, we will add the observed counts
    # to the appropriate 
    # For synonymous variants, the "length_codon" is 0
    
    pos = -1
    chunk_pos = 0
    type = ""
    name = ""
    mutation = ""
    len = 0
    
    mutant = df[index, ]$mutations
    count = df[index, ]$counts
    codon = df[index, ]$codon
    AA = df[index, ]$AA
    
    if (is.na(mutant)) {
      next
    }
    
    mutants = str_split(mutant, ";", simplify=TRUE)
    pos = as.numeric(str_match(df[index, ]$mutations, "^[A-Z]([0-9]+)[dA-Z_]")[2])
    
    if (length(mutants) == 1) {
      
      if (grepl("insdel", mutant)) {
        del_pos = as.numeric(str_match(mutant, "^[A-Z]([0-9]+)_[A-Z]([0-9]+)")[2:3])
        del_len = as.numeric(del_pos[2]) -  as.numeric(del_pos[1]) + 1
        
        # Account for the case where deletion of the 
        # form XYZA->X--- (just deletion) is being called as
        # XYZA->---X (deletion + mutation)
        # In this case, mutation start is off by 1, and deletion length is as well.
        insdel_residues = str_match(mutant, "^([A-Z])[0-9]+_[A-Z][0-9]+insdel[A-Z]*([A-Z]$)")[2:3]
        
        
        if(is.na(insdel_residues[1])) { next }
        
        if (insdel_residues[1] == insdel_residues[2]) {
          del_pos =  as.numeric(del_pos) + 1
          del_len =  del_len - 1
        }
        
        pos = del_pos[1]
        
        if (del_len == 1) {
          type = "D"
          mutation = "D_1"
          len = 1
        } else if (del_len == 2) {
          type = "D"
          mutation = "D_2"
          len = 2
        } else if (del_len == 3) {
          type = "D"
          mutation = "D_3"
          len = 3
        } else {
          type = "XD"
        }
        
        # Done with insdel case 
      }
      
      
      else if (grepl("fs", mutant)) {
        next
      }
      
      else if (grepl("insG$", mutant)) {
        type = "I"
        mutation = "I_1"
        len = 1
      }
      
      else if (grepl("insGS$", mutant)) {
        type = "I"
        mutation = "I_2"
        len = 2
        
      }
      
      else if (grepl("insGSG$", mutant)) {
        type = "I"
        mutation = "I_3"
        len = 3
      }
      
      # Deal with deletions
      
      else if (grepl("del$", mutant)) {
        
        if (grepl("^[A-Z]([0-9]+)_[A-Z]([0-9]+)del$", mutant)) {
          
          del_pos = str_match(mutant, "^[A-Z]([0-9]+)_[A-Z]([0-9]+)")[2:3]
          del_len = as.numeric(del_pos[2]) -  as.numeric(del_pos[1]) + 1
          
          if (del_len == 2) {
            type = "D"
            mutation = "D_2"
            len = 2
          }
          else if (del_len == 3) {
            type = "D"
            mutation = "D_3"
            len = 3
          }
          else {next}
        }
        
        else {
          type = "D"
          mutation = "D_1"
          len = 1
        }
        
      }
      
      # Identify synonymous mutants
      # AA is the field with AA mutation: we want to call ones with prefix S:
      else if (codon == "") {
        split_AA = strsplit(AA, ", ")[[1]]
        
        if (length(split_AA == 1)) {
          if (substr(split_AA[1], 0, 1) == "S") {
            type = "S"
            pos = as.numeric(substr(codon, 0, 1))
            len = 0
            name = codon
          }
        }
      }
      
      else {
        type = "M"
        mutation = substring(mutant,nchar(mutant))
        len = 1
      }
      
      # Don't add variants that didn't match any case
      
      if (!is.na(pos)) {
        for (p in chunks) {
          if (pos > p[1] & pos < p[2]) {
            chunk_pos = as.numeric(pos - p[1] + 1)
            chunk = round(p[2]/(p[2]-p[1]))
          }
        }
        variants_df[nrow(variants_df)+1,] <- c(count, pos, chunk_pos, type, mutant, codon, chunk, mutation, len)
      }
      
    }
  }
  
  variants_df$count <- as.integer(variants_df$count)
  variants_df$pos <- as.integer(variants_df$pos)
  variants_df$chunk_pos <- as.integer(variants_df$chunk_pos)
  variants_df$chunk <- as.integer(variants_df$chunk)
  
  # Identify the position of the mutation within the relevant chunk
  
  return(variants_df %>% 
           group_by(name) %>% 
           summarise(counts = sum(count), chunk_pos = first(chunk_pos), chunk = first(chunk), pos=first(pos), type=first(type), mutation=first(mutation), len=first(len)))
}


# Parses an HGVS string, outputs a vector with variant info

parse_hgvs <- function(hgvs_string) {
  variant = ""
  pos = -1
  len = -1
  mutation_type = ""
  
  # WT case, as in enrich2 format
  if (str_detect(hgvs_string, "_wt")) {
    variant = "Z"
    pos = -1
    len = -1
    mutation_type = "X"
  }
  
  # M/S/N
  if (str_detect(hgvs_string, "[A-Z][0-9]+[A-Z]+")) {
    len = 1 
    match = str_match(hgvs_string, "([A-Z])([0-9]+)([A-Z]+)")
    pos = match[3]
    if (match[2] == match[4]) {
      mutation_type = "S"
      variant = match[4]
    } else if (match[4] == 'X') {
      mutation_type = "N"
      variant = match[4]
    } else {
      mutation_type = "M"
      variant = match[4]
    }
  }
  
  # D
  if (str_detect(hgvs_string, ".*del")) {
    mutation_type = "D"
    if (str_detect(hgvs_string, "[A-Z][0-9]+_[A-Z][0-9]+del")) {
      # D_2, D_3
      match = str_match(hgvs_string, "[A-Z]([0-9]+)_[A-Z]([0-9]+)del")
      pos = match[2]
      len = strtoi(match[3]) - strtoi(match[2]) + 1
      variant = paste("D_", as.character(len), sep="")
    } else {
      # D_1
      len = 1
      match = str_match(hgvs_string, "([A-Z])([0-9]+)del")
      pos = match[3]
      variant = "D_1"
    }
  }
  
  # I
  if (str_detect(hgvs_string, ".*ins.*")) {
    mutation_type = "I"
    match = str_match(hgvs_string, "[A-Z]([0-9]+)_[A-Z][0-9]+ins([A-Z]+)")
    len = nchar(match[3])
    pos = match[2]
    variant = paste("I_", as.character(len), sep="")
  }
  return(c(variant, pos, len, mutation_type))
}

## Turn a 3AA HGVS string into a 1AA HGVS string
# Input strings are in form "p.Lys116Lys"
# Output are in form "p.(K116K)"

convert_3AA_hgvs <- function(hgvs_string) {
  if (str_detect(hgvs_string, "_wt")) {
    hgvs_string_1x = "_wt"
  }
  if (str_detect(hgvs_string, "p\\.([A-Za-z]{3})([0-9]+)([A-Za-z]{3})")) {
    match = str_match(hgvs_string, "p\\.([A-Za-z]{3})([0-9]+)([A-Za-z]{3})")
    
    AA1 = toupper(match[2])
    AA2 = toupper(match[4])
    
    AA1 = aa.table[AA1, ]['aa1'][[1]]
    AA2 = aa.table[AA2, ]['aa1'][[1]]
    hgvs_string_1x = paste( 'p.(', AA1, match[3], AA2, ')', sep="")
    
  }
  
  return(hgvs_string_1x)
}

# Parses an HGVS string with 3x AA names, outputs a vector with variant info
# strings are in form "p.Lys116Lys"
# coerce to "p.(K116K)" and call parse_hgvs on it

parse_hgvs_2 <- function(hgvs_string) {

  if (str_detect(hgvs_string, "_wt")) {
    hgvs_string_1x = "_wt"
  }
  
  
  #Coerce the HGVS string.
  
  if (str_detect(hgvs_string, "p\\.([A-Za-z]{3})([0-9]+)([A-Za-z]{3})")) {
    match = str_match(hgvs_string, "p\\.([A-Za-z]{3})([0-9]+)([A-Za-z]{3})")

  AA1 = toupper(match[2])
  AA2 = toupper(match[4])
  
  AA1 = aa.table[AA1, ]['aa1'][[1]]
  AA2 = aa.table[AA2, ]['aa1'][[1]]
  hgvs_string_1x = paste( 'p.(', AA1, match[3], AA2, ')', sep="")

  }
  
  return(parse_hgvs(hgvs_string_1x))
}

# Takes a df with identifiers in HGVS format, outputs
# a df with four additional columns with variant info
process_hgvs_df <- function(df) {
  
  n = nrow(df)
  
  variant_list = setNames(
    as_tibble(
      data.frame(
        rep(NA, n)
      )), 
    c("variants"))
  pos_list = setNames(
    as_tibble(
      data.frame(
        rep(NA, n)
      )), 
    c("pos"))
  len_list = setNames(
    as_tibble(
      data.frame(
        rep(NA, n)
      )), 
    c("len"))
  mutation_type_list = setNames(
    as_tibble(
      data.frame(
        rep(NA, n)
      )), 
    c("mutation_type"))
  
  for(i in 1:n) {
    output = parse_hgvs(df$hgvs[i])
    variant_list[i,] = output[1]
    pos_list[i,] = strtoi(output[2])
    len_list[i,] = output[3]
    mutation_type_list[i,] = output[4]
  }
  
  
  df <- add_column(df, pos_list)
  df <- add_column(df, len_list)
  df <- add_column(df, mutation_type_list)
  df <- add_column(df, variant_list)
  
  return(df)
}

# Takes a df with identifiers in HGVS format with 3-character AA name, outputs
# a df with four additional columns with variant info.
# Note: this is for VatA, which only has substitutions and no indels. Length
# is dropped, and also the positions are already correctly set.

process_hgvs_df_2 <- function(df) {
  
  n = nrow(df)
  
  variant_list = setNames(
    as_tibble(
      data.frame(
        rep(NA, n)
      )), 
    c("variants"))
  pos_list = setNames(
    as_tibble(
      data.frame(
        rep(NA, n)
      )), 
    c("pos"))
  len_list = setNames(
    as_tibble(
      data.frame(
        rep(NA, n)
      )), 
    c("len"))
  mutation_type_list = setNames(
    as_tibble(
      data.frame(
        rep(NA, n)
      )), 
    c("mutation_type"))
  
  for(i in 1:n) {
    output = parse_hgvs_2(df$hgvs[i])
    variant_list[i,] = output[1]
    pos_list[i,] = strtoi(output[2])
    len_list[i,] = output[3]
    mutation_type_list[i,] = output[4]
  }
  
  
  df <- add_column(df, mutation_type_list)
  df <- add_column(df, variant_list)
  
  return(df)
}

# TODO: seq_len emits a warning about length.out:
# "Warning: first element used of 'length.out' argument"

map_scores_pdb <- function(input_pdb, mapping_scores, field, selection = NULL) {
  
  if (is.null(selection)) {
    
    selection = atom.select(input_pdb, "protein")
  }
  
  output_pdb = trim.pdb(input_pdb, selection)
  
  for (i in seq_len(dim(output_pdb$atom[1]))) {
    
    if (output_pdb$atom[i,]$resno > 0) {
      
      n = as.character(output_pdb$atom[i,]$resno)
      j = which(mapping_scores['pos'] == n)

      if (length(j) == 0) {
        score = 0
        
      } else {
        score = mapping_scores[j, field][[1]]
      }
      
      if (!is.na(score)) {
        
        output_pdb$atom[i,]$b = score
        
      } else {
        
        output_pdb$atom[i,]$b = 0
        
      }
    } else {
      output_pdb$atom[i,]$b = 0
    }
  }
  
  return(output_pdb)
}

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

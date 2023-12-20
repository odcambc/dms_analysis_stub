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
# Alphabetical order

alphabetical_order <- c("A", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L",
           "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
           "X",
           "D_1", "D_2", "D_3",
           "I_1", "I_2", "I_3")

alphabetical_variant_names <- c("A", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L",
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
                   "*",
                   "Del x1", "Del x2", "Del x3",
                   "Ins x1 (G)", "Ins x2 (GS)", "Ins x3 (GSG)") 

# Ordering by physical chemistry
physical_order <- c("G", "A", "I", "L", "V",
                    "C", "M",
                    "S", "T", 
                    "K", "R", "H",
                    "D", "E", "N", "Q",
                    "F", "W", "Y",
                    "P",
                    "X",
                    "D_1", "D_2", "D_3",
                    "I_1", "I_2", "I_3")

physical_variant_names <- c("G", "A", "I", "L", "V",
                    "C", "M",
                    "S", "T", 
                    "K", "R", "H",
                    "D", "E", "N", "Q",
                    "F", "W", "Y",
                    "P",
                    "*",
                    "Del x1", "Del x2", "Del x3",
                    "Ins x1 (G)", "Ins x2 (GS)", "Ins x3 (GSG)")

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
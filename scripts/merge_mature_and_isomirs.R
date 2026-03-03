#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3){
  stop("Usage: Rscript merge_miRNA_counts.R <isomir_counts.csv> <mature_mirbase.readcounts.tsv> <output.csv>")
}

isomir_file <- args[1]
mature_file <- args[2]
output_file <- args[3]

# ----------------------------
# Read isomiR file
# ----------------------------
isomir <- read_csv(isomir_file, show_col_types = FALSE)

# ----------------------------
# Read mature file
# ----------------------------
mature <- read_tsv(mature_file, show_col_types = FALSE) %>%
  select(-Length) %>%
  rename_with(~ str_remove(., "\\.mature$")) %>%
  rename(miRNA = 1)

# ----------------------------
# Align rows by miRNA
# ----------------------------
all_mirnas <- union(isomir$miRNA, mature$miRNA)

isomir_aligned <- isomir %>%
  complete(miRNA = all_mirnas) %>%
  replace(is.na(.), 0)

mature_aligned <- mature %>%
  complete(miRNA = all_mirnas) %>%
  replace(is.na(.), 0)

# Ensure same row order
isomir_aligned <- isomir_aligned %>% arrange(miRNA)
mature_aligned <- mature_aligned %>% arrange(miRNA)

# ----------------------------
# Add count matrices
# ----------------------------
count_cols <- setdiff(colnames(isomir_aligned), "miRNA")

summed_matrix <- isomir_aligned %>%
  mutate(across(all_of(count_cols),
                ~ . + mature_aligned[[cur_column()]]))

summed_matrix <- summed_matrix[-1,]

# ----------------------------
# Write output
# ----------------------------
write_csv(summed_matrix, output_file)
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: pivot_isomirs_longer.R <input.tsv> <output.csv>")
}

input_file  <- args[1]
output_file <- args[2]

# ---- read input ----
df <- readr::read_tsv(
  input_file,
  col_types = cols(.default = col_character())
)

# ---- clean sample column names ----
# Strip the suffix .seqcluster.hairpin.aln.srt
metadata_cols <- c("UID", "Read", "miRNA", "Variant", "iso_5p", "iso_3p", "iso_add3p", "iso_snp")

df <- df %>%
  rename_with(
    ~ str_remove(.x, "\\.seqcluster\\.hairpin\\.aln\\.srt$"),
    -all_of(metadata_cols)
  )

# ---- pivot longer ----
long_df <- df %>%
  pivot_longer(
    cols      = -all_of(metadata_cols),
    names_to  = "Sample_ID",
    values_to = "Counts"
  )

# ---- write output ----
readr::write_csv(long_df, output_file)
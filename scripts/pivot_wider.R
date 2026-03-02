#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# ---- command line args ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: pivot_wider_and_sum_mirna.R <input_long.csv> <output_wide.csv> <output_mirna_sum.csv>")
}

input_file        <- args[1]
output_wide_file  <- args[2]
output_sum_file   <- args[3]

# ---- read input ----
long_data <- readr::read_csv(
  input_file,
  col_types = cols(
    Counts = col_double(),
    .default = col_character()
  )
)

# ---- 1. pivot to wide (variant-level) ----
wide_variant <- long_data %>%
  tidyr::pivot_wider(
    names_from  = Sample_ID,
    values_from = Counts,
    values_fill = 0
  )

readr::write_csv(wide_variant, output_wide_file)

# ---- 2. sum variants per miRNA ----
mirna_summed <- long_data %>%
  group_by(miRNA, Sample_ID) %>%
  summarise(
    Counts = sum(Counts, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = Sample_ID,
    values_from = Counts,
    values_fill = 0
  ) %>%
  arrange(miRNA)

readr::write_csv(mirna_summed, output_sum_file)
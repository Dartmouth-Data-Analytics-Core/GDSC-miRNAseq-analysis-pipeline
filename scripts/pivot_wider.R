library(tidyverse)

#----- Set command line args
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: script.R <input_file> <outputDir>")
}

input_file <- args[1] # Long format isomiR master counts
outputDir <- args[2] # isomiR count dir

if (!dir.exists(outputDir)) {
  dir.create(outputDir)
}

#----- Read in the input file
data <- read.csv(input_file, na.strings = "")

#----- Collapse data by miRNA family
sum_by_miRNA <- function(df) {
  df %>%
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
}

#----- Run function and save
summed <- sum_by_miRNA(data)
write.csv(summed, file = paste0(outputDir, "raw_merged_canonical_and_all_isomirs.csv"), row.names = FALSE, quote = FALSE)

#----- Get just the canonical counts
canonical <- data[data$Variant == "NA",]
write.csv(canonical, file = paste0(outputDir, "raw_canonical_counts.csv"), row.names = FALSE, quote = FALSE)

#----- Get noncanonical
isomirs <- data[data$Variant != "NA",]

#----- Assign isomir class
assign_iso_class <- function(df) {
  df %>%
    mutate(
      iso_class = case_when(
        as.numeric(iso_5p) != 0 ~ "5p",
        as.numeric(iso_3p) != 0 ~ "3p",
        as.numeric(iso_add3p) != 0 ~ "non-templated",
        as.numeric(iso_snp) != 0 ~ "snp",
        TRUE ~ "canonical"
      )
    )
}

isomirs <- assign_iso_class(isomirs)

#----- Output all isomir data
write.csv(isomirs, file = paste0(outputDir, "raw_all_isomir_counts.csv"), row.names = FALSE, quote = FALSE)

#----- Split by each category
classes <- unique(isomirs$iso_class)
for (i in classes) {
  isoSub <- isomirs[isomirs$iso_class == i,]
  fileName <- paste0(i, "_isomir_counts.csv")
  write.csv(isoSub, file = paste0(outputDir, fileName), row.names = FALSE, quote = FALSE)
}


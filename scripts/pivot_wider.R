suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

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
#write.csv(canonical, file = paste0(outputDir, "raw_canonical_counts.csv"), row.names = FALSE, quote = FALSE)
readr::write_tsv(
  canonical,
  file = paste0(outputDir, "raw_canonical_counts.tsv")
)

#----- Get noncanonical
isomirs <- data[data$Variant != "NA",]

#----- Assign isomir class
assign_iso_class <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      iso_class = paste(
        c(
          if (as.numeric(iso_5p) != 0) "5p",
          if (as.numeric(iso_3p) != 0) "3p",
          if (as.numeric(iso_add3p) != 0) "non-templated",
          if (as.numeric(iso_snp) != 0) "snp"
        ),
        collapse = ","
      ),
      iso_class = ifelse(iso_class == "", "canonical", iso_class)
    ) %>%
    ungroup()
}
isomirs <- assign_iso_class(isomirs)

#----- Output all isomir data
#write.csv(isomirs, file = paste0(outputDir, "raw_all_isomir_counts.csv"), row.names = FALSE, quote = FALSE)
readr::write_tsv(
  canonical,
  file = paste0(outputDir, "raw_all_isomir_counts.tsv")
)

classes <- c("5p", "3p", "non-templated", "snp", "canonical")
for (i in classes) {
  
  isoSub <- isomirs %>%
    filter(str_detect(iso_class, paste0("\\b", i, "\\b")))
  
  fileName <- paste0(i, "_isomir_counts.tsv")
  
  write_tsv(
    isoSub,
    file = paste0(outputDir, fileName))
}
   


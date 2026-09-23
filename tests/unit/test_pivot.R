#!/usr/bin/env Rscript
# Unit tests for scripts/pivot_longer.R and scripts/pivot_wider.R
# Uses base R stopifnot() for lightweight assertions (no testthat dependency).

# Resolve the repo root from the --file= argument passed by Rscript
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
SCRIPTS <- if (length(file_arg)) {
  file.path(
    dirname(dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))),
    "scripts"
  )
} else {
  file.path(getwd(), "scripts")
}
if (!dir.exists(SCRIPTS)) SCRIPTS <- file.path(getwd(), "scripts")

tmp <- tempdir()

# ── Helper ────────────────────────────────────────────────────────────────────

write_isomir_tsv <- function(path) {
  content <- paste(
    "UID\tRead\tmiRNA\tVariant\tiso_5p\tiso_3p\tiso_add3p\tiso_snp\ts1.seqcluster.hairpin.aln.srt\ts2.seqcluster.hairpin.aln.srt",
    "uid1\tTAGCTT\tdre-miR-21\tNA\t0\t0\t0\t0\t120\t95",
    "uid2\tTGAGGT\tdre-let-7a\tNA\t0\t0\t0\t0\t340\t410",
    "uid3\tTGAGGT\tdre-let-7a\tiso_3p:1\t0\t1\t0\t0\t15\t8",
    sep = "\n"
  )
  writeLines(content, path)
}

# Source a script in-process with mocked commandArgs to avoid subprocess overhead
run_script <- function(script_path, args_vec) {
  commandArgs <- function(trailingOnly = FALSE) if (trailingOnly) args_vec else c("Rscript", args_vec)
  source(script_path, local = TRUE)
}

# ── Test pivot_longer.R ───────────────────────────────────────────────────────

cat("-- Test: pivot_longer.R produces long-format CSV --\n")
tsv_in <- file.path(tmp, "test_isomir.tsv")
csv_out <- file.path(tmp, "test_isomir_long.csv")
write_isomir_tsv(tsv_in)

source_script <- file.path(SCRIPTS, "pivot_longer.R")
if (!file.exists(source_script)) stop("pivot_longer.R not found at: ", source_script)

run_script(source_script, c(tsv_in, csv_out))

stopifnot("Output CSV not created" = file.exists(csv_out))

long_df <- read.csv(csv_out)

stopifnot("Sample_ID column missing" = "Sample_ID" %in% names(long_df))
stopifnot("Counts column missing" = "Counts" %in% names(long_df))
stopifnot("Wrong number of rows in long format" = nrow(long_df) == 6)
stopifnot("Suffix not stripped from sample names" =
  all(long_df$Sample_ID %in% c("s1", "s2")))

cat("PASS: pivot_longer.R\n\n")

# ── Test pivot_wider.R ────────────────────────────────────────────────────────

cat("-- Test: pivot_wider.R produces wide miRNA count table --\n")

source_script_wider <- file.path(SCRIPTS, "pivot_wider.R")
if (!file.exists(source_script_wider)) stop("pivot_wider.R not found at: ", source_script_wider)

out_dir <- file.path(tmp, "quant_out")
dir.create(out_dir, showWarnings = FALSE)

run_script(source_script_wider, c(csv_out, paste0(out_dir, "/")))

merged_csv <- file.path(out_dir, "raw_merged_canonical_and_all_isomirs.csv")
stopifnot("raw_merged_canonical_and_all_isomirs.csv not created" = file.exists(merged_csv))

wide_df <- read.csv(merged_csv)
stopifnot("miRNA column missing in wide table" = "miRNA" %in% names(wide_df))
stopifnot("Sample columns missing in wide table" =
  all(c("s1", "s2") %in% names(wide_df)))

row_21 <- wide_df[wide_df$miRNA == "dre-miR-21", ]
stopifnot("dre-miR-21 row missing" = nrow(row_21) == 1)
stopifnot("dre-miR-21 s1 count mismatch" = row_21$s1 == 120)

cat("PASS: pivot_wider.R\n\n")
cat("All R pivot tests passed.\n")

#!/usr/bin/env Rscript
# Unit tests for scripts/pivot_longer.R and scripts/pivot_wider.R
# Uses base R stopifnot() for lightweight assertions (no testthat dependency).

SCRIPTS <- file.path(dirname(dirname(dirname(normalizePath(sys.frame(1)$ofile,
  mustWork = FALSE)))), "scripts")
if (!nchar(SCRIPTS) || !dir.exists(SCRIPTS)) {
  SCRIPTS <- file.path(getwd(), "scripts")
}

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

cat("-- Test: pivot_longer.R produces long-format CSV --\n")
tsv_in <- file.path(tmp, "test_isomir.tsv")
csv_out <- file.path(tmp, "test_isomir_long.csv")
write_isomir_tsv(tsv_in)

source_script <- file.path(SCRIPTS, "pivot_longer.R")
if (!file.exists(source_script)) stop("pivot_longer.R not found at: ", source_script)

# Run via Rscript subprocess to avoid side effects
exit_code <- system2("Rscript", args = c(source_script, tsv_in, csv_out),
                     stdout = FALSE, stderr = FALSE)
stopifnot("pivot_longer.R exited with non-zero code" = exit_code == 0)
stopifnot("Output CSV not created" = file.exists(csv_out))

long_df <- read.csv(csv_out)

# Should have a Sample_ID column from pivoting
stopifnot("Sample_ID column missing" = "Sample_ID" %in% names(long_df))
stopifnot("Counts column missing" = "Counts" %in% names(long_df))

# 3 rows × 2 samples = 6 rows in long format
stopifnot("Wrong number of rows in long format" = nrow(long_df) == 6)

# Sample names should be stripped of the .seqcluster.hairpin.aln.srt suffix
stopifnot("Suffix not stripped from sample names" =
  all(long_df$Sample_ID %in% c("s1", "s2")))

cat("PASS: pivot_longer.R\n\n")

# ── Test pivot_wider.R ────────────────────────────────────────────────────────

cat("-- Test: pivot_wider.R produces wide miRNA count table --\n")

source_script_wider <- file.path(SCRIPTS, "pivot_wider.R")
if (!file.exists(source_script_wider)) stop("pivot_wider.R not found at: ", source_script_wider)

out_dir <- file.path(tmp, "quant_out")
dir.create(out_dir, showWarnings = FALSE)

exit_code2 <- system2("Rscript", args = c(source_script_wider, csv_out, paste0(out_dir, "/")),
                      stdout = FALSE, stderr = FALSE)
stopifnot("pivot_wider.R exited with non-zero code" = exit_code2 == 0)

merged_csv <- file.path(out_dir, "raw_merged_canonical_and_all_isomirs.csv")
stopifnot("raw_merged_canonical_and_all_isomirs.csv not created" = file.exists(merged_csv))

wide_df <- read.csv(merged_csv)
stopifnot("miRNA column missing in wide table" = "miRNA" %in% names(wide_df))
stopifnot("Sample columns missing in wide table" =
  all(c("s1", "s2") %in% names(wide_df)))

# dre-miR-21 canonical count for s1 should be 120
row_21 <- wide_df[wide_df$miRNA == "dre-miR-21", ]
stopifnot("dre-miR-21 row missing" = nrow(row_21) == 1)
stopifnot("dre-miR-21 s1 count mismatch" = row_21$s1 == 120)

cat("PASS: pivot_wider.R\n\n")
cat("All R pivot tests passed.\n")

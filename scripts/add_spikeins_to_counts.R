library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

mirna_file <- args[which(args == "--mirna") + 1]
spike_file <- args[which(args == "--spike") + 1]
out_file <- args[which(args =="--out") + 1]

# load miRNA counts
mir <- read_tsv(mirna_file)

# load spike-in stats (same parsing strategy used pela TAmiRNA)
spike <- read_delim(spike_file, "\t", skip = 3, col_types = cols_only("#Name" = col_character(), Reads = col_double())) %>%
  rename(mirbase_ID = "#Name", reads = Reads)

# calculate global lib size = total miRNA + spikeins
lib_size <- read_delim(spike_file, "\t", skip = 1, n_max = 1, col_names = FALSE)[[2]] %>%
  as.numeric()

# add normalized RPM
mir <- mir %>%
  mutate(rpm.lib = reads / lib_size * 1e6)

spike <- spike %>%
  mutate(rpm.lib = reads / lib_size * 1e6)

# combine
combined <- bind_rows(mir, spike)

write_tsv(combined, out_file)


library(tidyverse)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)
input_file <- args$'--mirbase'
spike_stats <- args$'--spike_stats'
spike_info <- args$'--spike_info'
sample_id <- args$'--sample'
out_slope <- args$'--out_slope'
out_norm <- args$'--out_norm'

data <- read_tsv(input_file)

spikes_obs <- read_delim(spike_stats, "\t", skip=3,
                         col_types = cols_only('#Name'=col_character(), Reads=col_double())) %>%
  rename(mirbase_ID = '#Name', reads = Reads)


spikes_info <- read_tsv(spike_info, show_col_types = FALSE) %>%
  rename(amol = amol)


df_spike <- spikes_obs %>%
  inner_join(spike_info, by = "mirbase_ID") %>%
  filter(reads > 0, amol > 0) %>%
  mutate(log_reads = log10(reads),
         log_amol  = log10(amol))

LOD <- min(df_spike$amol)


fit <- lm(log_reads ~ 0 + log_amol, data = df_spike)

coef(fit)

slope <- coef(fit)[2]
intercept <- coef(fit)[1]

r <- cor(df_spike$log_amol, df_spike$log_reads, method = "pearson")

write_tsv(tibble(slope = slope, intercept = intercept, LOD_amol_concentration = LOD, R_Pearson = r), out_slope)

mirna_norm <- data %>%
  mutate(amol_est = 10^(log10(reads) / slope))


write_tsv(mirna_norm, out_norm)
# -------------------------------------------------------------
# Spike-in QC and miRNA normalization pipeline
# -------------------------------------------------------------
# This script performs spike-in–based calibration and estimates
# miRNA concentrations from read counts.
#
# UNIT NOTE:
# All concentrations in this pipeline are expressed as molecules/µL.
# The calibration model predicts concentration in molecules/µL.
#
# -------------------------------------------------------------
# This script:
# 1. Loads spike-in and miRNA count tables
# 2. Annotates spike-in concentrations (molecules/µL)
# 3. Fits per-sample calibration models (molecules/µL ~ counts)
# 4. Computes QC metrics (R², detection limits, spike-ins detected)
# 5. Identifies miRNAs within the calibrated dynamic range
# 6. Normalizes miRNA counts to molecules/µL using the spike-in model
# 7. Generates per-sample QC plots
# 8. Outputs TSV files with QC metrics and normalized concentrations
# -------------------------------------------------------------


library(tidyverse)
library(ggplot2)
library(gridExtra)

if (!dir.exists("spikein_metrics/qc_plots")) {
  dir.create("spikein_metrics/qc_plots")
}

# -------------------------------------------------------------
# Load input data
# -------------------------------------------------------------
# Spike-in raw counts
spikes_obs <- read_tsv("spikein_counts/spikein.readcounts.tsv")
# miRNA raw counts (aligned to miRBase)
mirna_counts <- read_tsv("mirbase_counts/mirbase.readcounts.tsv")

# -------------------------------------------------------------
# Spike-in concentration metadata
# -------------------------------------------------------------
# Known amol concentrations of miND spike-ins (from provider)
spikes_info <- tibble(
  spikein_ID = c("miND-01", "miND-02", "miND-03", "miND-04", "miND-05", "miND-06", "miND-07"),
  amol_concentration = c(20, 5, 1.25, 0.3125, 0.075, 0.01, 0.005))

# Volumes used in spike-in mixture
amol <- 602214 # Avogadro constant used in provider's conversion (kept for clarity)
spikeInsVolume <- 1  # Spike-in input volume (µL)
finalVolume <- 5     # Final sample volume (µL)

# Adjust final spike-in concentration after dilution
spikes_info <- spikes_info %>%
  mutate(molecules_concentration = amol_concentration * spikeInsVolume * amol / finalVolume)


# -------------------------------------------------------------
# Reshape spike-in table and merge concentration info
# -------------------------------------------------------------
spikes_with_conc <- spikes_obs %>%
  pivot_longer(cols = -c(spikein_ID), names_to = "samples", values_to = "counts") %>%
  left_join(spikes_info, by = "spikein_ID")


# Spike-ins detected (counts > 0)
spikes_detected <- spikes_with_conc %>%
  filter(counts > 0)

# -------------------------------------------------------------
# Per-sample spike-in QC and calibration model
# -------------------------------------------------------------
# For each sample:
# - Check number of missing spike-ins
# - Fit linear model molecules_concentration ~ 0 + counts
# - Compute R², lower/upper detection limits, slope
spikeins_stats <- spikes_with_conc %>%
  group_by(samples) %>%
  do({
    df <- .
    # If ≥2 spike-ins missing, fail QC immediately
    spikein_missing <- df %>% filter(counts == 0)
    
    if (nrow(spikein_missing) >= 2) {
      tibble(
        sample_id = unique(df$samples),
        spikeins_detected = sum(df$counts > 0),
        mirnas_in_range = NA,
        spikein_lower_limit = NA,
        spikein_upper_limit = NA,
        intercept = NA,
        slope = NA,
        rsq = NA,
        qc = "FAILED (2 or more spike-ins missing)",
        spikeins_with_prediction = list(NULL)
      )
    } else {
      # Fit model using detected spike-ins only
      df2 <- df %>% filter(counts > 0)
      fit <- lm(molecules_concentration ~ 0 + counts, data = df2)
      # Prediction intervals for spike-ins
      pred_int <- suppressWarnings(predict(fit, interval = "prediction")) %>% as_tibble()
      df_with_pred <- df2 %>% bind_cols(pred_int)
      
      # Detection limits = smallest and largest observed spike-in counts
      lower <- min(df_with_pred$counts)
      upper <- max(df_with_pred$counts)
      rsq_value <- summary(fit)$r.squared
      
      tibble(
        sample_id = unique(df$samples),
        spikeins_detected = nrow(df2),
        mirnas_in_range = NA,
        spikein_lower_limit = lower,
        spikein_upper_limit = upper,
        intercept = 0,
        slope = coef(fit)[["counts"]],
        rsq = rsq_value,
        qc = ifelse(rsq_value < 0.95, "FAILED (R² < 0.95)", "OK"),
        spikeins_with_prediction = list(df_with_pred)
      )
    }
  }) %>%
  ungroup()

# -------------------------------------------------------------
# Save simple version of models for later use
# -------------------------------------------------------------
per_sample_models <- spikes_detected %>%
  group_split(samples) %>%
  map_df(~ {
    df <- .x
    fm <- lm(molecules_concentration ~ 0 + counts, data = df)
    tibble(
      sample_id = unique(df$samples),
      model = list(fm)
    )
  })

# -------------------------------------------------------------
# Annotate miRNAs with spike-in-derived detection range
# -------------------------------------------------------------
mirna_long <- mirna_counts %>%
  pivot_longer(
    cols = -c(mirbase_ID, Length),
    names_to = "sample_id",
    values_to = "reads"
  ) %>%
  filter(reads > 0) %>%
  left_join(
    spikeins_stats %>% select(sample_id, spikein_lower_limit, spikein_upper_limit),
    by = "sample_id"
  ) %>%
  mutate(
    in_range = reads >= spikein_lower_limit & reads <= spikein_upper_limit
  )

# Compute % of miRNAs within range per sample
mirnas_in_range_count <- mirna_long %>%
  group_by(sample_id) %>%
  summarise(
    mirnas_in_range = sum(in_range, na.rm = TRUE),
    percent_in_range = round((sum(in_range, na.rm = TRUE) / n()) * 100, 2),
    .groups = "drop"
  )

# Add to spike-in QC table
spikeins_stats <- spikeins_stats %>%
  select(-mirnas_in_range) %>%
  left_join(mirnas_in_range_count, by = "sample_id") %>%
  select(
    sample_id, spikeins_detected, spikein_lower_limit, spikein_upper_limit,
    mirnas_in_range, percent_in_range, slope, rsq, qc
  )

# -------------------------------------------------------------
# Compute Limit of Detection per sample
# -------------------------------------------------------------
LOD_by_sample <- spikes_detected %>%
  group_by(samples) %>%
  summarise(LOD = min(molecules_concentration), .groups = "drop") %>%
  rename(sample_id = samples)

# -------------------------------------------------------------
# Export spike-in QC metrics table
# -------------------------------------------------------------
spikeins_stats %>%
  left_join(LOD_by_sample, by = "sample_id") %>%
  rename(
    "Scale Factor/ Slope" = slope,
    "Limit of Detection molecules/uL (microliter)" = LOD,
    "Percent in range %" = percent_in_range,
    "Number of mirnas in Range" = mirnas_in_range,
    "Number of spike-ins detected" = spikeins_detected,
    "spike-in lower limit detection" = spikein_lower_limit,
    "spike-in upper limit detection" = spikein_upper_limit
  ) %>%
  select(-c(qc)) %>%
  write_tsv("spikein_metrics/spikein_detection_metrics.tsv")

# -------------------------------------------------------------
# Normalize miRNA counts into molecules/microliter using spike-in linear model
# -------------------------------------------------------------
mirna_normalized <- mirna_long %>%
  left_join(per_sample_models, by = "sample_id") %>%
  rowwise() %>%
  mutate(
    pred_data = list(
      predict(model,
              newdata = data.frame(counts = reads),
              interval = "prediction") %>% 
        as_tibble()
    )
  ) %>%
  unnest(pred_data) %>%
  mutate(
    molecules_concentration = round(fit, 2),
    interval = round(fit - lwr, 2)
  ) %>%
  ungroup() %>%
  select(mirbase_ID, Length, sample_id, reads, molecules_concentration)

# Export normalized miRNA table (wide format)
mirna_normalized %>%
  pivot_wider(
    names_from = sample_id,
    values_from = molecules_concentration,
    id_cols = c(mirbase_ID, Length),
    values_fill = 0
  ) %>%
  write_tsv("spikein_metrics/normalized_scalefactor_mirbase_counts.tsv")

# -------------------------------------------------------------
# Per-sample QC plots
# -------------------------------------------------------------
# For each sample:
# 1. Plot spike-in calibration curve (log–log) with R²
# 2. Plot miRNA reads with detection-range shading
for (sample in unique(mirna_long$sample_id)) {
  # Spike-ins for this sample
  spikes_this <- spikes_detected %>% filter(samples == sample)

  # Model for this sample
  fm <- per_sample_models %>%
    filter(sample_id == sample) %>%
    pull(model) %>%
    .[[1]]
  
  limits <- spikeins_stats %>%
    filter(sample_id == sample)
  
  spikein_lower <- limits$spikein_lower_limit[[1]]
  spikein_upper <- limits$spikein_upper_limit[[1]]
  rsq_value <- limits$rsq[[1]]
  percent_in_range <- limits$percent_in_range[[1]]
  
  spikes_this <- spikes_this %>%
    mutate(concentration = predict(fm, newdata = data.frame(counts = counts)))
  
  # ------------------- PLOT 1: Spike-in calibration -------------------
  p1 <- ggplot(spikes_this, aes(counts, concentration)) +
    geom_point(size = 3, color = "darkblue") +
    geom_line(aes(y = concentration), color = "blue", linetype = "dashed", linewidth = 0.5) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks(size = 0.2) +
    labs(
      title = "Spike-in Calibration",
      subtitle = paste0("R² = ", round(rsq_value, 4)),
      x = "Reads (log10)",
      y = "Concentration (log10 molecules/µL)"
    ) +
    theme_bw()
  
  # ------------------- PLOT 2: miRNA read distribution ----------------
  mirna_this <- mirna_long %>%
    filter(sample_id == sample) %>%
    mutate(
      range = ifelse(in_range, "in range",
                     ifelse(reads < spikein_lower, "too low", "too high"))
    )
  
  p2 <- ggplot(mirna_this, aes(x = sample_id, y = reads, color = range)) +
    # Shaded area for spike-in detection range
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = spikein_lower, ymax = spikein_upper,
             alpha = 0.15, fill = "gray80"
    )  +
    
    geom_point(alpha = 0.6, size = 2.5, position = position_jitter(width = 0.15))  +
    scale_y_log10() +
    scale_color_manual(
      values = c("in range" = "green", "too low" = "red", "too high" = "orange"),
      guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))
    ) +
    labs(
      title = "miRNA Read Distribution",
      subtitle = paste0(percent_in_range, "% miRNAs in range"),
      x = "Sample",
      y = "Reads (log10)",
      color = "Range"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.2)
    )
  # Save side-by-side plots
  g <- arrangeGrob(p1, p2, ncol = 2)
  ggsave(
    filename = paste0("spikein_metrics/qc_plots/", sample, "_spikein_qc.pdf"),
    plot = g,
    width = 13,
    height = 5.5
  )
  
}




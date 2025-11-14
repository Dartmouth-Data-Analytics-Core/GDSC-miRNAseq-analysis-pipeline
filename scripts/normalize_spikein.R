library(tidyverse)
library(ggplot2)
library(gridExtra)

if (!dir.exists("spikein_metrics/qc_plots")) {
  dir.create("spikein_metrics/qc_plots")
}


spikes_obs <- read_tsv("spikein_counts/spikein.readcounts.tsv")
mirna_counts <- read_tsv("mirbase_counts/mirbase.readcounts.tsv")


spikes_info <- tibble(
  spikein_ID = c("miND-01", "miND-02", "miND-03", "miND-04", "miND-05", "miND-06", "miND-07"),
  amol_concentration = c(20, 5, 1.25, 0.3125, 0.075, 0.01, 0.005))


amol <- 602214
spikeInsVolume <- 1
finalVolume <- 5


spikes_info <- spikes_info %>%
  mutate(amol_concentration = amol_concentration * spikeInsVolume / finalVolume)



spikes_with_conc <- spikes_obs %>%
  pivot_longer(cols = -c(spikein_ID), names_to = "samples", values_to = "counts") %>%
  left_join(spikes_info, by = "spikein_ID")

spikes_detected <- spikes_with_conc %>%
  filter(counts > 0)


spikeins_stats <- spikes_with_conc %>%
  group_by(samples) %>%
  do({
    df <- .
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
      df2 <- df %>% filter(counts > 0)
      fit <- lm(amol_concentration ~ 0 + counts, data = df2)
      pred_int <- suppressWarnings(predict(fit, interval = "prediction")) %>% as_tibble()
      df_with_pred <- df2 %>% bind_cols(pred_int)
      
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


per_sample_models <- spikes_detected %>%
  group_split(samples) %>%
  map_df(~ {
    df <- .x
    fm <- lm(amol_concentration ~ 0 + counts, data = df)
    tibble(
      sample_id = unique(df$samples),
      model = list(fm)
    )
  })


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


mirnas_in_range_count <- mirna_long %>%
  group_by(sample_id) %>%
  summarise(
    mirnas_in_range = sum(in_range, na.rm = TRUE),
    percent_in_range = round((sum(in_range, na.rm = TRUE) / n()) * 100, 2),
    .groups = "drop"
  )


spikeins_stats <- spikeins_stats %>%
  select(-mirnas_in_range) %>%
  left_join(mirnas_in_range_count, by = "sample_id") %>%
  select(
    sample_id, spikeins_detected, spikein_lower_limit, spikein_upper_limit,
    mirnas_in_range, percent_in_range, slope, rsq, qc
  )


LOD_by_sample <- spikes_detected %>%
  group_by(samples) %>%
  summarise(LOD = min(amol_concentration), .groups = "drop") %>%
  rename(sample_id = samples)


spikeins_stats %>%
  left_join(LOD_by_sample, by = "sample_id") %>%
  rename(
    "Scale Factor/ Slope" = slope,
    "Limit of Detection amol/uL (microliter)" = LOD,
    "Percent in range %" = percent_in_range,
    "Number of mirnas in Range" = mirnas_in_range,
    "Number of spike-ins detected" = spikeins_detected,
    "spike-in lower limit detection" = spikein_lower_limit,
    "spike-in upper limit detection" = spikein_upper_limit
  ) %>%
  select(-c(qc)) %>%
  write_tsv("spikein_metrics/spikein_detection_metrics.tsv")


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
    amol_concentration = round(fit, 2),
    interval = round(fit - lwr, 2)
  ) %>%
  ungroup() %>%
  select(mirbase_ID, Length, sample_id, reads, amol_concentration)


mirna_normalized %>%
  pivot_wider(
    names_from = sample_id,
    values_from = amol_concentration,
    id_cols = c(mirbase_ID, Length),
    values_fill = 0
  ) %>%
  write_tsv("spikein_metrics/normalized_scalefactor_mirbase_counts.tsv")


for (sample in unique(mirna_long$sample_id)) {
  
  spikes_this <- spikes_detected %>% filter(samples == sample)
  
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
  
  # PLOT 1
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
      y = "Concentration (log10 amol/µL)"
    ) +
    theme_bw()
  
  # PLOT 2
  mirna_this <- mirna_long %>%
    filter(sample_id == sample) %>%
    mutate(
      range = ifelse(in_range, "in range",
                     ifelse(reads < spikein_lower, "too low", "too high"))
    )
  
  p2 <- ggplot(mirna_this, aes(x = sample_id, y = reads, color = range)) +
    # Adicionar retângulo cinza para a região "in range"
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = spikein_lower, ymax = spikein_upper,
             alpha = 0.15, fill = "gray80"
    )  +
    # Pontos dos miRNAs
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
  
  g <- arrangeGrob(p1, p2, ncol = 2)
  ggsave(
    filename = paste0("spikein_metrics/qc_plots/", sample, "_spikein_qc.pdf"),
    plot = g,
    width = 13,
    height = 5.5
  )
  
}


# GDSC-miRNASeq-Analysis Pipeline v2 Outputs

#### Trimming Outputs (Default)

| Output Folder/File | Rule | Details |
|--------------------|------|---------|
| `trimming/{sample}.R1.trim.fastq.gz` | `trimming` | Trimmed fastq file with adapter sequences removed and minimum quality of Phred 30 |
| `triming/logs/{sample}.cutadapt.log` | `trimming` | CutAdapt trimming log |

#### Seqcluster Outputs (Default) 

| Output Folder/File | Rule | Details | 
|--------------------|------|---------|
| `collapsed/{sample}.seqcluster.fastq.gz` | `seqcluster` | Trimmed fastq reads collapsed to miRNA loci |
| `collapsed/logs/{sample}.seqcluster.log` | `seqcluster` | Seqcluster log file |

#### isomiR Mapping and Quantification Outputs (Default)

| Output Folder/File | Rule | Details | 
|--------------------|------|---------|
| `collapsed/{sample}.seqcluster.hairpin.aln.srt.bam` | `collapsed_hairpin_aln` | Seqcluster collapsed reads that have been aligned to miRbase hairpins with `Bowtie1` |
| `collapsed/{sample}.seqcluster.hairpin.aln.srt.bam.bai` | `collapsed_hairpin_aln` | Index file for seqcluster collapsed reads that have been aligned to miRbase hairpins |
| `collapsed/{sample}.unmapped.fastq` | `collapsed_hairpin_aln` | Seqcluster collapsed reads that failed to align to miRbase hairpin sequences |
| `alignment_logs/seqcluster/{sample}.bowtie1.hairpin.aln.log` | `collapsed_hairpin_aln` | Bowtie1 log for seqcluster collapsed reads to miRbase hairpin sequences |
| `collapsed/{sample}.seqcluster.hairpin.aln.srt.bam.idxstats` | `hairpin_stats` | Samtools idxstats for seqcluster collapsed reads hairpin alignment |
| `collapsed/{sample}.seqcluster.hairpin.aln.srt.bam.flagstat` | `hairpin_stats` | Samtools flagstats for seqcluster collapsed reads hairpin alignment |
| `mirtop/{sample}.hairpin.gff` | `miRtop` | miRTop generated isomiR annotation file in miR-gff3 format |
| `mirtop/{sample}.hairpin.tsv` | `miRtop` | miRTop generated isomiR annotation file in human-readable TSV format |
| `mirtop/logs/{sample}.mirtop.gff.log` | `miRtop` | miRTop gff command stderr/stdout log |
| `mirtop/logs/{sample}.mirtop.counts.log` | `miRtop` | miRTop counts command stderr/stdout log |
| `mirtop/mirtop_stats.log` | `mirtop_stats` | miRtop statistics |
| `mirtop/temp/{sample}.hairpin_long.csv` | `pivot_isomirs_longer` | IsomiR TSV file pivotted in long format to prepare for generating cohesive counts table |
| `mirtop/temp/master_counts_long.csv` | `collate_isomir_table` | Collated isomiR table in long format
| `miRNA_Quant/raw_merged_canonical_and_all_isomirs.csv` | `collate_isomiR_table` | miRNA counts (canonical + all isomiR classes grouped by miRNA) |
| `miRNA_Quant/raw_canonical_counts.tsv` | `collate_isomir_table` | miRNA counts (only canonical miRNAs) |
| `miRNA_Quant/raw_all_isomir_counts.tsv` | `collate_isomir_table` | miRNA isomiR counts, not grouped by miRNA (does not include canonical miRNAs) |
| `miRNA_Quant/<class>.isomir_counts.csv` | `collate_isomir_table` | miRNA isomiR counts for specific isomiR classes (not grouped by miRNA, does not include canonical) |

#### Genome Mapping and Quantification Outputs (Default)

| Output Folder/File | Rule | Details | 
|--------------------|------|---------|
| `mirbase_alignment/{sample}.mature.srt.bam` | `mirbase_padded_aln` | Trimmed reads to miRbase mature reference with 7bp padding on 3 and 5' ends with Bowtie2 for the purpose of getting unaligned reads (mapped reads are filtered for length, mismatches and quality score) |
| `mirbase_alignment/{sample}.mature.srt.bam.bai` | `mirbase_padded_aln` | Index file for reads that aligned to miRbase mature padded reference |
| `mirbase_alignment/{sample}.mature.unalign.fastq` | `mirbase_padded_aln` | Reads that did not align to miRbase mature padded reference that will be used as input for genome mapping |
| `alignment_logs/mirbase_mature_padded/{sample}.bowtie2.mature.log` | `mirbase_padded_aln` | Bowtie2 alignment log for miRbase mature padded alignment |
| `mirbase_alignment/{sample}.mature.srt.bam.idxstats` | `mature_mirbase_stats` | Samtools idxstats for filtered miRbase mature padded alignment |
| `mirbase_alignment/{sample}.mature.srt.bam.flagstat` | `mature_mirbase_stats` | Samtools flagstats for filtered miRbase mature padded alignemnt |
| `genome_alignment/{sample}.genome.srt.filt.bam` | `genome_alignment` | Reads that did not align to miRbase mature padded that aligned to the full genome with Bowtie2 |
| `alignment_logs/genome_alignment/{sample}.bowtie2.genome.log` | `genome_alignment` | Bowtie2 alignment log for genome alignment |
| `genome_alignment/{sample}.genome.srt.filt.bam.idxstats` | `genome_stats` | Samtools idxstats for genome alignment |
| `mirbase_alignment/{sample}.genome.srt.filt.bam.flagstat` | `genome_stats` | Samtools flagstats for genome alignment |
| `genome_counts/featurecounts.tsv` | `genome_featureCounts` | Raw unannotated genome + smRNA counts output (should be considered a temp file) |
| `genome_counts/featurecounts.readcounts.tsv` | `genome_featureCounts` | Raw, formatted genome + smRNA counts output (should be considered a temp file) |
| `genome_counts/featurecounts.readcounts.biotype.tsv` | `genome_featureCounts` | Raw, formatted and annotated genome + smRNA counts output with smRNA biotype information added (primary output) |
| `genome_counts/featurecounts.readcounts.biotype_rpkm.tsv` | `genome_featureCounts` | RPKM-normalized annotated genome + smRNA counts |
| `genome_counts/featurecounts.readcounts.biotype_tpm.tsv` | `genome_featureCounts` | TPM-normalized annotated genome + smRNA counts |
| `genome_counts/featurecounts.tsv.summary` | `genome_featureCounts` | FeatureCounts log |
| `metrics/mirna_genome_alignment_metrics.tsv` | `alignment_metrics_counts` | miRNA and genome alignment metrics in TSV format |
| `metrics/mirna_genome_alignment_metrics.xlsx` | `alignment_metrics_counts` | miRNA and genome alignment metrics in xlsx format |
| `metrics/smRNA_metrics.tsv` | `alignment_metrics_counts` | smRNA (contamination) metrics in TSV format |
| `metrics/smRNA_mqc.tsv` | `alignment_metrics_counts` | smRNA (contamination) metrics formatted for use in MultiQC report |

#### Principal Component Analysis (Default)

| Output Folder/File | Rule | Details | 
|--------------------|------|---------|
| `plots/PCA_top_PC1_vs_PC2.png` | `pca_plots` | PCA calculated on the top N-most highly variable miRNAs (dynamically calculated by variance plateau) |
| `plots/PCA_top_PCA_variance_bar.png` | `pca_plots` | Percent variance for each principal component (calculated on the top N-most highly variable miRNAs) |
| `plots/miRNA_Variance_Plot.png` | `pca_plots` | Variance profile of miRNAs + threshold where variance plateaus |
| `plots/PCA_all_PC1_vs_PC2.png` | `pca_plots` | PCA calculated on all profiled miRNAs |
| `plots/PCA_all_PCA_variance_bar.png` | `pca_plots` | Percent variance for each principal component (calculated on all miRNAs) |
| `plots/pca_hvg_log.txt` | `pca_plots` | Number of total miRNAs, number of top highly variable miRNAs, threshold information |
| `plots/Top_miRNAs_Heatmap.png` | `pca_plots` | Top N-most highly variable miRNAs visualized as a normalized, Z-scaled heatmap |

#### Spike-in Outputs (Optional)

If Spike-ins are added during library prep, the optional spike-in rules can be included. 7 spike-ins are included. Any sample where 2 or more spike-ins are not-detected is considered a QC failure. This will be represented in the `spikein_detection_metrics.tsv` file as well as the associated QC plot for that sample (which will show up as grey and "NA".
If any samples have 0 counts across all spike-ins, this sample is removed and included in the `removed_samples_no_spikeins.txt` file. 
The merged canonical and isomiR counts are used as the input data for spike-in normalization. If using spike-ins, the associated output, (`normalized_scalefactor_canon_and_isomir_counts.tsv`) should be used for downstream analysis.

| Output Folder/File | Rule | Details | 
|--------------------|------|---------|
| `spikein_alignment/{sample}.stats` | `spikein_bbduk` | Spike-in alignment statistics |
| `spikein_alignment/{sample}.unmapped.fastq.gz` | `spikein_bbduk` | Reads that did not align to spike-ins |
| `spikein_counts/spikein.readcounts.tsv` | `spikein_counts` | Counts per spike-in for each sample |
| `spikein_metrics/spikein_detection_metrics.tsv` | `normalize_data_spikein` | Detection metrics for spike-ins |
| `spikein_metrics/normalized_scalefactor_canon_and_isomir_counts.tsv` | `normalize_data_spikein` | Spike-in normalized data for the merged canonical and isomiR counts |
| `spikein_metrics/qc_plots/{sample}_spikein_qc.pdf` | `normalize_data_spikein` | QC plots for spike-ins per sample |
| `spikein_metrics/removed_samples_no_spikeins.txt` | `normalize_data_spikein` | Samples removed due to no spike-ins |

#!/usr/bin/env Rscript
# ============================================================
# 01_data_overview.R
# Kcnq1ot1 & Rmst in Chronic Liver Disease - Dry Analysis
# Data loading, QC, and initial exploration
# ============================================================

library(ggplot2)
library(pheatmap)
library(DESeq2)
library(edgeR)

data_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/data"
results_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/results"
fig_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/figures"

# ============================================================
# 1. Load GSE246221 (STZ+HFD MASLD→HCC time course)
# ============================================================
cat("Loading GSE246221...\n")
gse246 <- read.delim(
  gzfile(file.path(data_dir, "GSE246221_rawcounts_allsamples_mouse.txt.gz")),
  row.names = 1, check.names = FALSE
)
cat(sprintf("  Dimensions: %d genes x %d samples\n", nrow(gse246), ncol(gse246)))

# Parse sample metadata from column names
parse_gse246_meta <- function(sample_names) {
  meta <- data.frame(
    sample = sample_names,
    stringsAsFactors = FALSE
  )
  # Extract batch, condition, and replicate
  meta$batch <- gsub("^(Batch\\d+)_.*", "\\1", sample_names)
  meta$raw_condition <- gsub("^Batch\\d+_", "", sample_names)

  # Classify by disease stage
  meta$stage <- NA
  meta$timepoint_weeks <- NA
  meta$is_tumor <- grepl("_T\\d+$", sample_names)

  # Batch1: main time course
  meta$stage[grepl("Control07w", sample_names)] <- "Healthy"
  meta$timepoint_weeks[grepl("Control07w", sample_names)] <- 7

  meta$stage[grepl("STZSCD08w", sample_names)] <- "Acute_STZ"
  meta$timepoint_weeks[grepl("STZSCD08w", sample_names)] <- 8

  meta$stage[grepl("STZHFD14w", sample_names)] <- "Early_MASLD"
  meta$timepoint_weeks[grepl("STZHFD14w", sample_names)] <- 14

  meta$stage[grepl("STZHFD20w", sample_names) & meta$batch == "Batch1"] <- "MASH"
  meta$timepoint_weeks[grepl("STZHFD20w", sample_names) & meta$batch == "Batch1"] <- 20

  meta$stage[grepl("STZHFD32w", sample_names) & meta$batch == "Batch1"] <- "Fibrosis"
  meta$timepoint_weeks[grepl("STZHFD32w", sample_names) & meta$batch == "Batch1"] <- 32

  meta$stage[grepl("STZHFD44w", sample_names) & meta$batch == "Batch1" & !meta$is_tumor] <- "Early_HCC_NonTumor"
  meta$stage[grepl("STZHFD44w", sample_names) & meta$batch == "Batch1" & meta$is_tumor] <- "Early_HCC_Tumor"
  meta$timepoint_weeks[grepl("STZHFD44w", sample_names) & meta$batch == "Batch1"] <- 44

  meta$stage[grepl("STZHFD50w", sample_names) & meta$batch == "Batch1" & !meta$is_tumor] <- "Mid_HCC_NonTumor"
  meta$stage[grepl("STZHFD50w", sample_names) & meta$batch == "Batch1" & meta$is_tumor] <- "Mid_HCC_Tumor"
  meta$timepoint_weeks[grepl("STZHFD50w", sample_names) & meta$batch == "Batch1"] <- 50

  meta$stage[grepl("STZHFD56w", sample_names) & meta$batch == "Batch1" & !meta$is_tumor] <- "Late_HCC_NonTumor"
  meta$stage[grepl("STZHFD56w", sample_names) & meta$batch == "Batch1" & meta$is_tumor] <- "Late_HCC_Tumor"
  meta$timepoint_weeks[grepl("STZHFD56w", sample_names) & meta$batch == "Batch1"] <- 56

  # Batch2: diet comparison at 20w
  meta$stage[grepl("HFDonly20w", sample_names)] <- "HFD_only"
  meta$timepoint_weeks[grepl("HFDonly20w", sample_names)] <- 20

  meta$stage[grepl("SCDonly20w", sample_names)] <- "SCD_only"
  meta$timepoint_weeks[grepl("SCDonly20w", sample_names)] <- 20

  meta$stage[grepl("STZSCD20w", sample_names)] <- "STZ_SCD"
  meta$timepoint_weeks[grepl("STZSCD20w", sample_names)] <- 20

  # Batch4-5: diet change experiments
  meta$stage[grepl("Dietchange", sample_names)] <- "Diet_Reversal"
  meta$timepoint_weeks[grepl("Dietchange26w", sample_names)] <- 26
  meta$timepoint_weeks[grepl("Dietchange50w", sample_names)] <- 50

  meta$stage[grepl("STZHFD26w", sample_names)] <- "STZHFD_26w"
  meta$timepoint_weeks[grepl("STZHFD26w", sample_names)] <- 26

  meta$stage[grepl("STZHFD50w", sample_names) & meta$batch == "Batch5"] <- "STZHFD_50w"
  meta$timepoint_weeks[grepl("STZHFD50w", sample_names) & meta$batch == "Batch5"] <- 50

  # Batch6-8: Tirzepatide treatment
  meta$stage[grepl("Tirzepatide", sample_names)] <- "Tirzepatide"
  meta$stage[grepl("Vehicle", sample_names)] <- "Vehicle"
  meta$timepoint_weeks[grepl("STZHFD32w_", sample_names) & meta$batch == "Batch6"] <- 32
  meta$timepoint_weeks[grepl("STZHFD38w_", sample_names)] <- 38
  meta$timepoint_weeks[grepl("STZHFD52w_", sample_names)] <- 52

  # Create main analysis group (Batch1 disease progression only)
  meta$main_timecourse <- meta$batch == "Batch1"

  # Disease category for broader grouping
  meta$disease_category <- "Other"
  meta$disease_category[meta$stage == "Healthy"] <- "1_Healthy"
  meta$disease_category[meta$stage == "Acute_STZ"] <- "2_Acute"
  meta$disease_category[meta$stage == "Early_MASLD"] <- "3_MASLD"
  meta$disease_category[meta$stage == "MASH"] <- "4_MASH"
  meta$disease_category[meta$stage == "Fibrosis"] <- "5_Fibrosis"
  meta$disease_category[meta$stage %in% c("Early_HCC_NonTumor", "Mid_HCC_NonTumor", "Late_HCC_NonTumor")] <- "6_HCC_NonTumor"
  meta$disease_category[meta$stage %in% c("Early_HCC_Tumor", "Mid_HCC_Tumor", "Late_HCC_Tumor")] <- "7_HCC_Tumor"

  rownames(meta) <- meta$sample
  return(meta)
}

meta246 <- parse_gse246_meta(colnames(gse246))
cat("\nGSE246221 sample distribution:\n")
print(table(meta246$stage))
cat("\nBatch1 time course samples:\n")
print(table(meta246$stage[meta246$main_timecourse]))

# ============================================================
# 2. Load GSE207855 (CCl4 liver fibrosis)
# ============================================================
cat("\nLoading GSE207855...\n")
gse207 <- read.delim(
  gzfile(file.path(data_dir, "GSE207855_count.tsv.gz")),
  row.names = 1, check.names = FALSE
)
cat(sprintf("  Dimensions: %d genes x %d samples\n", nrow(gse207), ncol(gse207)))

# Parse metadata
meta207 <- data.frame(
  sample = colnames(gse207),
  condition = ifelse(grepl("^CCl4", colnames(gse207)), "CCl4", "Oil"),
  stringsAsFactors = FALSE
)
rownames(meta207) <- meta207$sample
cat("\nGSE207855 sample distribution:\n")
print(table(meta207$condition))

# ============================================================
# 3. Extract Kcnq1ot1 and Rmst expression
# ============================================================
cat("\n=== Target lncRNA Expression ===\n")

# GSE246221
kcnq1ot1_246 <- as.numeric(gse246["ENSMUSG00000101609_Kcnq1ot1", ])
rmst_246 <- as.numeric(gse246["ENSMUSG00000112117_Rmst", ])
names(kcnq1ot1_246) <- colnames(gse246)
names(rmst_246) <- colnames(gse246)

cat("\nGSE246221 - Kcnq1ot1 summary:\n")
cat(sprintf("  Range: %.1f - %.1f\n", min(kcnq1ot1_246), max(kcnq1ot1_246)))
cat(sprintf("  Mean: %.1f, Median: %.1f\n", mean(kcnq1ot1_246), median(kcnq1ot1_246)))

cat("\nGSE246221 - Rmst summary:\n")
cat(sprintf("  Range: %d - %d\n", min(rmst_246), max(rmst_246)))
cat(sprintf("  Mean: %.2f (essentially undetectable)\n", mean(rmst_246)))

# GSE207855
kcnq1ot1_207 <- as.numeric(gse207["63830", ])
rmst_207 <- as.numeric(gse207["74148", ])
names(kcnq1ot1_207) <- colnames(gse207)
names(rmst_207) <- colnames(gse207)

cat("\nGSE207855 - Kcnq1ot1 (Entrez 63830):\n")
cat(sprintf("  CCl4 mean: %.1f, Oil mean: %.1f\n",
    mean(kcnq1ot1_207[meta207$condition == "CCl4"]),
    mean(kcnq1ot1_207[meta207$condition == "Oil"])))

cat("\nGSE207855 - Rmst (Entrez 74148):\n")
cat(sprintf("  CCl4 mean: %.1f, Oil mean: %.1f\n",
    mean(rmst_207[meta207$condition == "CCl4"]),
    mean(rmst_207[meta207$condition == "Oil"])))

# ============================================================
# 4. Kcnq1ot1 expression across disease stages (GSE246221)
# ============================================================

# Focus on Batch1 main time course
tc_idx <- meta246$main_timecourse
tc_meta <- meta246[tc_idx, ]
tc_counts <- gse246[, tc_idx]

# Normalize with edgeR (TMM)
cat("\nNormalizing GSE246221 (Batch1, TMM)...\n")
y_tc <- DGEList(counts = round(tc_counts))  # round fractional counts
y_tc <- calcNormFactors(y_tc)
tc_cpm <- cpm(y_tc, log = TRUE)

# Extract Kcnq1ot1 CPM values
kcnq1ot1_cpm <- tc_cpm["ENSMUSG00000101609_Kcnq1ot1", ]

# Create plot data
plot_data <- data.frame(
  sample = names(kcnq1ot1_cpm),
  logCPM = kcnq1ot1_cpm,
  stage = tc_meta[names(kcnq1ot1_cpm), "stage"],
  disease_category = tc_meta[names(kcnq1ot1_cpm), "disease_category"],
  timepoint = tc_meta[names(kcnq1ot1_cpm), "timepoint_weeks"],
  is_tumor = tc_meta[names(kcnq1ot1_cpm), "is_tumor"]
)
plot_data <- plot_data[plot_data$disease_category != "Other", ]

# Order factor levels
plot_data$disease_category <- factor(plot_data$disease_category,
  levels = c("1_Healthy", "2_Acute", "3_MASLD", "4_MASH",
             "5_Fibrosis", "6_HCC_NonTumor", "7_HCC_Tumor"))

cat("\nKcnq1ot1 log2CPM by disease stage:\n")
tapply(plot_data$logCPM, plot_data$disease_category, function(x) {
  sprintf("mean=%.2f, sd=%.2f, n=%d", mean(x), sd(x), length(x))
})

# ============================================================
# 5. Save results
# ============================================================
# Save metadata
write.csv(meta246, file.path(results_dir, "GSE246221_metadata.csv"), row.names = FALSE)
write.csv(meta207, file.path(results_dir, "GSE207855_metadata.csv"), row.names = FALSE)

# Save Kcnq1ot1 plot data
write.csv(plot_data, file.path(results_dir, "Kcnq1ot1_expression_by_stage.csv"), row.names = FALSE)

# ============================================================
# 6. Visualization: Kcnq1ot1 expression across disease stages
# ============================================================

# Color palette
stage_colors <- c(
  "1_Healthy" = "#2166AC",
  "2_Acute" = "#67A9CF",
  "3_MASLD" = "#D1E5F0",
  "4_MASH" = "#FDDBC7",
  "5_Fibrosis" = "#EF8A62",
  "6_HCC_NonTumor" = "#B2182B",
  "7_HCC_Tumor" = "#67001F"
)

stage_labels <- c(
  "1_Healthy" = "Healthy\n(7w)",
  "2_Acute" = "Acute STZ\n(8w)",
  "3_MASLD" = "MASLD\n(14w)",
  "4_MASH" = "MASH\n(20w)",
  "5_Fibrosis" = "Fibrosis\n(32w)",
  "6_HCC_NonTumor" = "HCC\nNon-Tumor",
  "7_HCC_Tumor" = "HCC\nTumor"
)

p1 <- ggplot(plot_data, aes(x = disease_category, y = logCPM, fill = disease_category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = stage_colors) +
  scale_x_discrete(labels = stage_labels) +
  labs(
    title = "Kcnq1ot1 Expression Across Disease Progression (GSE246221)",
    subtitle = "STZ+HFD Mouse Model: Healthy → MASLD → MASH → Fibrosis → HCC",
    x = "Disease Stage",
    y = "log2(CPM)",
    caption = "Data: GSE246221 (Batch1 time course, TMM-normalized)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 9)
  )

ggsave(file.path(fig_dir, "Fig1_Kcnq1ot1_disease_progression.pdf"), p1,
       width = 10, height = 6)
ggsave(file.path(fig_dir, "Fig1_Kcnq1ot1_disease_progression.png"), p1,
       width = 10, height = 6, dpi = 300)

cat("\nFigure saved: Fig1_Kcnq1ot1_disease_progression\n")

# ============================================================
# 7. GSE207855: CCl4 vs Oil comparison
# ============================================================
cat("\n=== GSE207855: DEG analysis (CCl4 vs Oil) ===\n")

# DESeq2 analysis
dds207 <- DESeqDataSetFromMatrix(
  countData = gse207,
  colData = meta207,
  design = ~ condition
)
dds207$condition <- relevel(dds207$condition, ref = "Oil")
dds207 <- DESeq(dds207)
res207 <- results(dds207, contrast = c("condition", "CCl4", "Oil"))

# Kcnq1ot1 and Rmst results
cat("\nKcnq1ot1 (63830) - DESeq2 results:\n")
print(res207["63830", ])
cat("\nRmst (74148) - DESeq2 results:\n")
print(res207["74148", ])

# Normalize for visualization
y207 <- DGEList(counts = gse207)
y207 <- calcNormFactors(y207)
cpm207 <- cpm(y207, log = TRUE)

# Plot CCl4 vs Oil for both lncRNAs
plot_207 <- data.frame(
  sample = rep(colnames(gse207), 2),
  gene = rep(c("Kcnq1ot1", "Rmst"), each = ncol(gse207)),
  logCPM = c(cpm207["63830", ], cpm207["74148", ]),
  condition = rep(meta207$condition, 2)
)

p2 <- ggplot(plot_207, aes(x = condition, y = logCPM, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = c("Oil" = "#2166AC", "CCl4" = "#B2182B")) +
  labs(
    title = "Kcnq1ot1 & Rmst Expression in CCl4-Induced Liver Fibrosis (GSE207855)",
    x = "Treatment",
    y = "log2(CPM)",
    caption = "Data: GSE207855 (6-week CCl4 IP injection, TMM-normalized)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "italic", size = 12)
  )

ggsave(file.path(fig_dir, "Fig2_CCl4_Kcnq1ot1_Rmst.pdf"), p2,
       width = 8, height = 5)
ggsave(file.path(fig_dir, "Fig2_CCl4_Kcnq1ot1_Rmst.png"), p2,
       width = 8, height = 5, dpi = 300)

cat("\nFigure saved: Fig2_CCl4_Kcnq1ot1_Rmst\n")

# ============================================================
# 8. Summary statistics
# ============================================================
cat("\n============================================================\n")
cat("SUMMARY OF INITIAL DATA EXPLORATION\n")
cat("============================================================\n")
cat("\nDatasets:\n")
cat("  GSE246221: 55,536 genes x 106 samples (STZ+HFD model)\n")
cat("  GSE207855: 36,528 genes x 12 samples (CCl4 fibrosis)\n")
cat("\nTarget lncRNA availability:\n")
cat("  Kcnq1ot1: Detected in BOTH datasets\n")
cat("  Rmst: Detected in GSE207855 only (undetectable in GSE246221)\n")

# Save DESeq2 results for GSE207855
res207_df <- as.data.frame(res207)
res207_df$gene_id <- rownames(res207_df)
write.csv(res207_df, file.path(results_dir, "GSE207855_DESeq2_CCl4_vs_Oil.csv"),
          row.names = FALSE)

cat("\nAll results saved to:", results_dir, "\n")
cat("All figures saved to:", fig_dir, "\n")
cat("\nDone!\n")

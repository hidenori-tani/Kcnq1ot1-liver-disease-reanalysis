#!/usr/bin/env Rscript
# ============================================================
# 04_human_NAFLD_severity.R
# Step 2: Human NAFLD severity analysis for KCNQ1OT1
# ============================================================

library(ggplot2)
library(edgeR)
library(DESeq2)

data_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/data"
results_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/results"
fig_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/figures"

# ============================================================
# 1. Load human data from GSE246221
# ============================================================
cat("Loading human NAFLD data...\n")
human <- read.delim(
  gzfile(file.path(data_dir, "GSE246221_rawcounts_allsamples_human.txt.gz")),
  row.names = 1, check.names = FALSE
)
cat(sprintf("  Dimensions: %d genes x %d samples\n", nrow(human), ncol(human)))
cat(sprintf("  Samples: %s\n", paste(colnames(human), collapse = ", ")))

# ============================================================
# 2. Fetch sample metadata from GEO series matrix
# ============================================================
# The human samples: NAFLD_01-28 + Healthy_01-04
# From the paper (Nature Comms 2024): patients had varying NAS and fibrosis stages
# We need to parse the sample numbers to infer severity grouping

# Based on the GSE246221 publication:
# Samples are ordered by disease severity
# NAFLD patients: NAS scores and fibrosis stages vary
# Let's extract what we can from the data itself

meta_human <- data.frame(
  sample = colnames(human),
  group = ifelse(grepl("Healthy", colnames(human)), "Healthy", "NAFLD"),
  sample_num = as.numeric(gsub("Human_(NAFLD|Healthy)_", "", colnames(human))),
  stringsAsFactors = FALSE
)
rownames(meta_human) <- meta_human$sample

# ============================================================
# 3. Normalize and analyze KCNQ1OT1 + RMST
# ============================================================
cat("\nNormalizing human data...\n")
y_h <- DGEList(counts = round(human))
y_h <- calcNormFactors(y_h)
cpm_h <- cpm(y_h, log = TRUE)

kcnq_h <- cpm_h["ENSG00000269821_KCNQ1OT1", ]
rmst_h <- cpm_h["ENSG00000255794_RMST", ]

# ============================================================
# 4. DESeq2: NAFLD vs Healthy
# ============================================================
cat("\nDESeq2: NAFLD vs Healthy...\n")

dds_h <- DESeqDataSetFromMatrix(
  countData = round(human),
  colData = meta_human,
  design = ~ group
)
# Filter low counts
keep <- rowSums(counts(dds_h) >= 10) >= 3
dds_h <- dds_h[keep, ]
dds_h$group <- relevel(factor(dds_h$group), ref = "Healthy")
dds_h <- DESeq(dds_h)
res_h <- results(dds_h, contrast = c("group", "NAFLD", "Healthy"))

cat("\nKCNQ1OT1 - DESeq2 result:\n")
print(res_h["ENSG00000269821_KCNQ1OT1", ])

cat("\nRMST - DESeq2 result:\n")
if ("ENSG00000255794_RMST" %in% rownames(res_h)) {
  print(res_h["ENSG00000255794_RMST", ])
} else {
  cat("  Filtered out (low expression)\n")
}

# ============================================================
# 5. Full DEG analysis: Human NAFLD
# ============================================================
res_h_df <- as.data.frame(res_h)
res_h_df$gene_id <- rownames(res_h_df)
res_h_df$gene_symbol <- gsub("^ENSG\\d+_", "", res_h_df$gene_id)

sig_h <- res_h_df[!is.na(res_h_df$padj) & res_h_df$padj < 0.05, ]
cat(sprintf("\nHuman NAFLD vs Healthy: %d significant DEGs (padj<0.05)\n", nrow(sig_h)))
cat(sprintf("  Upregulated: %d\n", nrow(sig_h[sig_h$log2FoldChange > 0, ])))
cat(sprintf("  Downregulated: %d\n", nrow(sig_h[sig_h$log2FoldChange < 0, ])))

# Check known liver disease lncRNAs
known_human_lnc <- c("KCNQ1OT1", "RMST", "NEAT1", "GAS5", "HOTAIR", "MALAT1", "H19", "MEG3")
cat("\nKnown lncRNAs in Human NAFLD vs Healthy:\n")
cat(sprintf("%-12s %10s %10s %12s %12s\n", "Gene", "log2FC", "baseMean", "pvalue", "padj"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (lnc in known_human_lnc) {
  lnc_row <- grep(paste0("_", lnc, "$"), res_h_df$gene_id, value = TRUE)
  if (length(lnc_row) > 0) {
    r <- res_h_df[res_h_df$gene_id == lnc_row[1], ]
    padj_val <- ifelse(is.na(r$padj), "NA", sprintf("%.2e", r$padj))
    cat(sprintf("%-12s %10.4f %10.1f %12.2e %12s\n",
        lnc, r$log2FoldChange, r$baseMean, r$pvalue, padj_val))
  }
}

write.csv(res_h_df, file.path(results_dir, "Human_NAFLD_vs_Healthy_DESeq2.csv"),
          row.names = FALSE)

# ============================================================
# 6. Unsupervised clustering to infer severity
# ============================================================
cat("\n=== Unsupervised severity analysis ===\n")

# Cluster NAFLD samples by overall expression pattern
nafld_idx <- meta_human$group == "NAFLD"
nafld_cpm <- cpm_h[, nafld_idx]

# Use top variable genes for clustering
rv <- apply(nafld_cpm, 1, var)
top_var <- names(sort(rv, decreasing = TRUE))[1:2000]
nafld_top <- nafld_cpm[top_var, ]

# Hierarchical clustering
d <- dist(t(nafld_top))
hc <- hclust(d, method = "ward.D2")
clusters <- cutree(hc, k = 3)

# KCNQ1OT1 by cluster
cat("\nKCNQ1OT1 expression by unsupervised cluster:\n")
for (cl in sort(unique(clusters))) {
  samples_in <- names(clusters[clusters == cl])
  vals <- kcnq_h[samples_in]
  cat(sprintf("  Cluster %d (n=%d): mean logCPM=%.3f, sd=%.3f\n",
      cl, length(vals), mean(vals), sd(vals)))
}

# ANOVA
kcnq_nafld <- kcnq_h[nafld_idx]
anova_res <- summary(aov(kcnq_nafld ~ factor(clusters)))
cat(sprintf("  ANOVA p-value: %.4f\n", anova_res[[1]]$`Pr(>F)`[1]))

# ============================================================
# 7. Correlation with NAS-related gene signatures
# ============================================================
cat("\n=== KCNQ1OT1 correlation with disease severity markers ===\n")

# Known severity markers in NAFLD
severity_markers <- c(
  # Inflammation
  "TNF", "IL1B", "IL6", "CCL2",
  # Fibrosis
  "COL1A1", "COL1A2", "COL3A1", "ACTA2", "TGFB1",
  # Steatosis
  "FASN", "SREBF1", "PPARG", "CD36",
  # Hepatocyte damage
  "CYP2E1", "ALB",
  # Mitochondrial
  "NDUFA9", "SDHA", "SDHB", "ATP5F1A", "UQCRC2"
)

cat(sprintf("%-12s %10s %10s\n", "Gene", "Correlation", "p-value"))
cat(paste(rep("-", 36), collapse = ""), "\n")

severity_cors <- data.frame()
for (marker in severity_markers) {
  marker_row <- grep(paste0("_", marker, "$"), rownames(cpm_h), value = TRUE)
  if (length(marker_row) > 0) {
    ct <- cor.test(kcnq_h, cpm_h[marker_row[1], ], method = "pearson")
    cat(sprintf("%-12s %10.3f %10.4f %s\n",
        marker, ct$estimate, ct$p.value,
        ifelse(ct$p.value < 0.05, "*", "")))
    severity_cors <- rbind(severity_cors, data.frame(
      gene = marker, correlation = ct$estimate, pvalue = ct$p.value
    ))
  }
}

write.csv(severity_cors, file.path(results_dir, "KCNQ1OT1_severity_marker_correlations.csv"),
          row.names = FALSE)

# ============================================================
# 8. Figure: Human KCNQ1OT1 expression
# ============================================================

# Box plot: Healthy vs NAFLD
plot_h <- data.frame(
  sample = names(kcnq_h),
  logCPM = kcnq_h,
  group = meta_human[names(kcnq_h), "group"]
)

p_human <- ggplot(plot_h, aes(x = group, y = logCPM, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.7) +
  scale_fill_manual(values = c("Healthy" = "#2166AC", "NAFLD" = "#B2182B")) +
  annotate("text", x = 1.5, y = max(plot_h$logCPM) + 0.1,
           label = sprintf("log2FC=%.2f\np=%.4f (DESeq2)",
             as.numeric(res_h["ENSG00000269821_KCNQ1OT1", "log2FoldChange"]),
             as.numeric(res_h["ENSG00000269821_KCNQ1OT1", "padj"])),
           size = 3.5) +
  labs(
    title = "KCNQ1OT1 Expression in Human NAFLD",
    subtitle = "GSE246221 (4 Healthy donors + 28 NAFLD patients)",
    x = "",
    y = "KCNQ1OT1 log2(CPM)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(fig_dir, "Fig8_Human_KCNQ1OT1.pdf"), p_human, width = 5, height = 5)
ggsave(file.path(fig_dir, "Fig8_Human_KCNQ1OT1.png"), p_human, width = 5, height = 5, dpi = 300)
cat("\n  Fig8 saved.\n")

# ============================================================
# 9. Figure: Severity marker correlation heatmap
# ============================================================

if (nrow(severity_cors) > 0) {
  severity_cors$sig <- ifelse(severity_cors$pvalue < 0.05, "*", "")
  severity_cors$category <- NA
  severity_cors$category[severity_cors$gene %in% c("TNF", "IL1B", "IL6", "CCL2")] <- "Inflammation"
  severity_cors$category[severity_cors$gene %in% c("COL1A1", "COL1A2", "COL3A1", "ACTA2", "TGFB1")] <- "Fibrosis"
  severity_cors$category[severity_cors$gene %in% c("FASN", "SREBF1", "PPARG", "CD36")] <- "Steatosis"
  severity_cors$category[severity_cors$gene %in% c("CYP2E1", "ALB")] <- "Hepatocyte"
  severity_cors$category[severity_cors$gene %in% c("NDUFA9", "SDHA", "SDHB", "ATP5F1A", "UQCRC2")] <- "Mitochondria"

  severity_cors$gene <- factor(severity_cors$gene, levels = rev(severity_cors$gene))

  p_cor <- ggplot(severity_cors, aes(x = correlation, y = gene, fill = category)) +
    geom_col(alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_text(aes(label = sig, x = correlation + sign(correlation) * 0.02),
              size = 5, vjust = 0.5) +
    scale_fill_manual(values = c(
      "Inflammation" = "#E41A1C", "Fibrosis" = "#FF7F00",
      "Steatosis" = "#984EA3", "Hepatocyte" = "#4DAF4A",
      "Mitochondria" = "#377EB8"
    )) +
    labs(
      title = "KCNQ1OT1 Correlation with Disease Severity Markers",
      subtitle = "Human NAFLD + Healthy (n=32)",
      x = "Pearson Correlation with KCNQ1OT1",
      y = "",
      fill = "Category"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(fig_dir, "Fig9_KCNQ1OT1_severity_markers.pdf"), p_cor,
         width = 8, height = 6)
  ggsave(file.path(fig_dir, "Fig9_KCNQ1OT1_severity_markers.png"), p_cor,
         width = 8, height = 6, dpi = 300)
  cat("  Fig9 saved.\n")
}

cat("\nStep 2 (Human NAFLD severity analysis) complete!\n")

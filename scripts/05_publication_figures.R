#!/usr/bin/env Rscript
# ============================================================
# 05_publication_figures.R
# Step 3: Publication-quality figures
# ============================================================

library(ggplot2)
library(edgeR)
library(EnhancedVolcano)
library(DESeq2)
library(pheatmap)
# library(gridExtra)  # not needed

data_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/data"
results_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/results"
fig_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/figures"

# ============================================================
# Load all data
# ============================================================
cat("Loading data...\n")
gse246 <- read.delim(gzfile(file.path(data_dir, "GSE246221_rawcounts_allsamples_mouse.txt.gz")),
                     row.names = 1, check.names = FALSE)
human <- read.delim(gzfile(file.path(data_dir, "GSE246221_rawcounts_allsamples_human.txt.gz")),
                    row.names = 1, check.names = FALSE)

# Load DEG results
res_mash <- read.csv(file.path(results_dir, "DEG_MASH_vs_Healthy.csv"))
res_fib <- read.csv(file.path(results_dir, "DEG_Fibrosis_vs_Healthy.csv"))
res_hcc <- read.csv(file.path(results_dir, "DEG_HCC_Tumor_vs_Healthy.csv"))
res_human <- read.csv(file.path(results_dir, "Human_NAFLD_vs_Healthy_DESeq2.csv"))

# ============================================================
# Figure A: Volcano plot - Mouse MASH vs Healthy
# ============================================================
cat("Creating Volcano plots...\n")

# Prepare for EnhancedVolcano
prep_volcano <- function(df) {
  rownames(df) <- df$gene_id
  df$gene_label <- df$gene_symbol
  return(df)
}

res_mash_v <- prep_volcano(res_mash)
res_hcc_v <- prep_volcano(res_hcc)
res_human_v <- prep_volcano(res_human)

# Mark Kcnq1ot1
highlight_mouse <- "ENSMUSG00000101609_Kcnq1ot1"
highlight_human <- "ENSG00000269821_KCNQ1OT1"

# Known lncRNAs to label
mouse_lnc <- c("ENSMUSG00000101609_Kcnq1ot1", "ENSMUSG00000092274_Neat1",
               "ENSMUSG00000053332_Gas5")
human_lnc <- c("ENSG00000269821_KCNQ1OT1", "ENSG00000245532_NEAT1",
               "ENSG00000234741_GAS5")

# Volcano: MASH vs Healthy
p_v1 <- EnhancedVolcano(res_mash_v,
  lab = res_mash_v$gene_label,
  selectLab = gsub("^ENSMUSG\\d+_", "", mouse_lnc),
  x = "log2FoldChange", y = "padj",
  title = "Mouse MASH vs Healthy",
  subtitle = "GSE246221 (STZ+HFD, 20w)",
  pCutoff = 0.05, FCcutoff = 1,
  pointSize = 1.5, labSize = 3.5,
  col = c("grey80", "#2166AC", "#D6604D", "#B2182B"),
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  maxoverlapsConnectors = 20
) + theme(plot.title = element_text(size = 12))

ggsave(file.path(fig_dir, "FigA_Volcano_MASH.pdf"), p_v1, width = 8, height = 6)
ggsave(file.path(fig_dir, "FigA_Volcano_MASH.png"), p_v1, width = 8, height = 6, dpi = 300)
cat("  FigA saved.\n")

# Volcano: HCC Tumor vs Healthy
p_v2 <- EnhancedVolcano(res_hcc_v,
  lab = res_hcc_v$gene_label,
  selectLab = gsub("^ENSMUSG\\d+_", "", mouse_lnc),
  x = "log2FoldChange", y = "padj",
  title = "Mouse HCC Tumor vs Healthy",
  subtitle = "GSE246221 (STZ+HFD, 44-56w tumors)",
  pCutoff = 0.05, FCcutoff = 1,
  pointSize = 1.5, labSize = 3.5,
  col = c("grey80", "#2166AC", "#D6604D", "#B2182B"),
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  maxoverlapsConnectors = 20
) + theme(plot.title = element_text(size = 12))

ggsave(file.path(fig_dir, "FigB_Volcano_HCC.pdf"), p_v2, width = 8, height = 6)
ggsave(file.path(fig_dir, "FigB_Volcano_HCC.png"), p_v2, width = 8, height = 6, dpi = 300)
cat("  FigB saved.\n")

# Volcano: Human NAFLD vs Healthy
p_v3 <- EnhancedVolcano(res_human_v,
  lab = res_human_v$gene_symbol,
  selectLab = gsub("^ENSG\\d+_", "", human_lnc),
  x = "log2FoldChange", y = "padj",
  title = "Human NAFLD vs Healthy",
  subtitle = "GSE246221 (28 NAFLD + 4 Healthy)",
  pCutoff = 0.05, FCcutoff = 1,
  pointSize = 1.5, labSize = 3.5,
  col = c("grey80", "#2166AC", "#D6604D", "#B2182B"),
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  maxoverlapsConnectors = 20
) + theme(plot.title = element_text(size = 12))

ggsave(file.path(fig_dir, "FigC_Volcano_Human_NAFLD.pdf"), p_v3, width = 8, height = 6)
ggsave(file.path(fig_dir, "FigC_Volcano_Human_NAFLD.png"), p_v3, width = 8, height = 6, dpi = 300)
cat("  FigC saved.\n")

# ============================================================
# Figure D: Integrated summary - Kcnq1ot1/KCNQ1OT1 across all contexts
# ============================================================
cat("Creating integrated summary figure...\n")

# Collect all comparisons
summary_data <- data.frame(
  comparison = c(
    "MASLD vs\nHealthy", "MASH vs\nHealthy", "Fibrosis vs\nHealthy",
    "HCC NonTumor\nvs Healthy", "HCC Tumor\nvs Healthy",
    "Human NAFLD\nvs Healthy"
  ),
  log2FC = c(
    res_mash_v["ENSMUSG00000101609_Kcnq1ot1", "log2FoldChange"],
    res_mash_v["ENSMUSG00000101609_Kcnq1ot1", "log2FoldChange"],  # placeholder
    res_fib[res_fib$gene_id == "ENSMUSG00000101609_Kcnq1ot1", "log2FoldChange"],
    NA, NA, NA
  ),
  padj = NA,
  species = c(rep("Mouse", 5), "Human"),
  stringsAsFactors = FALSE
)

# Re-read from saved results properly
res_masld <- read.csv(file.path(results_dir, "DEG_MASLD_vs_Healthy.csv"))
res_hcc_nt <- read.csv(file.path(results_dir, "DEG_HCC_NonTumor_vs_Healthy.csv"))

get_kcnq <- function(df, id = "ENSMUSG00000101609_Kcnq1ot1") {
  row <- df[df$gene_id == id, ]
  if (nrow(row) > 0) return(c(row$log2FoldChange[1], row$padj[1]))
  return(c(NA, NA))
}
get_kcnq_h <- function(df, id = "ENSG00000269821_KCNQ1OT1") {
  row <- df[df$gene_id == id, ]
  if (nrow(row) > 0) return(c(row$log2FoldChange[1], row$padj[1]))
  return(c(NA, NA))
}

summary_data <- data.frame(
  comparison = c("MASLD\n(14w)", "MASH\n(20w)", "Fibrosis\n(32w)",
                 "HCC Non-Tumor\n(44-56w)", "HCC Tumor\n(44-56w)",
                 "Human\nNAFLD"),
  log2FC = c(get_kcnq(res_masld)[1], get_kcnq(res_mash)[1],
             get_kcnq(res_fib)[1], get_kcnq(res_hcc_nt)[1],
             get_kcnq(res_hcc)[1], get_kcnq_h(res_human)[1]),
  padj = c(get_kcnq(res_masld)[2], get_kcnq(res_mash)[2],
           get_kcnq(res_fib)[2], get_kcnq(res_hcc_nt)[2],
           get_kcnq(res_hcc)[2], get_kcnq_h(res_human)[2]),
  species = c(rep("Mouse (GSE246221)", 5), "Human (GSE246221)"),
  stringsAsFactors = FALSE
)

summary_data$comparison <- factor(summary_data$comparison,
  levels = c("MASLD\n(14w)", "MASH\n(20w)", "Fibrosis\n(32w)",
             "HCC Non-Tumor\n(44-56w)", "HCC Tumor\n(44-56w)",
             "Human\nNAFLD"))

summary_data$sig <- ifelse(!is.na(summary_data$padj) & summary_data$padj < 0.05, "*", "")
summary_data$sig_label <- ifelse(summary_data$sig == "*",
  sprintf("%.2f*", summary_data$log2FC),
  sprintf("%.2f", summary_data$log2FC))

p_summary <- ggplot(summary_data, aes(x = comparison, y = log2FC, fill = species)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sig_label),
            vjust = ifelse(summary_data$log2FC >= 0, -0.5, 1.5), size = 3.5) +
  scale_fill_manual(values = c("Mouse (GSE246221)" = "#4393C3", "Human (GSE246221)" = "#D6604D")) +
  labs(
    title = "Kcnq1ot1/KCNQ1OT1 Expression Changes Across Liver Disease Contexts",
    subtitle = "log2 fold change vs Healthy controls (* padj < 0.05, DESeq2)",
    x = "",
    y = "log2 Fold Change",
    fill = "Dataset"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 9),
    legend.position = "top"
  )

ggsave(file.path(fig_dir, "FigD_Integrated_Kcnq1ot1_summary.pdf"), p_summary,
       width = 9, height = 5)
ggsave(file.path(fig_dir, "FigD_Integrated_Kcnq1ot1_summary.png"), p_summary,
       width = 9, height = 5, dpi = 300)
cat("  FigD saved.\n")

# ============================================================
# Figure E: Model comparison schematic data
# ============================================================

# Original paper results (from the publication)
original_data <- data.frame(
  comparison = c("MASH\n(WD+CCl4)", "Fibrosis\n(CCl4)"),
  fold_change = c(0.6, 0.1),
  log2FC = c(log2(0.6), log2(0.1)),
  model = "Original Paper\n(WD+CCl4 / CCl4)",
  stringsAsFactors = FALSE
)

# Current study results
current_mouse <- data.frame(
  comparison = c("MASH\n(STZ+HFD)", "Fibrosis\n(STZ+HFD)"),
  fold_change = c(2^get_kcnq(res_mash)[1], 2^get_kcnq(res_fib)[1]),
  log2FC = c(get_kcnq(res_mash)[1], get_kcnq(res_fib)[1]),
  model = "This Study\n(STZ+HFD)",
  stringsAsFactors = FALSE
)

model_compare <- rbind(original_data, current_mouse)
model_compare$comparison <- factor(model_compare$comparison,
  levels = c("MASH\n(WD+CCl4)", "MASH\n(STZ+HFD)", "Fibrosis\n(CCl4)", "Fibrosis\n(STZ+HFD)"))

p_model <- ggplot(model_compare, aes(x = comparison, y = log2FC, fill = model)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sprintf("%.2f", fold_change)),
            vjust = ifelse(model_compare$log2FC >= 0, -0.5, 1.5), size = 3.5) +
  scale_fill_manual(values = c(
    "Original Paper\n(WD+CCl4 / CCl4)" = "#EF8A62",
    "This Study\n(STZ+HFD)" = "#67A9CF"
  )) +
  labs(
    title = "Kcnq1ot1 Response: Model-Dependent Regulation",
    subtitle = "Chemical injury (CCl4) vs Metabolic model (STZ+HFD)",
    x = "",
    y = "log2 Fold Change (vs Control)",
    fill = "Study",
    caption = "Values on bars = linear fold change"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "top"
  )

ggsave(file.path(fig_dir, "FigE_Model_comparison.pdf"), p_model, width = 7, height = 5)
ggsave(file.path(fig_dir, "FigE_Model_comparison.png"), p_model, width = 7, height = 5, dpi = 300)
cat("  FigE saved.\n")

cat("\nStep 3 (Publication figures) complete!\n")
cat("\nGenerated figures:\n")
cat("  FigA: Volcano plot - Mouse MASH vs Healthy\n")
cat("  FigB: Volcano plot - Mouse HCC Tumor vs Healthy\n")
cat("  FigC: Volcano plot - Human NAFLD vs Healthy\n")
cat("  FigD: Integrated Kcnq1ot1/KCNQ1OT1 summary across contexts\n")
cat("  FigE: Model comparison (Original paper vs This study)\n")
cat("  (plus previously generated Fig1-9)\n")

#!/usr/bin/env Rscript
# ============================================================
# 06_main_figures.R
# Multi-panel Main Figures for publication (Fig. 1-6)
# ============================================================

library(ggplot2)
library(edgeR)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(patchwork)
library(cowplot)

data_dir   <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/data"
results_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/results"
fig_dir    <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/figures/main"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Common theme
# ============================================================
theme_pub <- theme_bw(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    axis.title    = element_text(size = 9),
    axis.text     = element_text(size = 8),
    legend.title  = element_text(size = 9),
    legend.text   = element_text(size = 8),
    strip.text    = element_text(size = 9),
    plot.tag      = element_text(face = "bold", size = 14)
  )

stage_colors <- c(
  "Healthy" = "#2166AC", "Acute_STZ" = "#67A9CF", "MASLD" = "#D1E5F0",
  "MASH" = "#FDDBC7", "Fibrosis" = "#EF8A62",
  "HCC_NonTumor" = "#B2182B", "HCC_Tumor" = "#67001F"
)

# ============================================================
# Load all data
# ============================================================
cat("Loading data...\n")
gse246 <- read.delim(gzfile(file.path(data_dir, "GSE246221_rawcounts_allsamples_mouse.txt.gz")),
                     row.names = 1, check.names = FALSE)
human <- read.delim(gzfile(file.path(data_dir, "GSE246221_rawcounts_allsamples_human.txt.gz")),
                    row.names = 1, check.names = FALSE)

# Load DEG results
res_masld  <- read.csv(file.path(results_dir, "DEG_MASLD_vs_Healthy.csv"))
res_mash   <- read.csv(file.path(results_dir, "DEG_MASH_vs_Healthy.csv"))
res_fib    <- read.csv(file.path(results_dir, "DEG_Fibrosis_vs_Healthy.csv"))
res_hcc_nt <- read.csv(file.path(results_dir, "DEG_HCC_NonTumor_vs_Healthy.csv"))
res_hcc    <- read.csv(file.path(results_dir, "DEG_HCC_Tumor_vs_Healthy.csv"))
res_human  <- read.csv(file.path(results_dir, "Human_NAFLD_vs_Healthy_DESeq2.csv"))
cor_df     <- read.csv(file.path(results_dir, "Kcnq1ot1_correlations_all.csv"))
mito_cor   <- read.csv(file.path(results_dir, "Kcnq1ot1_mitochondrial_correlations.csv"))
sev_cor    <- read.csv(file.path(results_dir, "KCNQ1OT1_severity_marker_correlations.csv"))

# Normalize mouse Batch1
batch1_cols <- grep("^Batch1_", colnames(gse246), value = TRUE)
gse246_b1 <- gse246[, batch1_cols]
y246 <- DGEList(counts = round(gse246_b1))
y246 <- calcNormFactors(y246)
cpm246 <- cpm(y246, log = TRUE)

# Stage assignment
stage <- rep(NA, length(batch1_cols))
stage[grepl("Control07w", batch1_cols)] <- "Healthy"
stage[grepl("STZSCD08w", batch1_cols)]  <- "Acute_STZ"
stage[grepl("STZHFD14w", batch1_cols)]  <- "MASLD"
stage[grepl("STZHFD20w", batch1_cols)]  <- "MASH"
stage[grepl("STZHFD32w", batch1_cols)]  <- "Fibrosis"
stage[grepl("STZHFD44w|STZHFD50w|STZHFD56w", batch1_cols) & !grepl("_T\\d+$", batch1_cols)] <- "HCC_NonTumor"
stage[grepl("_T\\d+$", batch1_cols)] <- "HCC_Tumor"
names(stage) <- batch1_cols

kcnq_expr <- cpm246["ENSMUSG00000101609_Kcnq1ot1", ]

# Normalize human
y_h <- DGEList(counts = round(human))
y_h <- calcNormFactors(y_h)
cpm_h <- cpm(y_h, log = TRUE)
kcnq_h <- cpm_h["ENSG00000269821_KCNQ1OT1", ]

# ============================================================
# Figure 1: Kcnq1ot1 across disease progression (A: boxplot, B: trajectory)
# ============================================================
cat("Creating Figure 1...\n")

# Panel A: Boxplot by stage
df1a <- data.frame(
  logCPM = kcnq_expr,
  stage = factor(stage[names(kcnq_expr)],
    levels = c("Healthy", "Acute_STZ", "MASLD", "MASH", "Fibrosis", "HCC_NonTumor", "HCC_Tumor"))
)

p1a <- ggplot(df1a, aes(x = stage, y = logCPM, fill = stage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = stage_colors) +
  labs(title = "Kcnq1ot1 expression across disease stages",
       subtitle = "Mouse STZ+HFD model (GSE246221 Batch1, n=52)",
       x = "", y = "Kcnq1ot1 log2(CPM)") +
  theme_pub + theme(legend.position = "none",
                    axis.text.x = element_text(angle = 35, hjust = 1, size = 7))

# Panel B: Trajectory (mean +/- SE by timepoint)
timepoint <- rep(NA, length(batch1_cols))
timepoint[grepl("Control07w", batch1_cols)] <- 7
timepoint[grepl("STZSCD08w", batch1_cols)]  <- 8
timepoint[grepl("STZHFD14w", batch1_cols)]  <- 14
timepoint[grepl("STZHFD20w", batch1_cols)]  <- 20
timepoint[grepl("STZHFD32w", batch1_cols)]  <- 32
timepoint[grepl("STZHFD44w", batch1_cols)]  <- 44
timepoint[grepl("STZHFD50w", batch1_cols)]  <- 50
timepoint[grepl("STZHFD56w", batch1_cols)]  <- 56
names(timepoint) <- batch1_cols

is_tumor <- grepl("_T\\d+$", batch1_cols)
names(is_tumor) <- batch1_cols

df1b <- data.frame(
  logCPM = kcnq_expr,
  week = timepoint[names(kcnq_expr)],
  tissue = ifelse(is_tumor[names(kcnq_expr)], "Tumor", "Non-Tumor/Liver"),
  stringsAsFactors = FALSE
)
df1b <- df1b[!is.na(df1b$week), ]

# Compute mean/se
agg <- aggregate(logCPM ~ week + tissue, data = df1b, FUN = function(x)
  c(mean = mean(x), se = sd(x) / sqrt(length(x))))
agg <- data.frame(week = agg$week, tissue = agg$tissue,
                  mean = agg$logCPM[, "mean"], se = agg$logCPM[, "se"])

p1b <- ggplot(agg, aes(x = week, y = mean, color = tissue, group = tissue)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 1.5) +
  geom_vline(xintercept = 8, linetype = "dotted", color = "grey50") +
  annotate("text", x = 9, y = max(agg$mean) + 0.15, label = "STZ", size = 2.5, color = "grey50") +
  annotate("text", x = 15, y = min(agg$mean) - 0.15, label = "HFD start", size = 2.5, color = "grey50") +
  scale_color_manual(values = c("Non-Tumor/Liver" = "#4393C3", "Tumor" = "#B2182B")) +
  labs(title = "Kcnq1ot1 expression trajectory",
       subtitle = "Mean \u00b1 SE by timepoint",
       x = "Weeks", y = "Kcnq1ot1 log2(CPM)", color = "") +
  theme_pub + theme(legend.position = c(0.3, 0.85),
                    legend.background = element_rect(fill = alpha("white", 0.8)))

fig1 <- (p1a + labs(tag = "A")) + (p1b + labs(tag = "B")) + plot_layout(widths = c(1.2, 1))
ggsave(file.path(fig_dir, "Fig1_Kcnq1ot1_progression.pdf"), fig1, width = 12, height = 5)
ggsave(file.path(fig_dir, "Fig1_Kcnq1ot1_progression.png"), fig1, width = 12, height = 5, dpi = 300)
cat("  Fig1 saved.\n")

# ============================================================
# Figure 2: Volcano plots (A: MASH, B: HCC Tumor, C: Human NAFLD)
# ============================================================
cat("Creating Figure 2...\n")

prep_volcano <- function(df) {
  rownames(df) <- df$gene_id
  df$gene_label <- df$gene_symbol
  return(df)
}

mouse_lnc <- c("ENSMUSG00000101609_Kcnq1ot1", "ENSMUSG00000092274_Neat1",
               "ENSMUSG00000053332_Gas5")
human_lnc <- c("ENSG00000269821_KCNQ1OT1", "ENSG00000245532_NEAT1",
               "ENSG00000234741_GAS5")

res_mash_v  <- prep_volcano(res_mash)
res_hcc_v   <- prep_volcano(res_hcc)
res_human_v <- prep_volcano(res_human)

p2a <- EnhancedVolcano(res_mash_v,
  lab = res_mash_v$gene_label,
  selectLab = gsub("^ENSMUSG\\d+_", "", mouse_lnc),
  x = "log2FoldChange", y = "padj",
  title = "MASH vs Healthy", subtitle = "(Mouse STZ+HFD, 20w)",
  pCutoff = 0.05, FCcutoff = 1,
  pointSize = 1.2, labSize = 3,
  col = c("grey80", "#2166AC", "#D6604D", "#B2182B"),
  drawConnectors = TRUE, widthConnectors = 0.4, maxoverlapsConnectors = 20,
  legendPosition = "none"
) + theme_pub + theme(plot.title = element_text(size = 10))

p2b <- EnhancedVolcano(res_hcc_v,
  lab = res_hcc_v$gene_label,
  selectLab = gsub("^ENSMUSG\\d+_", "", mouse_lnc),
  x = "log2FoldChange", y = "padj",
  title = "HCC Tumor vs Healthy", subtitle = "(Mouse STZ+HFD, 44-56w)",
  pCutoff = 0.05, FCcutoff = 1,
  pointSize = 1.2, labSize = 3,
  col = c("grey80", "#2166AC", "#D6604D", "#B2182B"),
  drawConnectors = TRUE, widthConnectors = 0.4, maxoverlapsConnectors = 20,
  legendPosition = "none"
) + theme_pub + theme(plot.title = element_text(size = 10))

p2c <- EnhancedVolcano(res_human_v,
  lab = res_human_v$gene_symbol,
  selectLab = gsub("^ENSG\\d+_", "", human_lnc),
  x = "log2FoldChange", y = "padj",
  title = "NAFLD vs Healthy", subtitle = "(Human, 28 NAFLD + 4 Healthy)",
  pCutoff = 0.05, FCcutoff = 1,
  pointSize = 1.2, labSize = 3,
  col = c("grey80", "#2166AC", "#D6604D", "#B2182B"),
  drawConnectors = TRUE, widthConnectors = 0.4, maxoverlapsConnectors = 20,
  legendPosition = "none"
) + theme_pub + theme(plot.title = element_text(size = 10))

fig2 <- (p2a + labs(tag = "A")) + (p2b + labs(tag = "B")) + (p2c + labs(tag = "C")) +
  plot_layout(ncol = 3)
ggsave(file.path(fig_dir, "Fig2_Volcano_panels.pdf"), fig2, width = 16, height = 5.5)
ggsave(file.path(fig_dir, "Fig2_Volcano_panels.png"), fig2, width = 16, height = 5.5, dpi = 300)
cat("  Fig2 saved.\n")

# ============================================================
# Figure 3: KEGG pathway enrichment
# ============================================================
cat("Creating Figure 3...\n")

# Reload KEGG data for dotplot
library(clusterProfiler)
library(org.Mm.eg.db)

neg_genes <- cor_df$gene_symbol[cor_df$correlation < -0.5]
neg_entrez <- bitr(neg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
all_entrez <- bitr(cor_df$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

kegg_neg <- enrichKEGG(
  gene = neg_entrez$ENTREZID,
  universe = all_entrez$ENTREZID,
  organism = "mmu", pAdjustMethod = "BH",
  pvalueCutoff = 0.05, qvalueCutoff = 0.1
)

p3 <- dotplot(kegg_neg, showCategory = 15,
  title = "KEGG Pathways Enriched Among Genes Negatively\nCorrelated with Kcnq1ot1 (r < -0.5)") +
  theme_pub +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 11))

ggsave(file.path(fig_dir, "Fig3_KEGG_enrichment.pdf"), p3, width = 9, height = 7)
ggsave(file.path(fig_dir, "Fig3_KEGG_enrichment.png"), p3, width = 9, height = 7, dpi = 300)
cat("  Fig3 saved.\n")

# ============================================================
# Figure 4: Mitochondrial analysis (A: complex barplot, B: scatter 2x2)
# ============================================================
cat("Creating Figure 4...\n")

# Panel A: Complex-level mean correlation
mito_cor$complex <- factor(mito_cor$complex,
  levels = c("Complex_I", "Complex_II", "Complex_III", "Complex_IV", "Complex_V",
             "TCA_cycle", "FAO"))

complex_means <- aggregate(correlation ~ complex, data = mito_cor, mean)
complex_means$se <- aggregate(correlation ~ complex, data = mito_cor,
  function(x) sd(x)/sqrt(length(x)))$correlation
complex_means$n <- aggregate(correlation ~ complex, data = mito_cor, length)$correlation

p4a <- ggplot(complex_means, aes(x = complex, y = correlation, fill = correlation)) +
  geom_col(alpha = 0.85, width = 0.65) +
  geom_errorbar(aes(ymin = correlation - se, ymax = correlation + se), width = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dotted", color = "red", alpha = 0.6) +
  geom_text(aes(label = sprintf("n=%d", n)), y = -0.05, size = 2.5, color = "grey40") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Mean correlation with Kcnq1ot1\nby mitochondrial complex",
       x = "", y = "Mean Pearson r") +
  theme_pub + theme(legend.position = "none",
                    axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

# Panel B: Scatter plots (2x2)
top_mito <- c("Sdha", "Sdhb", "Ndufa9", "Atp5a1")
scatter_data <- data.frame()
for (g in top_mito) {
  g_row <- grep(paste0("_", g, "$"), rownames(cpm246), value = TRUE)
  if (length(g_row) > 0) {
    r_val <- cor(kcnq_expr, cpm246[g_row[1], ])
    scatter_data <- rbind(scatter_data, data.frame(
      Kcnq1ot1 = kcnq_expr, mito_gene = cpm246[g_row[1], ],
      gene = paste0(g, " (r=", round(r_val, 3), ")"),
      stage = stage[names(kcnq_expr)], stringsAsFactors = FALSE
    ))
  }
}

p4b <- ggplot(scatter_data, aes(x = Kcnq1ot1, y = mito_gene, color = stage)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.4) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  scale_color_manual(values = stage_colors) +
  labs(title = "Kcnq1ot1 vs key mitochondrial genes",
       x = "Kcnq1ot1 log2(CPM)", y = "Mito gene log2(CPM)", color = "Stage") +
  theme_pub + theme(strip.text = element_text(face = "italic", size = 8),
                    legend.position = "bottom",
                    legend.key.size = unit(0.3, "cm"))

fig4 <- (p4a + labs(tag = "A")) + (p4b + labs(tag = "B")) +
  plot_layout(widths = c(0.8, 1.2))
ggsave(file.path(fig_dir, "Fig4_Mitochondria.pdf"), fig4, width = 14, height = 6)
ggsave(file.path(fig_dir, "Fig4_Mitochondria.png"), fig4, width = 14, height = 6, dpi = 300)
cat("  Fig4 saved.\n")

# ============================================================
# Figure 5: Human KCNQ1OT1 (A: boxplot, B: severity markers)
# ============================================================
cat("Creating Figure 5...\n")

# Panel A: Boxplot
meta_human <- data.frame(
  sample = colnames(human),
  group = ifelse(grepl("Healthy", colnames(human)), "Healthy", "NAFLD"),
  stringsAsFactors = FALSE
)
rownames(meta_human) <- meta_human$sample

# DESeq2 result for annotation
dds_h <- DESeqDataSetFromMatrix(round(human), meta_human, ~ group)
keep <- rowSums(counts(dds_h) >= 10) >= 3
dds_h <- dds_h[keep, ]
dds_h$group <- relevel(factor(dds_h$group), ref = "Healthy")
dds_h <- DESeq(dds_h)
res_h <- results(dds_h, contrast = c("group", "NAFLD", "Healthy"))

plot_h <- data.frame(
  sample = names(kcnq_h), logCPM = kcnq_h,
  group = meta_human[names(kcnq_h), "group"]
)

p5a <- ggplot(plot_h, aes(x = group, y = logCPM, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("Healthy" = "#2166AC", "NAFLD" = "#B2182B")) +
  annotate("text", x = 1.5, y = max(plot_h$logCPM) + 0.15,
           label = sprintf("log2FC = +%.2f\npadj = %.1e",
             as.numeric(res_h["ENSG00000269821_KCNQ1OT1", "log2FoldChange"]),
             as.numeric(res_h["ENSG00000269821_KCNQ1OT1", "padj"])),
           size = 3) +
  labs(title = "KCNQ1OT1 in Human NAFLD",
       subtitle = "GSE246221 (4 Healthy + 28 NAFLD)",
       x = "", y = "KCNQ1OT1 log2(CPM)") +
  theme_pub + theme(legend.position = "none")

# Panel B: Severity marker correlation
sev_cor$sig <- ifelse(sev_cor$pvalue < 0.05, "*", "")
sev_cor$category <- NA
sev_cor$category[sev_cor$gene %in% c("TNF", "IL1B", "IL6", "CCL2")] <- "Inflammation"
sev_cor$category[sev_cor$gene %in% c("COL1A1", "COL1A2", "COL3A1", "ACTA2", "TGFB1")] <- "Fibrosis"
sev_cor$category[sev_cor$gene %in% c("FASN", "SREBF1", "PPARG", "CD36")] <- "Steatosis"
sev_cor$category[sev_cor$gene %in% c("CYP2E1", "ALB")] <- "Hepatocyte"
sev_cor$category[sev_cor$gene %in% c("NDUFA9", "SDHA", "SDHB", "ATP5F1A", "UQCRC2")] <- "Mitochondria"
sev_cor$gene <- factor(sev_cor$gene, levels = rev(sev_cor$gene))

p5b <- ggplot(sev_cor, aes(x = correlation, y = gene, fill = category)) +
  geom_col(alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sig, x = correlation + sign(correlation) * 0.02),
            size = 5, vjust = 0.5) +
  scale_fill_manual(values = c(
    "Inflammation" = "#E41A1C", "Fibrosis" = "#FF7F00",
    "Steatosis" = "#984EA3", "Hepatocyte" = "#4DAF4A",
    "Mitochondria" = "#377EB8"
  )) +
  labs(title = "KCNQ1OT1 correlation with\ndisease severity markers",
       subtitle = "Human (n=32, Pearson)",
       x = "Pearson r", y = "", fill = "Category") +
  theme_pub + theme(legend.position = "right",
                    legend.key.size = unit(0.35, "cm"))

fig5 <- (p5a + labs(tag = "A")) + (p5b + labs(tag = "B")) +
  plot_layout(widths = c(0.8, 1.2))
ggsave(file.path(fig_dir, "Fig5_Human_KCNQ1OT1.pdf"), fig5, width = 12, height = 5.5)
ggsave(file.path(fig_dir, "Fig5_Human_KCNQ1OT1.png"), fig5, width = 12, height = 5.5, dpi = 300)
cat("  Fig5 saved.\n")

# ============================================================
# Figure 6: Integrated summary (A: all comparisons, B: model comparison)
# ============================================================
cat("Creating Figure 6...\n")

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

# Panel A: All comparisons
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
  species = c(rep("Mouse (STZ+HFD)", 5), "Human (NAFLD)"),
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

p6a <- ggplot(summary_data, aes(x = comparison, y = log2FC, fill = species)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sig_label),
            vjust = ifelse(summary_data$log2FC >= 0, -0.5, 1.5), size = 3) +
  scale_fill_manual(values = c("Mouse (STZ+HFD)" = "#4393C3", "Human (NAFLD)" = "#D6604D")) +
  labs(title = "Kcnq1ot1/KCNQ1OT1 log2FC across all comparisons",
       subtitle = "vs Healthy (* padj < 0.05, DESeq2)",
       x = "", y = "log2 Fold Change", fill = "") +
  theme_pub + theme(legend.position = "top", axis.text.x = element_text(size = 7))

# Panel B: Model comparison
original_data <- data.frame(
  comparison = c("MASH\n(WD+CCl4)", "Fibrosis\n(CCl4)"),
  fold_change = c(0.6, 0.1),
  log2FC = c(log2(0.6), log2(0.1)),
  model = "Original (CCl4-based)",
  stringsAsFactors = FALSE
)
current_mouse <- data.frame(
  comparison = c("MASH\n(STZ+HFD)", "Fibrosis\n(STZ+HFD)"),
  fold_change = c(2^get_kcnq(res_mash)[1], 2^get_kcnq(res_fib)[1]),
  log2FC = c(get_kcnq(res_mash)[1], get_kcnq(res_fib)[1]),
  model = "This Study (Metabolic)",
  stringsAsFactors = FALSE
)
model_compare <- rbind(original_data, current_mouse)
model_compare$comparison <- factor(model_compare$comparison,
  levels = c("MASH\n(WD+CCl4)", "MASH\n(STZ+HFD)", "Fibrosis\n(CCl4)", "Fibrosis\n(STZ+HFD)"))

p6b <- ggplot(model_compare, aes(x = comparison, y = log2FC, fill = model)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sprintf("%.2fx", fold_change)),
            vjust = ifelse(model_compare$log2FC >= 0, -0.5, 1.5), size = 3) +
  scale_fill_manual(values = c("Original (CCl4-based)" = "#EF8A62", "This Study (Metabolic)" = "#67A9CF")) +
  labs(title = "Model-dependent Kcnq1ot1 regulation",
       subtitle = "Chemical injury vs Metabolic model",
       x = "", y = "log2 Fold Change", fill = "") +
  theme_pub + theme(legend.position = "top")

fig6 <- (p6a + labs(tag = "A")) + (p6b + labs(tag = "B")) +
  plot_layout(widths = c(1.3, 1))
ggsave(file.path(fig_dir, "Fig6_Integrated_summary.pdf"), fig6, width = 13, height = 5.5)
ggsave(file.path(fig_dir, "Fig6_Integrated_summary.png"), fig6, width = 13, height = 5.5, dpi = 300)
cat("  Fig6 saved.\n")

cat("\n=== All 6 Main Figures generated ===\n")
cat("  Fig1: Kcnq1ot1 disease progression (A: boxplot, B: trajectory)\n")
cat("  Fig2: Volcano plots (A: MASH, B: HCC Tumor, C: Human NAFLD)\n")
cat("  Fig3: KEGG enrichment dotplot\n")
cat("  Fig4: Mitochondrial analysis (A: complex barplot, B: scatter)\n")
cat("  Fig5: Human KCNQ1OT1 (A: boxplot, B: severity markers)\n")
cat("  Fig6: Integrated summary (A: all comparisons, B: model comparison)\n")

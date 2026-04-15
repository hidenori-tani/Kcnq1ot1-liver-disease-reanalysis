#!/usr/bin/env Rscript
# ============================================================
# 02_Kcnq1ot1_comprehensive.R
# Comprehensive analysis of Kcnq1ot1 in liver disease progression
# Using GSE246221 (STZ+HFD) and GSE207855 (CCl4)
# ============================================================

library(ggplot2)
library(DESeq2)
library(edgeR)
library(pheatmap)
library(clusterProfiler)
library(sva)

data_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/data"
results_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/results"
fig_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/figures"

# ============================================================
# 1. Load data
# ============================================================
cat("Loading datasets...\n")
gse246 <- read.delim(
  gzfile(file.path(data_dir, "GSE246221_rawcounts_allsamples_mouse.txt.gz")),
  row.names = 1, check.names = FALSE
)

gse207 <- read.delim(
  gzfile(file.path(data_dir, "GSE207855_count.tsv.gz")),
  row.names = 1, check.names = FALSE
)

# ============================================================
# 2. GSE246221: Focus on Batch1 disease progression
# ============================================================
cat("\n=== GSE246221: Disease Progression Analysis ===\n")

# Select Batch1 samples only (main time course)
batch1_cols <- grep("^Batch1_", colnames(gse246), value = TRUE)
gse246_b1 <- gse246[, batch1_cols]

# Create metadata
meta_b1 <- data.frame(
  sample = batch1_cols,
  stringsAsFactors = FALSE
)
rownames(meta_b1) <- batch1_cols

meta_b1$is_tumor <- grepl("_T\\d+$", batch1_cols)
meta_b1$stage <- NA
meta_b1$stage[grepl("Control07w", batch1_cols)] <- "Healthy"
meta_b1$stage[grepl("STZSCD08w", batch1_cols)] <- "Acute_STZ"
meta_b1$stage[grepl("STZHFD14w", batch1_cols)] <- "MASLD"
meta_b1$stage[grepl("STZHFD20w", batch1_cols)] <- "MASH"
meta_b1$stage[grepl("STZHFD32w", batch1_cols)] <- "Fibrosis"
meta_b1$stage[grepl("STZHFD44w|STZHFD50w|STZHFD56w", batch1_cols) & !meta_b1$is_tumor] <- "HCC_NonTumor"
meta_b1$stage[grepl("STZHFD44w|STZHFD50w|STZHFD56w", batch1_cols) & meta_b1$is_tumor] <- "HCC_Tumor"

meta_b1$stage <- factor(meta_b1$stage,
  levels = c("Healthy", "Acute_STZ", "MASLD", "MASH", "Fibrosis", "HCC_NonTumor", "HCC_Tumor"))

# Extract gene symbols from EnsemblID_Symbol format
gene_symbols <- gsub("^ENSMUSG\\d+_", "", rownames(gse246_b1))
ensembl_ids <- gsub("_.*$", "", rownames(gse246_b1))

# ============================================================
# 3. DESeq2 analysis: Each stage vs Healthy
# ============================================================
cat("Running DESeq2 on GSE246221 Batch1...\n")

# Round fractional counts
gse246_b1_int <- round(gse246_b1)

# Remove genes with very low counts
keep <- rowSums(gse246_b1_int >= 10) >= 3
gse246_b1_filt <- gse246_b1_int[keep, ]
cat(sprintf("  Genes after filtering: %d (from %d)\n", nrow(gse246_b1_filt), nrow(gse246_b1_int)))

dds246 <- DESeqDataSetFromMatrix(
  countData = gse246_b1_filt,
  colData = meta_b1,
  design = ~ stage
)
dds246 <- DESeq(dds246)

# Extract results for each comparison vs Healthy
comparisons <- c("MASLD", "MASH", "Fibrosis", "HCC_NonTumor", "HCC_Tumor")
deg_results <- list()
kcnq1ot1_row <- "ENSMUSG00000101609_Kcnq1ot1"

cat("\nKcnq1ot1 differential expression (vs Healthy):\n")
cat(sprintf("%-15s %10s %10s %12s %12s\n", "Comparison", "log2FC", "lfcSE", "pvalue", "padj"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (comp in comparisons) {
  res <- results(dds246, contrast = c("stage", comp, "Healthy"))
  deg_results[[comp]] <- as.data.frame(res)
  deg_results[[comp]]$gene_id <- rownames(res)
  deg_results[[comp]]$gene_symbol <- gsub("^ENSMUSG\\d+_", "", rownames(res))

  if (kcnq1ot1_row %in% rownames(res)) {
    r <- res[kcnq1ot1_row, ]
    cat(sprintf("%-15s %10.4f %10.4f %12.2e %12.2e\n",
                comp, r$log2FoldChange, r$lfcSE, r$pvalue, r$padj))
  }
}

# Save all DEG results
for (comp in names(deg_results)) {
  write.csv(deg_results[[comp]],
    file.path(results_dir, paste0("DEG_", comp, "_vs_Healthy.csv")),
    row.names = FALSE)
}

# ============================================================
# 4. Kcnq1ot1 correlated genes (MASH stage)
# ============================================================
cat("\n=== Kcnq1ot1 Co-expression Analysis ===\n")

# Use TMM-normalized CPM across all Batch1 samples
y246 <- DGEList(counts = gse246_b1_int)
y246 <- calcNormFactors(y246)
cpm246 <- cpm(y246, log = TRUE)

# Filter low-expression genes for correlation
keep_cor <- rowMeans(cpm246) > 1
cpm246_filt <- cpm246[keep_cor, ]
cat(sprintf("  Genes for correlation analysis: %d\n", nrow(cpm246_filt)))

# Compute correlation with Kcnq1ot1
kcnq1ot1_expr <- cpm246_filt[kcnq1ot1_row, ]
cors <- apply(cpm246_filt, 1, function(x) cor(x, kcnq1ot1_expr, method = "pearson"))
cors <- cors[!is.na(cors)]
cors <- cors[names(cors) != kcnq1ot1_row]  # remove self

# Top positively and negatively correlated
top_pos <- sort(cors, decreasing = TRUE)[1:100]
top_neg <- sort(cors)[1:100]

top_pos_symbols <- gsub("^ENSMUSG\\d+_", "", names(top_pos))
top_neg_symbols <- gsub("^ENSMUSG\\d+_", "", names(top_neg))

cat(sprintf("\nTop 10 positively correlated genes with Kcnq1ot1:\n"))
for (i in 1:10) {
  cat(sprintf("  %s: r = %.3f\n", top_pos_symbols[i], top_pos[i]))
}

cat(sprintf("\nTop 10 negatively correlated genes with Kcnq1ot1:\n"))
for (i in 1:10) {
  cat(sprintf("  %s: r = %.3f\n", top_neg_symbols[i], top_neg[i]))
}

# Save correlation results
cor_df <- data.frame(
  gene_id = names(cors),
  gene_symbol = gsub("^ENSMUSG\\d+_", "", names(cors)),
  correlation = cors,
  abs_correlation = abs(cors)
)
cor_df <- cor_df[order(-cor_df$abs_correlation), ]
write.csv(cor_df, file.path(results_dir, "Kcnq1ot1_correlations_all.csv"), row.names = FALSE)

# ============================================================
# 5. GO/KEGG Enrichment of Kcnq1ot1-correlated genes
# ============================================================
cat("\n=== Pathway Enrichment of Kcnq1ot1-Correlated Genes ===\n")

# Positively correlated genes (r > 0.5)
pos_genes <- gsub("^ENSMUSG\\d+_", "", names(cors[cors > 0.5]))
neg_genes <- gsub("^ENSMUSG\\d+_", "", names(cors[cors < -0.5]))
cat(sprintf("  Positively correlated (r>0.5): %d genes\n", length(pos_genes)))
cat(sprintf("  Negatively correlated (r<-0.5): %d genes\n", length(neg_genes)))

# All gene symbols as universe
all_symbols <- gsub("^ENSMUSG\\d+_", "", rownames(cpm246_filt))

# GO enrichment for positively correlated genes
if (length(pos_genes) >= 10) {
  ego_pos <- enrichGO(
    gene = pos_genes,
    universe = all_symbols,
    OrgDb = "org.Mm.eg.db",
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )
  if (nrow(ego_pos@result[ego_pos@result$p.adjust < 0.05, ]) > 0) {
    cat("\nTop GO BP terms (positively correlated with Kcnq1ot1):\n")
    top_go_pos <- head(ego_pos@result[ego_pos@result$p.adjust < 0.05, ], 15)
    for (i in 1:min(15, nrow(top_go_pos))) {
      cat(sprintf("  %s (p.adj=%.2e, count=%d)\n",
                  top_go_pos$Description[i], top_go_pos$p.adjust[i], top_go_pos$Count[i]))
    }
    write.csv(ego_pos@result, file.path(results_dir, "GO_BP_Kcnq1ot1_positive_corr.csv"), row.names = FALSE)
  } else {
    cat("  No significant GO terms for positively correlated genes.\n")
  }
}

# GO enrichment for negatively correlated genes
if (length(neg_genes) >= 10) {
  ego_neg <- enrichGO(
    gene = neg_genes,
    universe = all_symbols,
    OrgDb = "org.Mm.eg.db",
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )
  if (nrow(ego_neg@result[ego_neg@result$p.adjust < 0.05, ]) > 0) {
    cat("\nTop GO BP terms (negatively correlated with Kcnq1ot1):\n")
    top_go_neg <- head(ego_neg@result[ego_neg@result$p.adjust < 0.05, ], 15)
    for (i in 1:min(15, nrow(top_go_neg))) {
      cat(sprintf("  %s (p.adj=%.2e, count=%d)\n",
                  top_go_neg$Description[i], top_go_neg$p.adjust[i], top_go_neg$Count[i]))
    }
    write.csv(ego_neg@result, file.path(results_dir, "GO_BP_Kcnq1ot1_negative_corr.csv"), row.names = FALSE)
  } else {
    cat("  No significant GO terms for negatively correlated genes.\n")
  }
}

# ============================================================
# 6. DEG analysis for MASH vs Healthy (full result)
# ============================================================
cat("\n=== MASH vs Healthy: Full DEG Summary ===\n")
res_mash <- deg_results[["MASH"]]
sig_mash <- res_mash[!is.na(res_mash$padj) & res_mash$padj < 0.05, ]
sig_mash_up <- sig_mash[sig_mash$log2FoldChange > 0, ]
sig_mash_down <- sig_mash[sig_mash$log2FoldChange < 0, ]
cat(sprintf("  Significant DEGs (padj<0.05): %d total\n", nrow(sig_mash)))
cat(sprintf("    Upregulated: %d\n", nrow(sig_mash_up)))
cat(sprintf("    Downregulated: %d\n", nrow(sig_mash_down)))

# Check if other known lncRNAs are differentially expressed
known_lnc <- c("Neat1", "Gas5", "Mir17hg", "Snhg4", "Snhg8", "Snhg17", "Peg13")
cat("\nKnown lncRNAs from original paper (MASH vs Healthy):\n")
cat(sprintf("%-12s %10s %10s %12s %12s\n", "Gene", "log2FC", "baseMean", "pvalue", "padj"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (lnc in known_lnc) {
  lnc_row <- grep(paste0("_", lnc, "$"), res_mash$gene_id, value = TRUE)
  if (length(lnc_row) > 0) {
    r <- res_mash[res_mash$gene_id == lnc_row[1], ]
    padj_val <- ifelse(is.na(r$padj), "NA", sprintf("%.2e", r$padj))
    pval_val <- ifelse(is.na(r$pvalue), "NA", sprintf("%.2e", r$pvalue))
    cat(sprintf("%-12s %10.4f %10.1f %12s %12s\n",
                lnc, r$log2FoldChange, r$baseMean, pval_val, padj_val))
  }
}

# Kcnq1ot1 in the same table
r_kcnq <- res_mash[res_mash$gene_id == "ENSMUSG00000101609_Kcnq1ot1", ]
cat(sprintf("%-12s %10.4f %10.1f %12.2e %12.2e\n",
            "Kcnq1ot1", r_kcnq$log2FoldChange, r_kcnq$baseMean,
            r_kcnq$pvalue, r_kcnq$padj))

# ============================================================
# 7. DEG analysis for Fibrosis vs Healthy
# ============================================================
cat("\n=== Fibrosis vs Healthy: Full DEG Summary ===\n")
res_fib <- deg_results[["Fibrosis"]]
sig_fib <- res_fib[!is.na(res_fib$padj) & res_fib$padj < 0.05, ]
sig_fib_up <- sig_fib[sig_fib$log2FoldChange > 0, ]
sig_fib_down <- sig_fib[sig_fib$log2FoldChange < 0, ]
cat(sprintf("  Significant DEGs (padj<0.05): %d total\n", nrow(sig_fib)))
cat(sprintf("    Upregulated: %d\n", nrow(sig_fib_up)))
cat(sprintf("    Downregulated: %d\n", nrow(sig_fib_down)))

# Known lncRNAs in Fibrosis
cat("\nKnown lncRNAs (Fibrosis vs Healthy):\n")
cat(sprintf("%-12s %10s %10s %12s %12s\n", "Gene", "log2FC", "baseMean", "pvalue", "padj"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (lnc in c(known_lnc, "Kcnq1ot1")) {
  lnc_row <- grep(paste0("_", lnc, "$"), res_fib$gene_id, value = TRUE)
  if (length(lnc_row) > 0) {
    r <- res_fib[res_fib$gene_id == lnc_row[1], ]
    padj_val <- ifelse(is.na(r$padj), "NA", sprintf("%.2e", r$padj))
    pval_val <- ifelse(is.na(r$pvalue), "NA", sprintf("%.2e", r$pvalue))
    cat(sprintf("%-12s %10.4f %10.1f %12s %12s\n",
                lnc, r$log2FoldChange, r$baseMean, pval_val, padj_val))
  }
}

# ============================================================
# 8. HCC Tumor vs Non-Tumor
# ============================================================
cat("\n=== HCC Tumor vs Non-Tumor: Full DEG Summary ===\n")
res_hcc <- deg_results[["HCC_Tumor"]]
sig_hcc <- res_hcc[!is.na(res_hcc$padj) & res_hcc$padj < 0.05, ]
cat(sprintf("  Significant DEGs (padj<0.05): %d total\n", nrow(sig_hcc)))
cat(sprintf("    Upregulated: %d\n", nrow(sig_hcc[sig_hcc$log2FoldChange > 0, ])))
cat(sprintf("    Downregulated: %d\n", nrow(sig_hcc[sig_hcc$log2FoldChange < 0, ])))

cat("\nKcnq1ot1 in HCC Tumor vs Healthy:\n")
r_hcc <- res_hcc[res_hcc$gene_id == "ENSMUSG00000101609_Kcnq1ot1", ]
cat(sprintf("  log2FC=%.4f, padj=%.2e\n", r_hcc$log2FoldChange, r_hcc$padj))

# ============================================================
# 9. GSE207855: Corrected analysis
# ============================================================
cat("\n=== GSE207855: Corrected CCl4 Analysis ===\n")

meta207 <- data.frame(
  sample = colnames(gse207),
  condition = ifelse(grepl("^CCl4", colnames(gse207)), "CCl4", "Oil"),
  stringsAsFactors = FALSE
)
rownames(meta207) <- meta207$sample

dds207 <- DESeqDataSetFromMatrix(
  countData = gse207,
  colData = meta207,
  design = ~ condition
)
dds207$condition <- relevel(dds207$condition, ref = "Oil")
dds207 <- DESeq(dds207)
res207 <- results(dds207, contrast = c("condition", "CCl4", "Oil"))

cat("\nKcnq1ot1 (Entrez 63830) in CCl4 vs Oil:\n")
print(res207["63830", ])

cat("\nRmst (Entrez 110333) in CCl4 vs Oil:\n")
if ("110333" %in% rownames(res207)) {
  print(res207["110333", ])
} else {
  cat("  Not in filtered results (too low expression)\n")
}

# Total DEGs in GSE207855
res207_df <- as.data.frame(res207)
sig207 <- res207_df[!is.na(res207_df$padj) & res207_df$padj < 0.05, ]
cat(sprintf("\nGSE207855 total DEGs (padj<0.05): %d\n", nrow(sig207)))
cat(sprintf("  Upregulated: %d\n", nrow(sig207[sig207$log2FoldChange > 0, ])))
cat(sprintf("  Downregulated: %d\n", nrow(sig207[sig207$log2FoldChange < 0, ])))

# ============================================================
# 10. Visualization: Comprehensive figure
# ============================================================
cat("\n=== Generating Comprehensive Figures ===\n")

# --- Figure 3: Heatmap of top DEGs in MASH ---
sig_mash_top <- sig_mash[order(sig_mash$padj), ]
top50_genes <- head(sig_mash_top$gene_id, 50)
top50_labels <- gsub("^ENSMUSG\\d+_", "", top50_genes)

# Add Kcnq1ot1 if not in top50
if (!"ENSMUSG00000101609_Kcnq1ot1" %in% top50_genes) {
  top50_genes <- c(top50_genes, "ENSMUSG00000101609_Kcnq1ot1")
  top50_labels <- c(top50_labels, "Kcnq1ot1")
}

# Z-score normalize
heatmap_data <- cpm246[top50_genes, ]
heatmap_z <- t(scale(t(heatmap_data)))
rownames(heatmap_z) <- gsub("^ENSMUSG\\d+_", "", rownames(heatmap_z))

# Annotation
ann_col <- data.frame(
  Stage = meta_b1$stage,
  row.names = meta_b1$sample
)
ann_colors <- list(
  Stage = c(
    "Healthy" = "#2166AC", "Acute_STZ" = "#67A9CF", "MASLD" = "#D1E5F0",
    "MASH" = "#FDDBC7", "Fibrosis" = "#EF8A62",
    "HCC_NonTumor" = "#B2182B", "HCC_Tumor" = "#67001F"
  )
)

pdf(file.path(fig_dir, "Fig3_Heatmap_top_DEGs_MASH.pdf"), width = 14, height = 12)
pheatmap(heatmap_z,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 7,
  main = "Top 50 DEGs (MASH vs Healthy) + Kcnq1ot1\nGSE246221 Batch1, Z-score normalized"
)
dev.off()

png(file.path(fig_dir, "Fig3_Heatmap_top_DEGs_MASH.png"), width = 1400, height = 1200, res = 150)
pheatmap(heatmap_z,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 7,
  main = "Top 50 DEGs (MASH vs Healthy) + Kcnq1ot1\nGSE246221 Batch1, Z-score normalized"
)
dev.off()

cat("  Fig3 saved.\n")

# --- Figure 4: Kcnq1ot1 time-course trend line ---
kcnq1ot1_tc <- data.frame(
  sample = names(cpm246[kcnq1ot1_row, ]),
  logCPM = cpm246[kcnq1ot1_row, ],
  stage = meta_b1[names(cpm246[kcnq1ot1_row, ]), "stage"]
)
kcnq1ot1_tc <- kcnq1ot1_tc[!is.na(kcnq1ot1_tc$stage), ]

# Numeric x for trend
stage_num <- c("Healthy" = 1, "Acute_STZ" = 2, "MASLD" = 3, "MASH" = 4,
               "Fibrosis" = 5, "HCC_NonTumor" = 6, "HCC_Tumor" = 7)
kcnq1ot1_tc$stage_num <- stage_num[as.character(kcnq1ot1_tc$stage)]

stage_labels2 <- c("1" = "Healthy\n(7w)", "2" = "Acute\n(8w)", "3" = "MASLD\n(14w)",
                   "4" = "MASH\n(20w)", "5" = "Fibrosis\n(32w)",
                   "6" = "HCC\nNon-Tumor", "7" = "HCC\nTumor")

stage_colors2 <- c(
  "Healthy" = "#2166AC", "Acute_STZ" = "#67A9CF", "MASLD" = "#D1E5F0",
  "MASH" = "#FDDBC7", "Fibrosis" = "#EF8A62",
  "HCC_NonTumor" = "#B2182B", "HCC_Tumor" = "#67001F"
)

# Summary stats per stage
stage_summary <- aggregate(logCPM ~ stage + stage_num, data = kcnq1ot1_tc,
  FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
stage_summary <- do.call(data.frame, stage_summary)
colnames(stage_summary) <- c("stage", "stage_num", "mean", "se")

p4 <- ggplot() +
  geom_line(data = stage_summary, aes(x = stage_num, y = mean),
            color = "gray40", linewidth = 0.8) +
  geom_errorbar(data = stage_summary,
    aes(x = stage_num, ymin = mean - se, ymax = mean + se),
    width = 0.2, color = "gray40") +
  geom_jitter(data = kcnq1ot1_tc,
    aes(x = stage_num, y = logCPM, color = stage),
    width = 0.1, size = 2, alpha = 0.7) +
  geom_point(data = stage_summary,
    aes(x = stage_num, y = mean), size = 3, shape = 18) +
  scale_color_manual(values = stage_colors2) +
  scale_x_continuous(breaks = 1:7, labels = stage_labels2) +
  labs(
    title = "Kcnq1ot1 Expression Trajectory in Liver Disease Progression",
    subtitle = "GSE246221: STZ+HFD model (Batch1 time course)",
    x = "Disease Stage",
    y = "log2(CPM)",
    color = "Stage"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "none"
  )

ggsave(file.path(fig_dir, "Fig4_Kcnq1ot1_trajectory.pdf"), p4, width = 9, height = 5)
ggsave(file.path(fig_dir, "Fig4_Kcnq1ot1_trajectory.png"), p4, width = 9, height = 5, dpi = 300)
cat("  Fig4 saved.\n")

# ============================================================
# 11. Statistical tests for Kcnq1ot1 pairwise comparisons
# ============================================================
cat("\n=== Pairwise t-tests for Kcnq1ot1 ===\n")

stages_to_test <- c("MASLD", "MASH", "Fibrosis", "HCC_NonTumor", "HCC_Tumor")
healthy_vals <- kcnq1ot1_tc$logCPM[kcnq1ot1_tc$stage == "Healthy"]

for (st in stages_to_test) {
  vals <- kcnq1ot1_tc$logCPM[kcnq1ot1_tc$stage == st]
  tt <- t.test(vals, healthy_vals)
  fc <- mean(vals) - mean(healthy_vals)  # difference in log2 scale
  cat(sprintf("  %s vs Healthy: delta_logCPM=%.3f, p=%.4f (n=%d vs %d)\n",
              st, fc, tt$p.value, length(vals), length(healthy_vals)))
}

cat("\nAll analyses complete!\n")

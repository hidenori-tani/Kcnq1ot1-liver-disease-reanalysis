#!/usr/bin/env Rscript
# ============================================================
# 03_KEGG_and_network.R
# Step 1: KEGG pathway enrichment + mitochondrial gene details
# ============================================================

library(ggplot2)
library(edgeR)
library(clusterProfiler)
library(org.Mm.eg.db)

data_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/data"
results_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/results"
fig_dir <- "/Users/tanihidenori/claude-work/research/lncRNA_liver_dry/figures"

# ============================================================
# 1. Load pre-computed correlation data
# ============================================================
cat("Loading correlation data...\n")
cor_df <- read.csv(file.path(results_dir, "Kcnq1ot1_correlations_all.csv"))

neg_genes <- cor_df$gene_symbol[cor_df$correlation < -0.5]
pos_genes <- cor_df$gene_symbol[cor_df$correlation > 0.5]
all_genes <- cor_df$gene_symbol

cat(sprintf("  Negatively correlated (r < -0.5): %d genes\n", length(neg_genes)))
cat(sprintf("  Positively correlated (r > 0.5): %d genes\n", length(pos_genes)))

# ============================================================
# 2. Convert gene symbols to Entrez IDs
# ============================================================
cat("\nConverting gene symbols to Entrez IDs...\n")

neg_entrez <- bitr(neg_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
pos_entrez <- bitr(pos_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
all_entrez <- bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)

cat(sprintf("  Mapped: neg=%d, pos=%d, universe=%d\n",
    nrow(neg_entrez), nrow(pos_entrez), nrow(all_entrez)))

# ============================================================
# 3. KEGG enrichment - Negatively correlated genes
# ============================================================
cat("\n=== KEGG Enrichment: Negatively Correlated Genes ===\n")

kegg_neg <- enrichKEGG(
  gene = neg_entrez$ENTREZID,
  universe = all_entrez$ENTREZID,
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

if (!is.null(kegg_neg) && nrow(kegg_neg@result[kegg_neg@result$p.adjust < 0.05, ]) > 0) {
  kegg_neg_sig <- kegg_neg@result[kegg_neg@result$p.adjust < 0.05, ]
  cat(sprintf("\nSignificant KEGG pathways: %d\n", nrow(kegg_neg_sig)))
  cat("\nTop 20 KEGG pathways (negatively correlated with Kcnq1ot1):\n")
  cat(sprintf("%-50s %8s %8s %6s\n", "Pathway", "p.adjust", "Count", "Ratio"))
  cat(paste(rep("-", 78), collapse = ""), "\n")
  for (i in 1:min(20, nrow(kegg_neg_sig))) {
    cat(sprintf("%-50s %8.2e %8d %6s\n",
        kegg_neg_sig$Description[i],
        kegg_neg_sig$p.adjust[i],
        kegg_neg_sig$Count[i],
        kegg_neg_sig$GeneRatio[i]))
  }
  write.csv(kegg_neg@result, file.path(results_dir, "KEGG_Kcnq1ot1_negative_corr.csv"),
            row.names = FALSE)
}

# ============================================================
# 4. KEGG enrichment - Positively correlated genes
# ============================================================
cat("\n=== KEGG Enrichment: Positively Correlated Genes ===\n")

kegg_pos <- enrichKEGG(
  gene = pos_entrez$ENTREZID,
  universe = all_entrez$ENTREZID,
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

if (!is.null(kegg_pos) && nrow(kegg_pos@result[kegg_pos@result$p.adjust < 0.05, ]) > 0) {
  kegg_pos_sig <- kegg_pos@result[kegg_pos@result$p.adjust < 0.05, ]
  cat(sprintf("\nSignificant KEGG pathways: %d\n", nrow(kegg_pos_sig)))
  for (i in 1:min(20, nrow(kegg_pos_sig))) {
    cat(sprintf("  %s (p.adj=%.2e, count=%d)\n",
        kegg_pos_sig$Description[i],
        kegg_pos_sig$p.adjust[i],
        kegg_pos_sig$Count[i]))
  }
  write.csv(kegg_pos@result, file.path(results_dir, "KEGG_Kcnq1ot1_positive_corr.csv"),
            row.names = FALSE)
} else {
  cat("  No significant KEGG pathways.\n")
}

# ============================================================
# 5. Figure: KEGG dot plot
# ============================================================
cat("\n=== Generating KEGG Figures ===\n")

if (!is.null(kegg_neg) && nrow(kegg_neg@result[kegg_neg@result$p.adjust < 0.05, ]) > 0) {
  p_kegg <- dotplot(kegg_neg, showCategory = 20,
    title = "KEGG Pathways: Genes Negatively Correlated with Kcnq1ot1") +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 9)
    )
  ggsave(file.path(fig_dir, "Fig5_KEGG_negative_corr.pdf"), p_kegg,
         width = 10, height = 8)
  ggsave(file.path(fig_dir, "Fig5_KEGG_negative_corr.png"), p_kegg,
         width = 10, height = 8, dpi = 300)
  cat("  Fig5 saved.\n")
}

# ============================================================
# 6. Detailed mitochondrial gene analysis
# ============================================================
cat("\n=== Mitochondrial Gene Detail ===\n")

# Key mitochondrial complex genes
mito_complexes <- list(
  "Complex_I" = c("Ndufa9", "Ndufs1", "Ndufs2", "Ndufs3", "Ndufs7", "Ndufs8",
                  "Ndufv1", "Ndufv2", "Ndufa2", "Ndufa5", "Ndufa10", "Ndufb8", "Ndufb10"),
  "Complex_II" = c("Sdha", "Sdhb", "Sdhc", "Sdhd"),
  "Complex_III" = c("Uqcrc1", "Uqcrc2", "Uqcrb", "Uqcrfs1", "Uqcrh"),
  "Complex_IV" = c("Cox4i1", "Cox5a", "Cox5b", "Cox6a1", "Cox6b1", "Cox6c", "Cox7a2", "Cox7c"),
  "Complex_V" = c("Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5f1", "Atp5g1", "Atp5h", "Atp5j", "Atp5o"),
  "TCA_cycle" = c("Cs", "Aco2", "Idh2", "Idh3a", "Ogdh", "Sucla2", "Suclg1", "Sdha", "Fh1", "Mdh2"),
  "FAO" = c("Acadm", "Acadl", "Acadvl", "Hadha", "Hadhb", "Cpt1a", "Cpt2", "Acox1")
)

cat(sprintf("%-12s %-12s %10s\n", "Complex", "Gene", "Correlation"))
cat(paste(rep("-", 38), collapse = ""), "\n")

mito_results <- data.frame()
for (cx in names(mito_complexes)) {
  for (g in mito_complexes[[cx]]) {
    idx <- which(cor_df$gene_symbol == g)
    if (length(idx) > 0) {
      r <- cor_df$correlation[idx[1]]
      cat(sprintf("%-12s %-12s %10.3f %s\n", cx, g, r,
          ifelse(abs(r) > 0.5, "***", ifelse(abs(r) > 0.3, "**", ""))))
      mito_results <- rbind(mito_results, data.frame(
        complex = cx, gene = g, correlation = r
      ))
    }
  }
}

write.csv(mito_results, file.path(results_dir, "Kcnq1ot1_mitochondrial_correlations.csv"),
          row.names = FALSE)

# Summary by complex
cat("\n=== Mean Correlation by Mitochondrial Complex ===\n")
agg <- aggregate(correlation ~ complex, data = mito_results, FUN = function(x) {
  c(mean = mean(x), n = length(x))
})
for (i in 1:nrow(agg)) {
  cat(sprintf("  %-12s: mean r = %.3f (n=%d genes)\n",
      agg$complex[i], agg$correlation[i, "mean"], agg$correlation[i, "n"]))
}

# ============================================================
# 7. Figure: Mitochondrial complex correlation barplot
# ============================================================

mito_results$complex <- factor(mito_results$complex,
  levels = c("Complex_I", "Complex_II", "Complex_III", "Complex_IV", "Complex_V",
             "TCA_cycle", "FAO"))

complex_means <- aggregate(correlation ~ complex, data = mito_results, mean)
complex_means$se <- aggregate(correlation ~ complex, data = mito_results,
  function(x) sd(x)/sqrt(length(x)))$correlation

p_mito <- ggplot(complex_means, aes(x = complex, y = correlation, fill = correlation)) +
  geom_col(alpha = 0.8) +
  geom_errorbar(aes(ymin = correlation - se, ymax = correlation + se),
                width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dotted", color = "red") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(
    title = "Kcnq1ot1 Correlation with Mitochondrial Complex Genes",
    subtitle = "GSE246221 Batch1 (52 samples across disease stages)",
    x = "Mitochondrial Complex / Pathway",
    y = "Mean Pearson Correlation",
    caption = "Dotted red line: r = -0.5 threshold"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

ggsave(file.path(fig_dir, "Fig6_Mito_complex_correlation.pdf"), p_mito,
       width = 8, height = 5)
ggsave(file.path(fig_dir, "Fig6_Mito_complex_correlation.png"), p_mito,
       width = 8, height = 5, dpi = 300)
cat("\n  Fig6 saved.\n")

# ============================================================
# 8. Figure: Scatter plot - Kcnq1ot1 vs top mito genes
# ============================================================

# Reload expression data for scatter plots
gse246 <- read.delim(
  gzfile(file.path(data_dir, "GSE246221_rawcounts_allsamples_mouse.txt.gz")),
  row.names = 1, check.names = FALSE
)
batch1_cols <- grep("^Batch1_", colnames(gse246), value = TRUE)
gse246_b1 <- gse246[, batch1_cols]

y246 <- DGEList(counts = round(gse246_b1))
y246 <- calcNormFactors(y246)
cpm246 <- cpm(y246, log = TRUE)

# Disease stage for coloring
stage <- rep(NA, length(batch1_cols))
stage[grepl("Control07w", batch1_cols)] <- "Healthy"
stage[grepl("STZSCD08w", batch1_cols)] <- "Acute_STZ"
stage[grepl("STZHFD14w", batch1_cols)] <- "MASLD"
stage[grepl("STZHFD20w", batch1_cols)] <- "MASH"
stage[grepl("STZHFD32w", batch1_cols)] <- "Fibrosis"
stage[grepl("STZHFD44w|STZHFD50w|STZHFD56w", batch1_cols) & !grepl("_T\\d+$", batch1_cols)] <- "HCC_NonTumor"
stage[grepl("_T\\d+$", batch1_cols)] <- "HCC_Tumor"

kcnq_expr <- cpm246["ENSMUSG00000101609_Kcnq1ot1", ]

# Top 4 most negatively correlated mitochondrial genes
top_mito <- c("Sdha", "Sdhb", "Ndufa9", "Atp5a1")
scatter_data <- data.frame()
for (g in top_mito) {
  g_row <- grep(paste0("_", g, "$"), rownames(cpm246), value = TRUE)
  if (length(g_row) > 0) {
    r_val <- cor(kcnq_expr, cpm246[g_row[1], ])
    scatter_data <- rbind(scatter_data, data.frame(
      Kcnq1ot1 = kcnq_expr,
      mito_gene = cpm246[g_row[1], ],
      gene = paste0(g, " (r=", round(r_val, 3), ")"),
      stage = stage,
      stringsAsFactors = FALSE
    ))
  }
}

stage_colors <- c(
  "Healthy" = "#2166AC", "Acute_STZ" = "#67A9CF", "MASLD" = "#D1E5F0",
  "MASH" = "#FDDBC7", "Fibrosis" = "#EF8A62",
  "HCC_NonTumor" = "#B2182B", "HCC_Tumor" = "#67001F"
)

p_scatter <- ggplot(scatter_data, aes(x = Kcnq1ot1, y = mito_gene, color = stage)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  scale_color_manual(values = stage_colors) +
  labs(
    title = "Kcnq1ot1 vs Key Mitochondrial Genes",
    subtitle = "Negative correlation across liver disease progression",
    x = "Kcnq1ot1 log2(CPM)",
    y = "Mitochondrial gene log2(CPM)",
    color = "Disease Stage"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "italic")
  )

ggsave(file.path(fig_dir, "Fig7_Kcnq1ot1_vs_mito_scatter.pdf"), p_scatter,
       width = 10, height = 8)
ggsave(file.path(fig_dir, "Fig7_Kcnq1ot1_vs_mito_scatter.png"), p_scatter,
       width = 10, height = 8, dpi = 300)
cat("  Fig7 saved.\n")

cat("\nStep 1 (KEGG + Mitochondrial analysis) complete!\n")

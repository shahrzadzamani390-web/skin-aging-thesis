# ==============================================================================
# 03_figures.R
# Thesis figures for genome-wide (mRNA) expression profiling with age
# Dataset: GSE226189 (Tsitsipatis et al. 2023)
#
# Figures produced (styled after Tsitsipatis et al. 2023, Figure 1):
#   1. PCA / QC plot (VST, blind to design)
#   2. Volcano plot -- linear model (beta coefficient vs. -log10 p-value)
#   3. Heatmap -- top 10 elevated / top 10 reduced mRNAs (linear model)
#   4. Age-trend curves for top mRNAs detected only by the spline model
#      (non-monotonic trajectories the linear model misses)
# ==============================================================================

library(tidyverse)
library(DESeq2)
library(splines)

dds_linear <- readRDS("data/dds_linear.rds")
dds_spline <- readRDS("data/dds_spline.rds")
res_linear_df <- readRDS("data/res_linear.rds")
res_spline_df <- readRDS("data/res_spline.rds")
metadata_ord <- as.data.frame(colData(dds_linear))

# ------------------------------------------------------------------------
# 1. PCA / QC plot
# ------------------------------------------------------------------------
dds_qc <- DESeqDataSetFromMatrix(
  countData = counts(dds_linear),
  colData = metadata_ord,
  design = ~ sex
)
vsd_qc <- vst(dds_qc, blind = TRUE)
pca_data <- plotPCA(vsd_qc, intgroup = c("age", "sex"), returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))

fig_pca <- ggplot(pca_data, aes(PC1, PC2, color = age, shape = sex)) +
  geom_point(size = 2.5) +
  scale_color_viridis_c() +
  labs(x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)"),
       color = "Age", shape = "Sex",
       title = "PCA of mRNA expression (VST, blind)")
fig_pca

# NOTE: PC1 (~49% variance) separates samples into two clusters not clearly
# explained by age or sex -- likely reflects the unrecorded "collection
# batch" covariate the original study adjusted for but which is not
# available in public GEO metadata. Worth noting as a limitation.

# ------------------------------------------------------------------------
# 2. Volcano plot -- linear model
# ------------------------------------------------------------------------
res_linear_df <- res_linear_df |>
  mutate(sig = case_when(
    is.na(pvalue) ~ "NS",
    log2FoldChange > 0 & pvalue < 0.05 ~ "up",
    log2FoldChange < 0 & pvalue < 0.05 ~ "down",
    TRUE ~ "NS"
  ))

fig_volcano_linear <- ggplot(res_linear_df, aes(log2FoldChange, -log10(pvalue), color = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c(up = "firebrick", down = "steelblue", NS = "grey80")) +
  labs(x = "Beta coefficient (per year of age)", y = "-log10(p-value)",
       color = NULL, title = "mRNA differential expression with age (linear model)")
fig_volcano_linear

# ------------------------------------------------------------------------
# 3. Heatmap -- top 10 elevated / top 10 reduced mRNAs (linear model)
# ------------------------------------------------------------------------
top_up <- res_linear_df |> filter(log2FoldChange > 0) |> slice_min(pvalue, n = 10)
top_down <- res_linear_df |> filter(log2FoldChange < 0) |> slice_min(pvalue, n = 10)
heatmap_genes <- bind_rows(top_up, top_down)

vsd_linear <- vst(dds_linear, blind = FALSE)
mat <- assay(vsd_linear)[heatmap_genes$ensembl_gene_id, ]
mat_z <- t(scale(t(mat)))
rownames(mat_z) <- heatmap_genes$hgnc_symbol

heatmap_long <- as.data.frame(mat_z) |>
  rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "sample_id", values_to = "z") |>
  left_join(metadata_ord |> rownames_to_column("sample_id") |> dplyr::select(sample_id, age),
            by = "sample_id") |>
  mutate(
    gene = factor(gene, levels = rev(heatmap_genes$hgnc_symbol)),
    sample_id = factor(sample_id, levels = rownames(metadata_ord))
  )

fig_heatmap <- ggplot(heatmap_long, aes(sample_id, gene, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        limits = c(-2, 2), oob = scales::squish) +
  labs(x = "Samples (ordered by age)", y = NULL, fill = "Normalized\nRNA abundance",
       title = "Top 10 elevated / top 10 reduced mRNAs with age (linear model)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
fig_heatmap

# ------------------------------------------------------------------------
# 4. Age-trend curves for top mRNAs detected only by the spline model
# ------------------------------------------------------------------------
spline_specific <- res_spline_df |>
  filter(pvalue < 0.05) |>
  left_join(
    res_linear_df |> dplyr::select(ensembl_gene_id, pvalue_linear = pvalue),
    by = "ensembl_gene_id"
  ) |>
  filter(pvalue_linear > 0.05, !is.na(hgnc_symbol), hgnc_symbol != "") |>
  arrange(pvalue)

top_spline_genes <- spline_specific |> slice_min(pvalue, n = 6)

norm_counts <- counts(dds_spline, normalized = TRUE)
log_counts <- log2(norm_counts + 1)

spline_curve_data <- log_counts[top_spline_genes$ensembl_gene_id, ] |>
  as.data.frame() |>
  rownames_to_column("ensembl_gene_id") |>
  pivot_longer(-ensembl_gene_id, names_to = "sample_id", values_to = "log2_abundance") |>
  left_join(top_spline_genes |> dplyr::select(ensembl_gene_id, hgnc_symbol), by = "ensembl_gene_id") |>
  left_join(metadata_ord |> rownames_to_column("sample_id") |> dplyr::select(sample_id, age),
            by = "sample_id") |>
  mutate(hgnc_symbol = factor(hgnc_symbol, levels = top_spline_genes$hgnc_symbol))

fig_spline_curves <- ggplot(spline_curve_data, aes(age, log2_abundance)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ ns(x, df = 3), color = "firebrick", se = TRUE) +
  facet_wrap(~ hgnc_symbol, scales = "free_y", nrow = 2) +
  labs(x = "Age (years)", y = "mRNA abundance (Log2)",
       title = "Top mRNAs detected only by the spline model")
fig_spline_curves

# ------------------------------------------------------------------------
# Save figures
# ------------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE)
ggsave("figures/01_pca_qc.png", fig_pca, width = 6, height = 5, dpi = 300)
ggsave("figures/02_volcano_linear.png", fig_volcano_linear, width = 6, height = 5, dpi = 300)
ggsave("figures/03_heatmap_top_genes.png", fig_heatmap, width = 8, height = 5, dpi = 300)
ggsave("figures/04_spline_curves.png", fig_spline_curves, width = 8, height = 5, dpi = 300)

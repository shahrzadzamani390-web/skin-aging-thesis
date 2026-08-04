# ==============================================================================
# 02_qc_and_deseq2.R
# Genome-wide (mRNA) expression profiling with age in skin fibroblasts
# Dataset: GSE226189 (Tsitsipatis et al. 2023), n=82, ages 22-89
#
# Mirrors the original study's DESeq2 approach:
#   - Linear model: age (continuous, per year) + sex, Wald test
#   - Spline model: natural spline (df=3) on age + sex, LRT vs sex-only model
#   - Low-count filter: 82 samples split into 4 consecutive age-ordered
#     groups (~20/group); keep genes with >=10 counts in >=7 samples of
#     at least one group
#
# NOTE: the original study also adjusted for "collection batch," which is
# not available in the public GEO metadata. Models here adjust for sex only.
# ==============================================================================

library(tidyverse)
library(DESeq2)
library(splines)

count_matrix_mrna <- readRDS("data/count_matrix_mrna.rds")
metadata_clean <- readRDS("data/metadata_clean.rds")
gene_symbol_map <- readRDS("data/gene_symbol_map.rds")

# --- Order samples by age for the low-count filtering scheme ---
ord <- order(metadata_clean$age)
metadata_ord <- metadata_clean[ord, ]
metadata_ord$sex <- factor(metadata_ord$sex)
count_mrna_ord <- count_matrix_mrna[, ord]

# --- Low-count filter (paper's scheme) ---
n <- nrow(metadata_ord)
grp_size <- ceiling(n / 4)
filter_group <- ceiling(seq_len(n) / grp_size)
table(filter_group)  # ~20 samples/group, age-consecutive

keep <- sapply(seq_len(nrow(count_mrna_ord)), function(i) {
  x <- count_mrna_ord[i, ]
  any(tapply(x, filter_group, function(v) sum(v >= 10) >= 7))
})
table(keep)

count_mrna_filtered <- count_mrna_ord[keep, ]
dim(count_mrna_filtered)  # ~15910 x 82

# ------------------------------------------------------------------------
# QC: library sizes and PCA (on VST-transformed data)
# ------------------------------------------------------------------------
dds_qc <- DESeqDataSetFromMatrix(
  countData = count_mrna_filtered,
  colData = metadata_ord,
  design = ~ sex
)

lib_sizes <- colSums(count_mrna_filtered)
metadata_ord$lib_size <- lib_sizes

ggplot(metadata_ord, aes(age, lib_size)) +
  geom_point() +
  labs(x = "Age (years)", y = "Library size (total counts)")

vsd <- vst(dds_qc, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = c("age", "sex"), returnData = TRUE)

ggplot(pca_data, aes(PC1, PC2, color = age)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  labs(x = "PC1", y = "PC2", color = "Age")

# ------------------------------------------------------------------------
# Linear model: age (continuous) + sex, Wald test
# ------------------------------------------------------------------------
dds_linear <- DESeqDataSetFromMatrix(
  countData = count_mrna_filtered,
  colData = metadata_ord,
  design = ~ sex + age
)
dds_linear <- DESeq(dds_linear)
res_linear <- results(dds_linear, name = "age") |>
  as.data.frame() |>
  rownames_to_column("ensembl_gene_id") |>
  left_join(gene_symbol_map, by = "ensembl_gene_id") |>
  arrange(pvalue)

sum(res_linear$pvalue < 0.05, na.rm = TRUE)  # unadjusted p<0.05, cf. 1437 in paper

# ------------------------------------------------------------------------
# Spline model: natural spline (df=3) on age + sex, LRT
# ------------------------------------------------------------------------
dds_spline <- DESeqDataSetFromMatrix(
  countData = count_mrna_filtered,
  colData = metadata_ord,
  design = ~ sex + ns(age, df = 3)
)
dds_spline <- DESeq(dds_spline, test = "LRT", reduced = ~ sex)
res_spline <- results(dds_spline) |>
  as.data.frame() |>
  rownames_to_column("ensembl_gene_id") |>
  left_join(gene_symbol_map, by = "ensembl_gene_id") |>
  arrange(pvalue)

sum(res_spline$pvalue < 0.05, na.rm = TRUE)  # unadjusted p<0.05, cf. 1436 in paper

# --- Save results for downstream (volcano plots, NR3C2 deep-dive, etc.) ---
saveRDS(dds_linear, "data/dds_linear.rds")
saveRDS(dds_spline, "data/dds_spline.rds")
saveRDS(res_linear, "data/res_linear.rds")
saveRDS(res_spline, "data/res_spline.rds")

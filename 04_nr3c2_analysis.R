# ==============================================================================
# 04_nr3c2_analysis.R
# Focused analysis: NR3C2 (Mineralocorticoid Receptor) mRNA abundance vs. age
# Dataset: GSE226189 (Tsitsipatis et al. 2023)
# ==============================================================================

library(tidyverse)
library(DESeq2)
library(splines)

dds_linear <- readRDS("data/dds_linear.rds")
dds_spline <- readRDS("data/dds_spline.rds")
res_linear_df <- readRDS("data/res_linear.rds")
res_spline_df <- readRDS("data/res_spline.rds")
metadata_ord <- as.data.frame(colData(dds_linear))

nr3c2_id <- "ENSG00000151623"
stopifnot(nr3c2_id %in% rownames(dds_linear))

# ------------------------------------------------------------------------
# Significance in both models
# ------------------------------------------------------------------------
nr3c2_linear <- res_linear_df |> filter(ensembl_gene_id == nr3c2_id)
nr3c2_spline <- res_spline_df |> filter(ensembl_gene_id == nr3c2_id)

nr3c2_linear
nr3c2_spline

# Linear model: beta = 0.0109/year, unadjusted p = 0.0096, FDR padj = 0.537
# Spline model: LRT unadjusted p = 0.0013, FDR padj = 0.233
# -> NR3C2 shows a nominally significant, modestly increasing trend with age
#    in the linear model, and a *stronger* signal under the spline model,
#    suggesting a non-linear (accelerating) trajectory -- consistent with
#    the dip-then-rise pattern visible in the age-trend plot below.
#    Neither model survives genome-wide FDR correction (padj > 0.05), so
#    this should be reported as a candidate/nominal finding, not a
#    genome-wide-significant one.

# ------------------------------------------------------------------------
# Age-trend figure
# ------------------------------------------------------------------------
norm_counts_linear <- counts(dds_linear, normalized = TRUE)
nr3c2_log2 <- log2(norm_counts_linear[nr3c2_id, ] + 1)

nr3c2_df <- metadata_ord |>
  mutate(log2_abundance = nr3c2_log2[sample_id])

cor.test(nr3c2_df$age, nr3c2_df$log2_abundance)

fig_nr3c2 <- ggplot(nr3c2_df, aes(age, log2_abundance)) +
  geom_point(aes(color = age_group), size = 2, alpha = 0.8) +
  geom_smooth(aes(linetype = "Linear fit"), method = "lm", se = TRUE, color = "black") +
  geom_smooth(aes(linetype = "Spline fit (df=3)"), method = "lm",
              formula = y ~ ns(x, df = 3), se = TRUE, color = "firebrick") +
  scale_color_manual(values = c(young = "#2166ac", middle = "grey60", old = "#b2182b")) +
  labs(x = "Age (years)", y = "NR3C2 mRNA abundance (Log2)",
       color = "Age group", linetype = NULL,
       title = "NR3C2 mRNA abundance across the human lifespan") +
  theme_minimal()
fig_nr3c2

dir.create("figures", showWarnings = FALSE)
ggsave("figures/05_nr3c2_age_trend.png", fig_nr3c2, width = 7, height = 5, dpi = 300)

# ------------------------------------------------------------------------
# Does the NR3C2 age trend differ by sex?
# ------------------------------------------------------------------------
# Genome-wide interaction model (~ sex + age + sex:age), sharing dispersion
# estimates across all filtered genes, then extract the interaction term
# for NR3C2 specifically.
dds_interaction <- DESeqDataSetFromMatrix(
  countData = counts(dds_linear),
  colData = metadata_ord,
  design = ~ sex + age + sex:age
)
dds_interaction <- DESeq(dds_interaction)

nr3c2_interaction <- results(dds_interaction, name = "sexMale.age") |>
  as.data.frame() |>
  rownames_to_column("ensembl_gene_id") |>
  filter(ensembl_gene_id == nr3c2_id)
nr3c2_interaction
# beta = -0.0072/year, p = 0.42, padj = 1.0 -- no evidence the age slope
# differs by sex.

# Descriptive sex-stratified slopes/correlations (for context only; the
# interaction test above is the actual significance test)
nr3c2_df |>
  group_by(sex) |>
  summarise(
    n = n(),
    cor_age = cor(age, log2_abundance),
    slope = coef(lm(log2_abundance ~ age))[["age"]]
  )
# Female: r = 0.46, slope = 0.015/year; Male: r = 0.25, slope = 0.010/year
# -- descriptively steeper/tighter in females, but not statistically
# distinguishable given the interaction test above.

fig_nr3c2_sex <- ggplot(nr3c2_df, aes(age, log2_abundance, color = sex)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = c(Female = "#d6604d", Male = "#4393c3")) +
  labs(x = "Age (years)", y = "NR3C2 mRNA abundance (Log2)", color = "Sex",
       title = "NR3C2 mRNA abundance vs. age, stratified by sex",
       subtitle = "Sex \u00d7 age interaction: p = 0.42 (not significant)") +
  theme_minimal()
fig_nr3c2_sex

ggsave("figures/06_nr3c2_by_sex.png", fig_nr3c2_sex, width = 7, height = 5, dpi = 300)

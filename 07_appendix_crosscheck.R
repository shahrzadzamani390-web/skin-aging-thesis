# 07_appendix_crosscheck.R
# Cross-check NR3C2 against Tsitsipatis et al. 2023 Appendix S1
# (original paper's own supplementary differential expression gene lists)
# =============================================================================

library(readxl)
library(purrr)
library(dplyr)

# File saved locally in the project's data/ folder
# (downloaded manually from https://pmc.ncbi.nlm.nih.gov/articles/PMC10652340/
#  due to bot-protection blocking scripted downloads from the Wiley site)
appendix_path <- "data/ACEL-22-e13915-s004 (1).xlsx"

excel_sheets(appendix_path)
# [1] "coding transcripts-linear model" "coding-linear-females" "coding-linear-males"
# [4] "coding transcripts-spline model" "coding-spline-females" "coding-spline-males"
# [7] "Demographic data"

data_sheets <- c("coding transcripts-linear model", "coding-linear-females", "coding-linear-males",
                  "coding transcripts-spline model", "coding-spline-females", "coding-spline-males")

# read as text first (one sheet had a stray non-numeric value that broke type-guessing
# when binding rows across sheets), then convert the numeric columns explicitly
nr3c2_appendix <- map_dfr(data_sheets, function(s) {
  df <- read_excel(appendix_path, sheet = s, col_types = "text")
  df <- df[df$Gene_Name == "NR3C2", ]
  df$model_sheet <- s
  df
})

nr3c2_appendix$baseMean <- as.numeric(nr3c2_appendix$baseMean)
nr3c2_appendix$log2FoldChange <- as.numeric(nr3c2_appendix$log2FoldChange)
nr3c2_appendix$pvalue <- as.numeric(nr3c2_appendix$pvalue)
nr3c2_appendix$padj <- as.numeric(nr3c2_appendix$padj)

nr3c2_appendix_summary <- nr3c2_appendix[, c("model_sheet", "Gene_ID", "baseMean",
                                               "log2FoldChange", "pvalue", "padj")]
nr3c2_appendix_summary

write.csv(nr3c2_appendix_summary, "data/nr3c2_appendix_s1_crosscheck.csv", row.names = FALSE)

# RESULTS (all 6 model sheets):
#   Pooled linear:   p=0.0091, padj=0.470
#   Female linear:   p=0.0026, padj=0.405, log2FC=0.0164 (increasing)
#   Male linear:     p=0.136,  padj=1.000 (not significant)
#   Pooled spline:   p=0.0016, padj=0.292
#   Female spline:   p=0.0043, padj=0.064 (closest to FDR significance across entire analysis)
#   Male spline:     p=0.052,  padj=1.000
#
# Interpretation: strong independent validation of the reanalysis pipeline - pooled linear/spline
# p-values closely match this project's own results (p=0.0096 / p=0.0013). Female-specific results
# consistently stronger than male across both models, independently confirming this project's own
# sex-stratified finding (female r=0.46 vs. male r=0.25 age correlation).

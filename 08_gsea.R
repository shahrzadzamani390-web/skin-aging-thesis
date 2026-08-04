# 08_gsea.R
# Gene Set Enrichment Analysis (fgsea) on the four themes named in the proposal:
# ECM organization, cellular senescence, oxidative stress, autophagy
# References: Subramanian et al. 2005 PNAS (GSEA method); Korotkevich et al. 2021 (fgsea);
#             Liberzon et al. 2015 Cell Systems (MSigDB); Mubeen et al. 2022 (methods justification)
# Requires: res_linear_df (from 02_qc_and_deseq2.R)
# =============================================================================

library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Pull gene sets for the four proposal themes
# -----------------------------------------------------------------------------
# Note: no exact GO:BP match exists for generic "ECM organization" or "autophagy" -
# used REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION (direct match) and
# GOBP_MACROAUTOPHAGY (macroautophagy = the canonical/default form of autophagy).
go_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
reactome_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

ecm_genes_gsea        <- reactome_sets |> filter(gs_name == "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION") |> pull(gene_symbol)
senescence_genes_gsea <- go_sets |> filter(gs_name == "GOBP_CELLULAR_SENESCENCE") |> pull(gene_symbol)
oxidative_genes_gsea  <- go_sets |> filter(gs_name == "GOBP_RESPONSE_TO_OXIDATIVE_STRESS") |> pull(gene_symbol)
autophagy_genes_gsea  <- go_sets |> filter(gs_name == "GOBP_MACROAUTOPHAGY") |> pull(gene_symbol)

pathway_list <- list(
  ECM_ORGANIZATION = unique(ecm_genes_gsea),
  CELLULAR_SENESCENCE = unique(senescence_genes_gsea),
  OXIDATIVE_STRESS_RESPONSE = unique(oxidative_genes_gsea),
  MACROAUTOPHAGY = unique(autophagy_genes_gsea)
)
sapply(pathway_list, length)
# ECM_ORGANIZATION=321 (285 after filtering to expressed genes), CELLULAR_SENESCENCE=109 (92),
# OXIDATIVE_STRESS_RESPONSE=437 (394), MACROAUTOPHAGY=392 (375)

# -----------------------------------------------------------------------------
# 2. Build ranked gene list (Wald stat from linear age model)
# -----------------------------------------------------------------------------
ranked_genes <- res_linear_df |>
  filter(!is.na(hgnc_symbol), hgnc_symbol != "", !is.na(stat)) |>
  group_by(hgnc_symbol) |>
  summarise(stat = mean(stat)) |>
  arrange(desc(stat))

gene_ranks <- ranked_genes$stat
names(gene_ranks) <- ranked_genes$hgnc_symbol
# n=15,638 uniquely-named genes ranked

# -----------------------------------------------------------------------------
# 3. Run fgsea
# -----------------------------------------------------------------------------
set.seed(42)
gsea_results <- fgsea(pathways = pathway_list, stats = gene_ranks, minSize = 10, maxSize = 500)
gsea_results[order(pval), .(pathway, pval, padj, NES, size)]

# RESULTS:
#   MACROAUTOPHAGY:             pval=0.0001, padj=0.0004, NES=+1.57, n=375  ** SIGNIFICANT **
#   ECM_ORGANIZATION:           pval=0.088,  padj=0.176,  NES=-1.18, n=285
#   OXIDATIVE_STRESS_RESPONSE:  pval=0.192,  padj=0.256,  NES=+1.11, n=394
#   CELLULAR_SENESCENCE:        pval=0.789,  padj=0.789,  NES=-0.85, n=92
#
# Interpretation: macroautophagy is the only theme reaching FDR significance - genes in this
# set skew toward increasing expression with age (first FDR-significant pathway-level result
# in the entire thesis). The other three themes (ECM, oxidative stress, senescence) show no
# significant enrichment, consistent with the largely null/weak candidate-gene relationships
# found for these pathways throughout the rest of the analysis.

# -----------------------------------------------------------------------------
# 4. Enrichment plot for the significant pathway
# -----------------------------------------------------------------------------
fig_gsea_autophagy <- plotEnrichment(pathway_list[["MACROAUTOPHAGY"]], gene_ranks) +
  labs(title = "GSEA: Macroautophagy (GOBP_MACROAUTOPHAGY)",
       subtitle = "NES = 1.57, padj = 0.0004, n = 375 genes")

ggsave("figures/24_gsea_autophagy.png", fig_gsea_autophagy, width = 7, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 5. Combined 4-panel figure (all proposal themes, not just the significant one)
# -----------------------------------------------------------------------------
library(patchwork)

gsea_results_ordered <- gsea_results |>
  arrange(match(pathway, c("MACROAUTOPHAGY", "ECM_ORGANIZATION",
                            "OXIDATIVE_STRESS_RESPONSE", "CELLULAR_SENESCENCE")))

plot_list <- lapply(gsea_results_ordered$pathway, function(p) {
  row <- gsea_results_ordered[gsea_results_ordered$pathway == p, ]
  plotEnrichment(pathway_list[[p]], gene_ranks) +
    labs(title = p,
         subtitle = paste0("NES=", round(row$NES, 2), ", padj=", signif(row$padj, 2))) +
    theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 9))
})

fig_gsea_all <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(title = "GSEA across all four proposal themes")

ggsave("figures/24_gsea_all_four.png", fig_gsea_all, width = 11, height = 8, dpi = 300)

# 05_candidate_panel.R
# MR/GR signaling candidate gene panel: correlation network, category modules,
# IGF1/LOX/SGK1 checks, and HSD11B2 selectivity analysis
# Requires: dds_linear, gene_symbol_map, res_linear_df, res_spline_df (from 02_qc_and_deseq2.R)
# =============================================================================

library(tidyverse)
library(pheatmap)
library(igraph)
library(ggraph)

norm_counts_linear <- counts(dds_linear, normalized = TRUE)

# -----------------------------------------------------------------------------
# 1. Build the 34-gene candidate panel
# -----------------------------------------------------------------------------
gene_panel <- c("NR3C2","NR3C1","HSD11B2","HSD11B1","H6PD","SGK1","NEDD4L",
  "SCNN1A","SCNN1B","SCNN1G","ATP1A1","ATP1B1","SRC","EGFR","MAPK1","MAPK3",
  "MAPK8","MAPK9","AKT1","RAC1","NOX1","CYBB","NFKB1","RELA","TNF","IL6",
  "PTGS2","COL1A1","COL3A1","MMP1","MMP3","TIMP1","CDKN1A","CDKN2A","TP53")

panel_map <- gene_symbol_map |> filter(hgnc_symbol %in% gene_panel)
panel_map$in_filtered <- panel_map$ensembl_gene_id %in% rownames(count_mrna_filtered)
panel_final <- panel_map |> filter(in_filtered)  # SCNN1G dropped here (filtered out, low counts)

panel_expr <- norm_counts_linear[panel_final$ensembl_gene_id, ]
rownames(panel_expr) <- panel_final$hgnc_symbol
panel_expr_log2 <- log2(panel_expr + 1)

# -----------------------------------------------------------------------------
# 2. Gene-gene correlation matrix
# -----------------------------------------------------------------------------
panel_cor <- cor(t(panel_expr_log2), method = "pearson")

# -----------------------------------------------------------------------------
# 3. Correlation network graph (|r| > 0.5)
# -----------------------------------------------------------------------------
gene_category <- data.frame(
  hgnc_symbol = c("NR3C2","NR3C1","HSD11B2","HSD11B1","H6PD","SGK1",
    "NEDD4L","SCNN1A","SCNN1B","ATP1A1","ATP1B1",
    "SRC","EGFR","MAPK1","MAPK3","MAPK8","MAPK9","AKT1","RAC1",
    "NOX1","CYBB",
    "NFKB1","RELA","TNF","IL6","PTGS2",
    "COL1A1","COL3A1","MMP1","MMP3","TIMP1",
    "CDKN1A","CDKN2A","TP53"),
  category = c(rep("Receptor/Enzyme", 6),
    rep("Transport", 5),
    rep("Rapid signaling", 8),
    rep("Oxidative stress", 2),
    rep("Inflammation", 5),
    rep("ECM", 5),
    rep("Senescence/Stress", 3))
)
rownames(gene_category) <- gene_category$hgnc_symbol

cor_long <- panel_cor |>
  as.data.frame() |>
  rownames_to_column("gene1") |>
  pivot_longer(-gene1, names_to = "gene2", values_to = "r") |>
  filter(gene1 != gene2) |>
  rowwise() |>
  mutate(pair = paste(sort(c(gene1, gene2)), collapse = "_")) |>
  ungroup() |>
  distinct(pair, .keep_all = TRUE) |>
  dplyr::select(gene1, gene2, r)

edges <- cor_long |> filter(abs(r) > 0.5)
edges_for_graph <- edges
colnames(edges_for_graph) <- c("from", "to", "weight")

nodes <- data.frame(name = rownames(panel_cor)) |>
  left_join(gene_category, by = c("name" = "hgnc_symbol"))

g <- graph_from_data_frame(d = edges_for_graph, vertices = nodes, directed = FALSE)

set.seed(42)
layout_weights <- abs(E(g)$weight)

fig_network <- ggraph(g, layout = "fr", weights = layout_weights) +
  geom_edge_link(aes(edge_alpha = abs(weight), edge_color = weight), edge_width = 0.6) +
  scale_edge_color_gradient2(low = "steelblue", mid = "grey90", high = "firebrick", midpoint = 0) +
  scale_edge_alpha(range = c(0.15, 0.8)) +
  geom_node_point(aes(color = category), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.2, fontface = "bold") +
  labs(title = "MR/GR signaling panel: correlation network (|r| > 0.5)",
       edge_color = "r", edge_alpha = NULL) +
  theme_void()

ggsave("figures/16_candidate_network.png", fig_network, width = 8, height = 7, dpi = 300)

# Finding: NR3C2 is isolated in this network (no |r|>0.5 edges); NR3C1 is well-connected.

# -----------------------------------------------------------------------------
# 4. Per-category correlation heatmaps (NR3C2 + NR3C1 included in each)
# -----------------------------------------------------------------------------
category_heatmap <- function(cat_genes, cat_name, filename) {
  genes_to_plot <- unique(c("NR3C2", "NR3C1", cat_genes))
  sub_cor <- panel_cor[genes_to_plot, genes_to_plot]
  png(paste0("figures/", filename), width = 6, height = 5.5, units = "in", res = 300)
  pheatmap(sub_cor,
    color = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
    breaks = seq(-1, 1, length.out = 101),
    main = paste0(cat_name, " module (with NR3C2/NR3C1)"),
    display_numbers = TRUE, number_format = "%.2f",
    fontsize_number = 9, fontsize = 10)
  dev.off()
}

ecm_genes           <- c("COL1A1", "COL3A1", "MMP1", "MMP3", "TIMP1")
inflammation_genes   <- c("NFKB1", "RELA", "TNF", "IL6", "PTGS2")
oxidative_genes      <- c("NOX1", "CYBB")
transport_genes      <- c("NEDD4L", "SCNN1A", "SCNN1B", "ATP1A1", "ATP1B1")
rapid_signaling_genes <- c("SRC", "EGFR", "MAPK1", "MAPK3", "MAPK8", "MAPK9", "AKT1", "RAC1")
senescence_genes     <- c("CDKN1A", "CDKN2A", "TP53")

category_heatmap(ecm_genes, "ECM", "17_module_ecm.png")
category_heatmap(inflammation_genes, "Inflammation", "18_module_inflammation.png")
category_heatmap(oxidative_genes, "Oxidative stress", "19_module_oxidative.png")
category_heatmap(transport_genes, "Transport", "20_module_transport.png")
category_heatmap(rapid_signaling_genes, "Rapid signaling", "21_module_rapid_signaling.png")
category_heatmap(senescence_genes, "Senescence/Stress", "22_module_senescence.png")

# Key numbers (see thesis text for full interpretation):
# ECM: NR3C2 all |r|<=0.16; NR3C1 strongest w/ TIMP1 r=-0.45; COL1A1-COL3A1 r=0.55; MMP1-MMP3 r=0.61
# Inflammation: NR3C1 negative w/ RELA r=-0.62, NFKB1 r=-0.45, TNF r=-0.44 (GR represses NFkB)
# Oxidative stress: weak throughout, NOX1-CYBB r=0.04 (distinct isoforms, expected)
# Transport: NR3C1 mixed signs w/ ENaC/Na-K-ATPase; ATP1A1-ATP1B1 r=-0.65 (unexpected anomaly)
# Rapid signaling: two anti-correlated clusters (SRC/AKT1/MAPK3 vs MAPK1/8/9);
#   NR3C2 negative w/ SRC/AKT1 cluster (r=-0.32 to -0.38), positive w/ MAPK1/8/9 (r=0.29-0.42)
# Senescence: NR3C1 vs CDKN1A r=-0.77 (strongest single correlation in whole panel)

# -----------------------------------------------------------------------------
# 5. IGF1 and LOX individual checks (SGK1 already covered above)
# -----------------------------------------------------------------------------
igf1_id <- gene_symbol_map |> filter(hgnc_symbol == "IGF1") |> pull(ensembl_gene_id)
lox_id  <- gene_symbol_map |> filter(hgnc_symbol == "LOX")  |> pull(ensembl_gene_id)
sgk1_id <- gene_symbol_map |> filter(hgnc_symbol == "SGK1") |> pull(ensembl_gene_id)
nr3c2_id <- gene_symbol_map |> filter(hgnc_symbol == "NR3C2") |> pull(ensembl_gene_id)
nr3c1_id <- gene_symbol_map |> filter(hgnc_symbol == "NR3C1") |> pull(ensembl_gene_id)

res_linear_df |> filter(hgnc_symbol %in% c("IGF1", "LOX", "SGK1"))
res_spline_df |> filter(hgnc_symbol %in% c("IGF1", "LOX", "SGK1"))

igf1_log2 <- log2(norm_counts_linear[igf1_id, ] + 1)
lox_log2  <- log2(norm_counts_linear[lox_id, ] + 1)
nr3c2_log2 <- log2(norm_counts_linear[nr3c2_id, ] + 1)
nr3c1_log2 <- log2(norm_counts_linear[nr3c1_id, ] + 1)

candidate_summary <- data.frame(
  gene = rep(c("IGF1", "LOX", "SGK1"), each = 2),
  receptor = rep(c("NR3C2", "NR3C1"), times = 3),
  r = c(cor(igf1_log2, nr3c2_log2), cor(igf1_log2, nr3c1_log2),
        cor(lox_log2, nr3c2_log2),  cor(lox_log2, nr3c1_log2),
        panel_cor["SGK1", "NR3C2"], panel_cor["SGK1", "NR3C1"])
)

fig_candidate_summary <- ggplot(candidate_summary, aes(x = r, y = gene, color = receptor)) +
  geom_segment(aes(x = 0, xend = r, y = gene, yend = gene),
               position = position_dodge(width = 0.5), linewidth = 0.8) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c(NR3C2 = "#2166ac", NR3C1 = "#b2182b")) +
  labs(x = "Correlation with receptor (r)", y = NULL, color = "Receptor",
       title = "IGF1, LOX and SGK1: correlation with NR3C2 and NR3C1") +
  theme_minimal()

ggsave("figures/10_candidate_gene_summary.png", fig_candidate_summary, width = 7, height = 5, dpi = 300)

# Key numbers: IGF1 nominal age trend (linear p=0.036, spline p=0.004), LOX/SGK1 not significant.
# LOX strongest receptor correlation of the three: r=0.59 with NR3C1.

# -----------------------------------------------------------------------------
# 6. HSD11B2 / MR selectivity check
# -----------------------------------------------------------------------------
hsd11b2_id <- gene_symbol_map |> filter(hgnc_symbol == "HSD11B2") |> pull(ensembl_gene_id)

res_linear_df |> filter(hgnc_symbol == "HSD11B2")
round(panel_cor["HSD11B2", c("NR3C2", "NR3C1", "SGK1")], 3)

expr_summary <- data.frame(
  gene = c("HSD11B2", "NR3C2", "SGK1", "NR3C1"),
  mean_expr = c(
    mean(norm_counts_linear[hsd11b2_id, ]),
    mean(norm_counts_linear[nr3c2_id, ]),
    mean(norm_counts_linear[sgk1_id, ]),
    mean(norm_counts_linear[nr3c1_id, ])
  )
)
expr_summary$gene <- factor(expr_summary$gene, levels = expr_summary$gene)

fig_selectivity <- ggplot(expr_summary, aes(gene, mean_expr, fill = gene)) +
  geom_col() +
  scale_y_log10() +
  scale_fill_manual(values = c(HSD11B2 = "grey60", NR3C2 = "#2166ac", SGK1 = "#4daf4a", NR3C1 = "#b2182b")) +
  labs(x = NULL, y = "Mean normalized expression (log10 scale)",
       title = "HSD11B2 expression relative to MR/GR pathway genes") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("figures/11_hsd11b2_selectivity.png", fig_selectivity, width = 6, height = 5, dpi = 300)

# Finding: HSD11B2 mean expr = 26.8, ~17x lower than NR3C2 (454), ~2000x lower than NR3C1 (52,371)
# -> classical aldosterone-selectivity mechanism may be largely absent in dermal fibroblasts
# HSD11B2 vs SGK1 r=-0.39 (consistent with reduced cortisol inactivation -> more receptor activation)

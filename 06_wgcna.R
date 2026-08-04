# 06_wgcna.R
# Weighted Gene Co-expression Network Analysis (WGCNA)
# Requires: dds_linear, vsd_linear, metadata_ord, nr3c2_log2, nr3c1_log2 (from earlier scripts)
# References: Langfelder & Horvath 2008 (BMC Bioinformatics); Zhang & Horvath 2005 (Stat Appl Genet Mol Biol)
# =============================================================================

library(WGCNA)
options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 1. Prepare input matrix and QC
# -----------------------------------------------------------------------------
wgcna_input <- t(assay(vsd_linear))  # samples x genes: 82 x 15910

gsg <- goodSamplesGenes(wgcna_input, verbose = 3)
stopifnot(gsg$allOK)  # TRUE - no genes/samples flagged

# -----------------------------------------------------------------------------
# 2. Sample clustering / outlier check
# -----------------------------------------------------------------------------
sample_tree <- hclust(dist(wgcna_input), method = "average")

png("figures/07_wgcna_sample_dendrogram.png", width = 10, height = 6, units = "in", res = 300)
plot(sample_tree, main = "Sample clustering to detect outliers",
     xlab = "", sub = "", cex = 0.6, labels = FALSE)
dev.off()

# Investigated a 67/15 split at cutHeight=150 - confirmed NOT a discrete outlier group
# (PC1 ranges overlap heavily between the two clusters: continuous distribution, not batch).
# Decision: retain all 82 samples, no removal.
clust <- cutreeStatic(sample_tree, cutHeight = 150, minSize = 10)
table(clust)

# -----------------------------------------------------------------------------
# 3. Soft-threshold power selection
# -----------------------------------------------------------------------------
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(wgcna_input, powerVector = powers, verbose = 5)
sft_df <- sft$fitIndices

fig_sft_r2 <- ggplot(sft_df, aes(Power, -sign(slope) * SFT.R.sq)) +
  geom_point(color = "firebrick", size = 2) +
  ggrepel::geom_text_repel(aes(label = Power), size = 3.5) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40") +
  labs(x = "Soft threshold power", y = expression("Scale-free topology fit, signed " * R^2),
       title = "Scale independence") +
  theme_minimal()

fig_mean_k <- ggplot(sft_df, aes(Power, mean.k.)) +
  geom_point(color = "firebrick", size = 2) +
  ggrepel::geom_text_repel(aes(label = Power), size = 3.5) +
  labs(x = "Soft threshold power", y = "Mean connectivity",
       title = "Mean connectivity") +
  theme_minimal()

ggsave("figures/08_wgcna_sft_r2.png", fig_sft_r2, width = 7, height = 5, dpi = 300)
ggsave("figures/09_wgcna_mean_connectivity.png", fig_mean_k, width = 7, height = 5, dpi = 300)

# Power = 6 selected: R^2 ~= 0.8, slope ~= -0.92 (close to -1), mean connectivity = 514
# (WGCNA-recommended default power for unsigned networks)
soft_power <- 6

# -----------------------------------------------------------------------------
# 4. Network construction: adjacency, TOM, gene dendrogram
# -----------------------------------------------------------------------------
adjacency_matrix <- adjacency(wgcna_input, power = soft_power)

tom_matrix <- TOMsimilarity(adjacency_matrix)
dissTOM <- 1 - tom_matrix

gene_tree <- hclust(as.dist(dissTOM), method = "average")

min_module_size <- 30
dynamic_mods <- cutreeDynamic(dendro = gene_tree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = min_module_size)
table(dynamic_mods)  # 20 modules + grey (unassigned, n=53)

dynamic_colors <- labels2colors(dynamic_mods)
table(dynamic_colors)

png("figures/23_wgcna_gene_dendrogram.png", width = 10, height = 6, units = "in", res = 300)
plotDendroAndColors(gene_tree, dynamic_colors, "Module colors",
                     dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05,
                     main = "Gene dendrogram and module colors")
dev.off()

# -----------------------------------------------------------------------------
# 5. Module eigengenes and module-trait relationships
# -----------------------------------------------------------------------------
module_eigengenes <- moduleEigengenes(wgcna_input, colors = dynamic_colors)$eigengenes

trait_data <- data.frame(
  age = metadata_ord$age,
  sex_numeric = as.numeric(metadata_ord$sex == "Male"),
  NR3C2 = nr3c2_log2,
  NR3C1 = nr3c1_log2
)

module_trait_cor <- cor(module_eigengenes, trait_data, use = "p")
module_trait_pvalue <- corPvalueStudent(module_trait_cor, nSamples = nrow(wgcna_input))

text_matrix <- paste0(round(module_trait_cor, 2), "\n(",
                       formatC(module_trait_pvalue, format = "e", digits = 1), ")")
dim(text_matrix) <- dim(module_trait_cor)

png("figures/12_wgcna_module_trait_heatmap.png", width = 9, height = 10, units = "in", res = 300)
par(mar = c(6, 8, 3, 3))
labeledHeatmap(Matrix = module_trait_cor,
               xLabels = colnames(module_trait_cor),
               yLabels = rownames(module_trait_cor),
               ySymbols = rownames(module_trait_cor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = text_matrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.y = 0.7,
               zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

# Finding: no module significantly associated with age (all |r|<=0.16, all p>0.05) -
#   consistent with the diffuse, low-magnitude age effect at the individual-gene level.
# NR3C1 strongly correlated with several modules (e.g. MEred r=0.69 p=8.3e-13, MEblue r=-0.65 p=3.3e-11).
# NR3C2's strongest modules: MEturquoise (r=-0.38, p=4.6e-4), MEblue (r=-0.35, p=1.2e-3).

# -----------------------------------------------------------------------------
# 6. Module membership vs. gene significance (hub gene identification)
# -----------------------------------------------------------------------------
module_membership <- as.data.frame(cor(wgcna_input, module_eigengenes, use = "p"))
names(module_membership) <- paste0("MM_", names(module_eigengenes))

gene_significance_nr3c2 <- as.data.frame(cor(wgcna_input, trait_data$NR3C2, use = "p"))
names(gene_significance_nr3c2) <- "GS_NR3C2"

gene_significance_nr3c1 <- as.data.frame(cor(wgcna_input, trait_data$NR3C1, use = "p"))
names(gene_significance_nr3c1) <- "GS_NR3C1"

turquoise_genes <- dynamic_colors == "turquoise"
red_genes <- dynamic_colors == "red"
turquoise_gene_ids <- colnames(wgcna_input)[turquoise_genes]

df_turquoise <- data.frame(
  ensembl_gene_id = turquoise_gene_ids,
  MM = abs(module_membership$MM_MEturquoise[turquoise_genes]),
  GS = abs(gene_significance_nr3c2$GS_NR3C2[turquoise_genes])
)
# Exclude NR3C2's trivial self-correlation (GS with itself = 1.0)
turquoise_no_self <- df_turquoise |> dplyr::filter(ensembl_gene_id != nr3c2_id)

df_red <- data.frame(
  MM = abs(module_membership$MM_MEred[red_genes]),
  GS = abs(gene_significance_nr3c1$GS_NR3C1[red_genes])
)

fig_mm_gs_turquoise <- ggplot(turquoise_no_self, aes(MM, GS)) +
  geom_point(alpha = 0.4, color = "#1D9E75") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = "Module membership (turquoise)", y = "Gene significance for NR3C2",
       title = paste0("Turquoise module (n=", nrow(turquoise_no_self), " genes) vs. NR3C2")) +
  theme_minimal()

fig_mm_gs_red <- ggplot(df_red, aes(MM, GS)) +
  geom_point(alpha = 0.4, color = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = "Module membership (red)", y = "Gene significance for NR3C1",
       title = paste0("Red module (n=", sum(red_genes), " genes) vs. NR3C1")) +
  theme_minimal()

ggsave("figures/13_mm_gs_turquoise.png", fig_mm_gs_turquoise, width = 6, height = 5, dpi = 300)
ggsave("figures/14_mm_gs_red.png", fig_mm_gs_red, width = 6, height = 5, dpi = 300)

cor.test(turquoise_no_self$MM, turquoise_no_self$GS)  # r=0.648, p<2.2e-16
cor.test(df_red$MM, df_red$GS)                          # r=0.902, p<2.2e-16

# -----------------------------------------------------------------------------
# 7. TOM heatmap for the red module (NR3C1-associated)
# -----------------------------------------------------------------------------
red_indices <- which(colnames(wgcna_input) %in% colnames(wgcna_input)[red_genes])
tom_red <- tom_matrix[red_indices, red_indices]

plot_tom_red <- (1 - tom_red)^7
diag(plot_tom_red) <- NA
gene_tree_red <- hclust(as.dist(1 - tom_red), method = "average")

png("figures/15_tom_heatmap_red.png", width = 8, height = 8, units = "in", res = 300)
TOMplot(plot_tom_red, gene_tree_red,
        col = colorRampPalette(c("white", "firebrick"))(50),
        main = "TOM heatmap: red module (NR3C1-associated)")
dev.off()

# Finding: no sharp block-diagonal substructure - smooth continuous gradient instead.
# Consistent with the unexplained PC1 batch-like variance (~49%) flagged during QC;
# red module may partly reflect a broad continuous axis rather than discrete sub-programs.

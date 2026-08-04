# NR3C2 (Mineralocorticoid Receptor) Expression in Aging Human Skin Fibroblasts

An in silico bioinformatics analysis of NR3C2 (mineralocorticoid receptor, MR) expression
and its relationship to the aging transcriptome in human dermal fibroblasts.

**Dataset:** GSE226189 (Tsitsipatis et al., 2023, *Aging Cell*), n=82 primary skin fibroblast
samples, ages 22-89.

**Thesis:** Near East University, Medical Genetics and Biology Programme (Master with
Thesis). Supervisor: Seniye Targen; Co-supervisor: Tutku Yaras.

## Pipeline overview

Scripts are numbered in the order they should be run. Each depends on outputs from earlier
scripts (raw counts and metadata are cached as .rds files in data/ after the first run).

| Script | Description |
|---|---|
| 01_download_and_prep.R | Downloads GSE226189, aligns count matrix to metadata, filters to protein-coding mRNAs |
| 02_qc_and_deseq2.R | Low-count filtering, DESeq2 linear and spline (df=3) models |
| 03_figures.R | PCA/QC, volcano plot, top up/down heatmap, spline-specific gene trend curves |
| 04_nr3c2_analysis.R | NR3C2-specific age trend and sex-stratified analysis |
| 05_candidate_panel.R | 34-gene MR/GR signaling panel: correlation network, category modules, IGF1/LOX/SGK1, HSD11B2 selectivity |
| 06_wgcna.R | WGCNA: module detection, module-trait relationships, hub gene identification |
| 07_appendix_crosscheck.R | Cross-validates NR3C2 findings against the original study's supplementary DE gene lists |
| 08_gsea.R | GSEA on four proposal themes: ECM organization, senescence, oxidative stress, autophagy |

## Key findings

- NR3C2 shows a nominal, non-FDR-significant increase with age (linear model p=0.0096,
  spline model p=0.0013), stronger in female samples (r=0.46) than male (r=0.25).
- This trend was independently corroborated by the original study's own supplementary results.
- NR3C1 (glucocorticoid receptor) showed consistent, literature-aligned transcriptional
  relationships across a 34-gene candidate panel and WGCNA modules; NR3C2 showed comparatively
  weak and inconsistent relationships throughout.
- HSD11B2 is expressed at very low levels relative to both receptors, raising the possibility
  that MR selectivity for aldosterone over cortisol is not strongly enforced in dermal fibroblasts.
- GSEA identified macroautophagy as the only significantly enriched pathway among the four
  proposal themes (NES=1.57, padj=0.0004).

## Requirements

R packages: tidyverse, DESeq2, splines, biomaRt, pheatmap, igraph, ggraph, WGCNA, readxl,
fgsea, msigdbr, patchwork, ggrepel.

## Data availability

Raw counts and metadata are downloaded programmatically from GEO (GSE226189) in
01_download_and_prep.R and are not stored in this repository. The original study's
supplementary Appendix S1 file is available from the publisher (Tsitsipatis et al., 2023,
*Aging Cell*, DOI: 10.1111/acel.13915) and is not redistributed here.

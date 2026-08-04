# ==============================================================================
# 01_download_and_prep.R
# NR3C2 (Mineralocorticoid Receptor) expression in aging skin fibroblasts
# Dataset: GSE226189 (Tsitsipatis et al. 2023), n=82, ages 22-89
# ==============================================================================

library(GEOquery)
library(tidyverse)
library(DESeq2)

# --- Metadata (already downloaded) ---
metadata <- read.csv("data/GSE226189_pheno_data.csv", row.names = 1)
dim(metadata)
head(metadata)
# --- Locate raw count files ---
count_files <- list.files("data/GSE226189", pattern = "geneCOUNT.txt.gz$",
                          recursive = TRUE, full.names = TRUE)
length(count_files)  # should be 82
# --- Read and combine count files ---
read_one_sample <- function(path) {
  read.table(path, header = TRUE, sep = "\t", col.names = c("gene_id", "count"))
}

count_list <- lapply(count_files, read_one_sample)

sample_names <- basename(count_files) |>
  str_remove("^GSM[0-9]+_") |>
  str_remove("_geneCOUNT\\.txt\\.gz$")

gene_ids_check <- lapply(count_list, function(df) df$gene_id)
stopifnot(all(sapply(gene_ids_check, identical, gene_ids_check[[1]])))

count_matrix <- sapply(count_list, function(df) df$count)
rownames(count_matrix) <- count_list[[1]]$gene_id
colnames(count_matrix) <- sample_names

dim(count_matrix)  # should be 57773 x 82

# --- Clean & align metadata (metadata rows are GSM-indexed; count_matrix
# columns use the "title" naming scheme, so join on title) ---
metadata_clean <- metadata |>
  mutate(
    sample_id = title,
    sex = str_extract(characteristics_ch1.1, "Male|Female"),
    age = as.numeric(str_extract(characteristics_ch1.2, "[0-9]+")),
    # Young/old cutoffs match the validation cohorts used in Tsitsipatis et al. 2023
    age_group = case_when(
      age < 35 ~ "young",
      age > 70 ~ "old",
      TRUE ~ "middle"
    )
  ) |>
  dplyr::select(sample_id, sex, age, age_group)

stopifnot(setequal(metadata_clean$sample_id, colnames(count_matrix)))
metadata_clean <- metadata_clean[match(colnames(count_matrix), metadata_clean$sample_id), ]
rownames(metadata_clean) <- metadata_clean$sample_id
stopifnot(identical(metadata_clean$sample_id, colnames(count_matrix)))

# NOTE: parsed sex splits as 36 F / 46 M, vs. 35 F / 47 M reported in the
# paper's Table 1 -- a one-sample discrepancy not resolvable from GEO
# metadata alone. Proceeding with GEO-derived values.

# --- Restrict to protein-coding genes (mRNAs only) ---
# This dataset also contains lncRNA/circRNA quantifications (see
# GSE226189 raw files); we filter the gene-level count matrix down to
# protein-coding Ensembl gene IDs using biomaRt biotype annotation.
library(biomaRt)

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

biotype_map <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(count_matrix),
  mart = mart
)

gene_symbol_map <- biotype_map |>
  filter(gene_biotype == "protein_coding") |>
  distinct(ensembl_gene_id, hgnc_symbol)

protein_coding_ids <- gene_symbol_map$ensembl_gene_id

count_matrix_mrna <- count_matrix[rownames(count_matrix) %in% protein_coding_ids, ]
dim(count_matrix_mrna)  # 19225 x 82

saveRDS(count_matrix_mrna, "data/count_matrix_mrna.rds")
saveRDS(metadata_clean, "data/metadata_clean.rds")
saveRDS(gene_symbol_map, "data/gene_symbol_map.rds")
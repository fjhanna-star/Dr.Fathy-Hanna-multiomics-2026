# Module 5, Exercise 1, Step 1: rank genes, run fgsea on Hallmark
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(fgsea)
  library(msigdbr)
})

genes_by_outcome <- read_csv(here("data", "genes_by_outcome.csv"), show_col_types = FALSE)

# 1. Ranking vector: stat, named by Ensembl ID, missing values dropped
stat_named <- genes_by_outcome$stat
names(stat_named) <- genes_by_outcome$feature
stat_named <- stat_named[!is.na(stat_named)]
cat("Genes in the ranking:", length(stat_named), "\n")

# 2. Hallmark sets, Ensembl IDs
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
cat("msigdbr rows:", nrow(hallmark), " sets:", n_distinct(hallmark$gs_name), "\n")
sets <- split(hallmark$ensembl_gene, hallmark$gs_name)
cat("Genes in the ranking that appear in at least one set:",
    sum(names(stat_named) %in% unique(hallmark$ensembl_gene)), "\n")
cat("NA Ensembl IDs inside the sets:", sum(is.na(hallmark$ensembl_gene)), "\n")

# 3. fgsea
set.seed(1)
fg <- fgsea(sets, stat_named, minSize = 15, maxSize = 500) |>
  as_tibble() |>
  arrange(padj)
cat("Sets tested:", nrow(fg), "\n")
cat("Significant (padj < 0.05):", sum(fg$padj < 0.05), "\n")
cat("  higher in PRETERM (NES > 0):", sum(fg$padj < 0.05 & fg$NES > 0), "\n")
cat("  higher in TERM    (NES < 0):", sum(fg$padj < 0.05 & fg$NES < 0), "\n")

fg |>
  filter(padj < 0.05) |>
  select(pathway, NES, padj, size, leadingEdge) |>
  mutate(n_leading = map_int(leadingEdge, length)) |>
  select(-leadingEdge) |>
  print(n = 30)

saveRDS(fg, here("data", "m5_fgsea_hallmark_outcome.rds"))

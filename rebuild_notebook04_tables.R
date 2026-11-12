# Rebuilds the three result tables Module 5 needs, by re-running the
# analyses from notebook 04 (and the grouping column from notebook 02).
#
# Run: Rscript rebuild_notebook04_tables.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(DESeq2)
})

source(here("R/count_tables.R"))
d <- read_counts("instructional")

meta <- d$metadata
week18 <- meta$sample_id[meta$time_point == "week_18"]
cat("Week 18 samples:", length(week18), "\n")

# ---------------------------------------------------------------------------
# 1. The community grouping (notebook 02, Ex 2 step 5)
#    Lactobacillus fraction of reads per sample, from the FULL species table.
# ---------------------------------------------------------------------------
lacto_species <- d$taxonomy$taxon[d$taxonomy$Genus == "Lactobacillus"]
cat("Species counted as Lactobacillus:", length(lacto_species), "\n")

lacto_counts <- colSums(d$species[lacto_species, week18, drop = FALSE])
totals       <- colSums(d$species[, week18, drop = FALSE])
lacto_frac   <- lacto_counts / totals

sample_groups <- tibble(
  sample_id      = week18,
  lacto_fraction = as.numeric(lacto_frac),
  community      = ifelse(lacto_frac >= 0.5, "lactobacillus", "anaerobe")
)
write_csv(sample_groups, here("data", "my_sample_groups.csv"))
cat("Community counts:\n")
print(table(sample_groups$community))

# ---------------------------------------------------------------------------
# 2. Genes by community (notebook 04, Ex 1 step 1)
#    Week 18, genes with >=10 counts in >=25% of samples, ~ community.
# ---------------------------------------------------------------------------
m18 <- meta |> filter(time_point == "week_18") |>
  left_join(sample_groups, by = "sample_id")

counts_g <- d$genes[, m18$sample_id]
keep_g <- rowSums(counts_g >= 10) >= ceiling(ncol(counts_g) / 4)
cat("Genes passing the filter:", sum(keep_g), "\n")

dds_g <- DESeqDataSetFromMatrix(
  countData = counts_g[keep_g, ],
  colData   = m18,
  design    = ~ community
)
dds_g <- DESeq(dds_g)
cat("Community levels:", paste(levels(dds_g$community), collapse = ", "), "\n")

res_gc <- results(dds_g, contrast = c("community", "anaerobe", "lactobacillus"))
cat("padj < 0.05:", sum(res_gc$padj < 0.05, na.rm = TRUE), "\n")

genes_by_community <- as.data.frame(res_gc) |>
  rownames_to_column("feature") |>
  left_join(d$gene_info, by = c("feature" = "ensembl_gene_id")) |>
  select(feature, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
         symbol, gene_name)
write_csv(genes_by_community, here("data", "genes_by_community.csv"))
cat("Saved genes_by_community.csv:", nrow(genes_by_community), "rows\n")

# ---------------------------------------------------------------------------
# 3. Species by outcome (notebook 04, Ex 2 step 3, poscounts version)
#    Week 18, species present in >=10% of samples, ~ outcome, poscounts.
# ---------------------------------------------------------------------------
counts_s <- d$species[, m18$sample_id]
keep_s <- rowSums(counts_s > 0) >= 0.1 * ncol(counts_s)
cat("Species passing the 10% presence filter:", sum(keep_s), "\n")

dds_s <- DESeqDataSetFromMatrix(
  countData = counts_s[keep_s, ],
  colData   = m18,
  design    = ~ outcome
)
dds_s <- DESeq(dds_s, sfType = "poscounts")
res_s <- results(dds_s, contrast = c("outcome", "preterm", "term"))
cat("Species padj < 0.05:", sum(res_s$padj < 0.05, na.rm = TRUE), "\n")

species_by_outcome <- as.data.frame(res_s) |>
  rownames_to_column("feature") |>
  left_join(d$taxonomy, by = c("feature" = "taxon")) |>
  select(feature, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
         Genus, Family)
write_csv(species_by_outcome, here("data", "species_by_outcome.csv"))
cat("Saved species_by_outcome.csv:", nrow(species_by_outcome), "rows\n")

# ---------------------------------------------------------------------------
# 4. Genes by outcome (notebook 04, Ex 2 step 5)
#    Same gene filter as step 2, ~ outcome, preterm vs term.
# ---------------------------------------------------------------------------
dds_o <- DESeqDataSetFromMatrix(
  countData = counts_g[keep_g, ],
  colData   = m18,
  design    = ~ outcome
)
dds_o <- DESeq(dds_o)
res_o <- results(dds_o, contrast = c("outcome", "preterm", "term"))
cat("Genes padj < 0.05:", sum(res_o$padj < 0.05, na.rm = TRUE), "\n")

genes_by_outcome <- as.data.frame(res_o) |>
  rownames_to_column("feature") |>
  left_join(d$gene_info, by = c("feature" = "ensembl_gene_id")) |>
  select(feature, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
         symbol, gene_name)
write_csv(genes_by_outcome, here("data", "genes_by_outcome.csv"))
cat("Saved genes_by_outcome.csv:", nrow(genes_by_outcome), "rows\n")

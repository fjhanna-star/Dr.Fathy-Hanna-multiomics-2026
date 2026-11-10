# Module 3: Distances and Ordination
# Author: Dr. Fathy Hanna

library(phyloseq)
library(ggplot2)
library(vegan)

# Load Data
source("R/count_tables.R")
ps <- as_counts("instructional", table = "microbiome")

# 1. Table Anatomy
ntaxa(ps)
nsamples(ps)
summary(sample_sums(ps))

# 2. PCoA & PERMANOVA
ord_bray <- ordinate(ps, method = "PCoA", distance = "bray")
df_bray <- as.data.frame(ord_bray$vectors[, 1:2])
meta_df <- data.frame(sample_data(ps))
df_bray$outcome <- as.factor(meta_df$outcome)

dist_bray <- phyloseq::distance(ps, method = "bray")

print("--- PERMANOVA: Outcome ---")
print(adonis2(dist_bray ~ outcome, data = meta_df, permutations = 999))

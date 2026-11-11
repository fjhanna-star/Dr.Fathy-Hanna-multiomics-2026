# Study B: Quality Control (QC) Report
# Analyst: Fathy Jacob Fathy Hanna

# 1. Setup environment
library(phyloseq)
library(ggplot2)

# 2. Load Data
source("R/count_tables.R")
ps <- as_counts("group_b_contraception", table = "microbiome")
meta_df <- data.frame(sample_data(ps))

# 3. QC: Library Sizes
print("--- Library Sizes Summary ---")
print(summary(sample_sums(ps)))

# 4. QC: RNA Batch vs Treatment Arms
print("--- RNA Batch vs Arm ---")
print(table(meta_df$rnaseq_batch, meta_df$arm))

# 5. Exploratory Analysis: Alpha Diversity (Shannon)
alpha_div <- estimate_richness(ps, measures = "Shannon")
meta_df$Shannon <- alpha_div$Shannon

# 6. Generate the Plot
p_alpha <- ggplot(meta_df, aes(x = arm, y = Shannon, fill = arm)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Microbiome Alpha Diversity (Shannon Index) - Study B",
       x = "Contraceptive Arm",
       y = "Shannon Index") +
  theme_minimal() +
  scale_fill_manual(values = c("injectable" = "salmon", "pill" = "springgreen", "ring" = "cornflowerblue"))

print(p_alpha)

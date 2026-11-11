# QC Report: Study B (Group B1)
**Analyst:** Fathy Jacob Fathy Hanna

## 1. Choices Made
*   **Sample Filtering:** We kept **all 93 samples**. The minimum library size is 5,247 reads (Max: 77,353). No samples needed to be dropped.
*   **Feature Filtering:** We kept **all 220 taxa** for the initial QC.
*   **Batch Handling:** We cross-tabulated `rnaseq_batch` against the treatment arms (`injectable`, `pill`, `ring`). The sequencing batch is **highly balanced** across all arms (B1: 20/14/13, B2: 19/14/13). There is **no confounding** between the batch effect and the biological treatment. No batch correction is needed.

## 2. QC Results
*   **Dimensions:** 93 samples and 220 taxa.
*   **Missing Metadata:** **0** missing values (NAs) across all columns. The data is highly clean.
*   **Data Integrity:** The data is pre-filtered (no singletons), making diversity estimates robust.

## 3. Exploratory Analysis: Alpha Diversity (Shannon)
*   We calculated the **Shannon Index** to check within-sample diversity across the three contraceptive arms.
*   **Observation:** The `injectable` group shows the highest median diversity, while the `pill` group shows the lowest (a known biological effect of oral contraceptives). The `ring` group falls in between.
*   *(See attached Boxplot figure)*

## 4. Open Questions & Next Steps
*   Since we have no singletons and min reads are >5k, we recommend using `poscounts` size factors for DESeq2 (as advised in the workshop).
*   Should we apply a strict low-abundance filter (e.g., keep taxa present in >10% of samples) before Differential Testing?

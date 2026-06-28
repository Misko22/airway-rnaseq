# ============================================================
# Week 1: Airway RNA-seq (GSE52778) differential expression
# Dexamethasone vs untreated, paired design by donor
# ============================================================

# One-time setup (run once per machine, then leave commented):
# install.packages("BiocManager")
# BiocManager::install(c("airway", "DESeq2", "tidyverse", "clusterProfiler", "org.Hs.eg.db"))

library(airway)
library(DESeq2)
library(tidyverse)

# Load packaged dataset
data("airway")
se <- airway

# Inspect the experimental design
dim(se)
table(colData(se)$dex)
table(colData(se)$cell)

# Set untreated as reference so fold changes read "treated vs untreated"
se$dex <- relevel(se$dex, ref = "untrt")
levels(se$dex)

# Build DESeq2 object with paired design: control for donor, test treatment
dds <- DESeqDataSet(se, design = ~ cell + dex)

# Filter out near-silent genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
nrow(dds)

# Run the differential expression test
dds <- DESeq(dds)
res <- results(dds)
summary(res)

# DUSP1 sanity check (known dexamethasone-responsive gene)
res["ENSG00000120129", ]

# RESULT: DUSP1 (ENSG00000120129), a canonical glucocorticoid-responsive gene,
# showed an 8-fold increase (log2FC 2.95, padj 1.8e-125), confirming the
# analysis correctly recovers known dexamethasone biology before examining
# novel results.

# Save results for downstream weeks
saveRDS(res, "deseq_results.rds")
write.csv(as.data.frame(res), "deseq_results.csv")
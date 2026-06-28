install.packages("BiocManager")
BiocManager::install(c("airway", "DESeq2", "tidyverse", "clusterProfiler", "org.Hs.eg.db"))

library(airway)
library(DESeq2)
library(tidyverse)

data("airway")
se <- airway

dim(se)
table(colData(se)$dex)
table(colData(se)$cell)

se$dex <- relevel(se$dex, ref = "untrt")
levels(se$dex)

dds <- DESeqDataSet(se, design = ~ cell + dex)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
nrow(dds)


dds <- DESeq(dds)
res <- results(dds)
summary(res)

# DUSP1 sanity check (known dexamethasone-responsive gene)
res["ENSG00000120129", ]



# RESULT: DUSP1 (ENSG00000120129), a canonical glucocorticoid-responsive gene,
# showed an 8-fold increase (log2FC 2.95, padj 1.8e-125), confirming the
# analysis correctly recovers known dexamethasone biology before examining
# novel results.
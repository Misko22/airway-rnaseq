# Week 2 - Exploring DESeq2 results
# Airway RNA-seq differential expression analysis

library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

# --------------------------------------------------
# 1. Load DESeq2 results
# --------------------------------------------------

res <- readRDS("deseq_results.rds")
res_df <- as.data.frame(res)

head(res_df)


# --------------------------------------------------
# 2. Keep statistically significant genes
# padj < 0.1
# --------------------------------------------------

sig_res <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.1,
]


# --------------------------------------------------
# 3. Rank genes by log2 fold change
# --------------------------------------------------

top_up <- sig_res[
  order(sig_res$log2FoldChange, decreasing = TRUE),
]

top_down <- sig_res[
  order(sig_res$log2FoldChange),
]


# --------------------------------------------------
# 4. Add human-readable gene annotations
# --------------------------------------------------

top_up$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(top_up),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

top_up$gene_name <- mapIds(
  org.Hs.eg.db,
  keys = rownames(top_up),
  column = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)

head(
  top_up[, c(
    "gene_symbol",
    "gene_name",
    "baseMean",
    "log2FoldChange",
    "padj"
  )],
  10
)


# --------------------------------------------------
# 5. Inspect ALOX15B across paired samples
# --------------------------------------------------

alox_counts <- plotCounts(
  dds,
  gene = "ENSG00000179593",
  intgroup = c("dex", "cell"),
  returnData = TRUE
)

alox_counts

ggplot(
  alox_counts,
  aes(
    x = dex,
    y = count,
    group = cell
  )
) +
  geom_line() +
  geom_point(size = 3) +
  scale_y_log10() +
  labs(
    title = "ALOX15B expression after dexamethasone",
    x = "Treatment",
    y = "Normalized count",
    caption = "Paired donor-derived airway smooth muscle cell samples"
  ) +
  theme_minimal()


# --------------------------------------------------
# 6. Volcano plot
# --------------------------------------------------

volcano_df <- res_df

volcano_df$significance <- "Not significant"

volcano_df$significance[
  !is.na(volcano_df$padj) &
    volcano_df$padj < 0.1 &
    volcano_df$log2FoldChange > 0
] <- "Upregulated"

volcano_df$significance[
  !is.na(volcano_df$padj) &
    volcano_df$padj < 0.1 &
    volcano_df$log2FoldChange < 0
] <- "Downregulated"

ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = significance
  )
) +
  geom_point(
    alpha = 0.5,
    size = 1.5
  ) +
  labs(
    title = "Dexamethasone response in airway smooth muscle cells",
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    color = "Gene classification"
  ) +
  theme_minimal()


# --------------------------------------------------
# 7. Rank upregulated genes by statistical evidence
# --------------------------------------------------

top_evidence_up <- top_up[
  order(top_up$padj),
]

head(
  top_evidence_up[, c(
    "gene_symbol",
    "gene_name",
    "baseMean",
    "log2FoldChange",
    "padj"
  )],
  10
)
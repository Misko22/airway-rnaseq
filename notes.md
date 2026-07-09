# Airway RNA-seq (GSE52778) - working notes

## Week 1: setup and first result
- Dataset: airway smooth muscle, 4 donors, dexamethasone vs untreated, paired design
- Design: ~ cell + dex (controls for donor, tests treatment)
- 22,369 genes after filtering low counts
- At padj < 0.1: 2,610 up, 2,224 down

### Sanity check (the spine of the interpretation)
DUSP1, a canonical glucocorticoid-responsive gene, showed an 8-fold increase
(log2FC 2.95, padj 1.8e-125), confirming the analysis correctly recovers known
dexamethasone biology before examining novel results.



## Week 2 - Reading the differential expression results

I explored the DESeq2 results beyond the initial significant gene counts.

### Understanding the DESeq2 result table

I reviewed the main DESeq2 result columns:

- baseMean: average normalized RNA-seq signal across samples
- log2FoldChange: estimated direction and size of expression change
- lfcSE: uncertainty around the log2 fold-change estimate
- stat: Wald statistic, comparing estimated effect with its uncertainty
- pvalue: raw statistical evidence against the no-treatment-effect model
- padj: p-value adjusted for testing thousands of genes

A negative log2 fold change alone does not mean a gene is significantly downregulated. Statistical significance must also be evaluated using the adjusted p-value threshold.

### ALOX15B

ALOX15B had the largest positive log2 fold change among the significant genes:

- log2FoldChange: approximately 9.5
- estimated fold change: approximately 727-fold
- adjusted p-value: approximately 1.18e-17

Sample-level normalized counts were inspected using plotCounts().

All four donor-derived airway smooth muscle cell samples showed a strong increase in ALOX15B expression after dexamethasone treatment. The magnitude of the response varied between donors, but the direction of change was consistent.

### Volcano plot

I created a volcano plot to visualize differential expression.

The x-axis represents log2 fold change:
- negative values indicate decreased expression after treatment
- positive values indicate increased expression after treatment

The y-axis represents -log10 adjusted p-value. Genes higher on the plot have stronger statistical evidence for differential expression.

This analysis helped distinguish effect size from statistical evidence.

ALOX15B ranked highest by positive log2 fold change, representing the largest estimated increase in expression.

SPARCL1 ranked highest among upregulated genes when sorted by adjusted p-value, representing the strongest adjusted statistical evidence.

DUSP1 ranked third by adjusted p-value among the upregulated genes. Its known glucocorticoid responsiveness provides an additional biological sanity check for the analysis.

### Main lesson

The largest estimated expression change is not necessarily the gene with the strongest statistical evidence.

Differential expression results should be interpreted using effect size, uncertainty, statistical evidence, and sample-level expression patterns together.
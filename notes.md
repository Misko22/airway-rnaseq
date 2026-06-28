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
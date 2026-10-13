# ============================================================
# Week 3: What does "important" mean in gene expression data?
# Focus: scientific data thinking + communication
# ============================================================

library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ============================================================
# Visual 1: From thousands of genes to a story
# ============================================================

overview <- data.frame(
  stage = c(
    "Starting gene entries",
    "Analyzed after filtering",
    "Significantly changed",
    "Upregulated",
    "Downregulated"
  ),
  genes = c(
    63677,
    22369,
    4834,
    2610,
    2224
  )
)

overview


# ------------------------------------------------------------
# Visual 1: analysis story
# ------------------------------------------------------------

flow_df <- data.frame(
  x = c(0, 0, 0, -0.9, 0.9),
  y = c(5, 4, 3, 2, 2),
  
  number = c(
    "63,677",
    "22,369",
    "4,834",
    "2,610",
    "2,224"
  ),
  
  label = c(
    "starting gene entries",
    "retained for analysis",
    "significant at adjusted p < 0.1",
    "upregulated",
    "downregulated"
  )
)

ggplot(flow_df, aes(x = x, y = y)) +
  
  # Main vertical flow
  annotate(
    "segment",
    x = 0,
    xend = 0,
    y = 4.5,
    yend = 4.3,
    linewidth = 0.5
  ) +
  
  annotate(
    "segment",
    x = 0,
    xend = 0,
    y = 3.5,
    yend = 3.3,
    linewidth = 0.5
  ) +
  
  # Split into upregulated and downregulated groups
  annotate(
    "segment",
    x = 0,
    xend = -0.9,
    y = 2.55,
    yend = 2.3,
    linewidth = 0.5
  ) +
  
  annotate(
    "segment",
    x = 0,
    xend = 0.9,
    y = 2.55,
    yend = 2.3,
    linewidth = 0.5
  ) +
  
  # Large numbers
  geom_text(
    aes(label = number),
    size = 8,
    fontface = "bold"
  ) +
  
  # Smaller explanations
  geom_text(
    aes(
      label = label,
      y = y - 0.28
    ),
    size = 4
  ) +
  
  coord_cartesian(
    xlim = c(-1.7, 1.7),
    ylim = c(1.4, 5.5),
    clip = "off"
  ) +
  
  labs(
    title = "From 63,677 gene entries to 4,834 significant changes",
    subtitle = "Dexamethasone response in airway smooth muscle cells"
  ) +
  
  theme_void() +
  
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5,
      margin = margin(b = 20)
    ),
    
    plot.margin = margin(
      20,
      30,
      20,
      30
    )
  )


# ============================================================
# Visual 2: "Top" depends on the question
# ============================================================

# Load previously saved DESeq2 results
res <- readRDS("deseq_results.rds")
res_df <- as.data.frame(res)


# ------------------------------------------------------------
# Add human-readable gene symbols
# ------------------------------------------------------------

res_df$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res_df),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)


# ------------------------------------------------------------
# Keep significant upregulated genes
# ------------------------------------------------------------

up_genes <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.1 &
    res_df$log2FoldChange > 0,
]


# ------------------------------------------------------------
# Ask three different questions
# ------------------------------------------------------------

# 1. Which gene had the largest estimated increase?
largest_change <- up_genes[
  which.max(up_genes$log2FoldChange),
]


# 2. Which upregulated gene had the strongest statistical evidence?
strongest_evidence <- up_genes[
  which.min(up_genes$padj),
]


# 3. Which known dexamethasone-response gene can act as a biological check?
dusp1 <- res_df[
  !is.na(res_df$gene_symbol) &
    res_df$gene_symbol == "DUSP1",
]


# ------------------------------------------------------------
# Compare the three results
# ------------------------------------------------------------

comparison <- data.frame(
  reason = c(
    "Largest estimated increase",
    "Strongest statistical evidence",
    "Known biological response"
  ),
  
  gene = c(
    largest_change$gene_symbol,
    strongest_evidence$gene_symbol,
    dusp1$gene_symbol
  ),
  
  log2FoldChange = c(
    largest_change$log2FoldChange,
    strongest_evidence$log2FoldChange,
    dusp1$log2FoldChange
  ),
  
  padj = c(
    largest_change$padj,
    strongest_evidence$padj,
    dusp1$padj
  )
)

comparison


# ------------------------------------------------------------
# Convert log2 fold change into easier-to-read fold change
# ------------------------------------------------------------

comparison$fold_change <- 2^comparison$log2FoldChange

comparison


# ============================================================
# Visual 2: same data, different definition of "top"
# ============================================================

cards <- data.frame(
  x = c(1, 2, 3),
  
  title = c(
    "BIGGEST CHANGE",
    "STRONGEST EVIDENCE",
    "KNOWN RESPONSE"
  ),
  
  gene = comparison$gene,
  
  main_value = c(
    paste0(
      "~",
      round(comparison$fold_change[1]),
      "× increase"
    ),
    
    paste0(
      "padj = ",
      format(
        comparison$padj[2],
        scientific = TRUE,
        digits = 2
      )
    ),
    
    paste0(
      "~",
      round(comparison$fold_change[3], 1),
      "× increase"
    )
  ),
  
  detail = c(
    paste0(
      "log2FC = ",
      round(
        comparison$log2FoldChange[1],
        2
      )
    ),
    
    paste0(
      "log2FC = ",
      round(
        comparison$log2FoldChange[2],
        2
      )
    ),
    
    "known dexamethasone response"
  )
)


ggplot(cards, aes(x = x, y = 1)) +
  
  geom_tile(
    width = 0.85,
    height = 1.7,
    fill = "white",
    color = "grey70",
    linewidth = 0.6
  ) +
  
  geom_text(
    aes(
      label = title,
      y = 1.55
    ),
    size = 4,
    fontface = "bold"
  ) +
  
  geom_text(
    aes(
      label = gene,
      y = 1.20
    ),
    size = 7,
    fontface = "bold"
  ) +
  
  geom_text(
    aes(
      label = main_value,
      y = 0.85
    ),
    size = 5
  ) +
  
  geom_text(
    aes(
      label = detail,
      y = 0.55
    ),
    size = 3.8
  ) +
  
  annotate(
    "text",
    x = 2,
    y = -0.15,
    label = "Same dataset. Different question. Different \"top\" gene.",
    size = 4.5,
    fontface = "italic"
  ) +
  
  coord_cartesian(
    xlim = c(0.4, 3.6),
    ylim = c(-0.4, 2),
    clip = "off"
  ) +
  
  labs(
    title = "\"Top\" depends on the question",
    subtitle = "Three genes stand out for three different reasons"
  ) +
  
  theme_void() +
  
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5,
      margin = margin(b = 20)
    ),
    
    plot.margin = margin(
      20,
      30,
      20,
      30
    )
  )


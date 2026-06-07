# =============================================================================
# Bulk RNA-seq marker-gene expression heatmap
#
# Plots a curated, lineage-ordered panel of marker genes across bulk RNA-seq
# samples as a row-scaled (per-gene z-score) heatmap. Genes are grouped into
# lineage blocks separated by horizontal gaps.
# =============================================================================

library(edgeR)     # CPM normalisation
library(pheatmap)  # heatmap

# =============================================================================
# 1. LOAD BULK RNA-seq COUNTS
# =============================================================================
# Filtered gene-level count matrix (genes x samples), three replicates per
# condition. check.names = FALSE preserves the original sample column names
# (R would otherwise rewrite names containing spaces or starting with digits).

BulkRNA <- read.csv("data/BulkRNA_gene_counts_filtered_3repl.csv",
                    row.names = 1, check.names = FALSE)

# =============================================================================
# 2. NORMALISE: counts -> log2 CPM
# =============================================================================
# Library-size normalise to counts-per-million, then log2-transform with a
# pseudocount of 1 to stabilise variance and avoid log(0).

BulkRNA_cpm <- cpm(BulkRNA)
BulkRNA_log <- log2(BulkRNA_cpm + 1)

# =============================================================================
# 3. CURATED MARKER PANEL (lineage-ordered)
# =============================================================================
# Gene order defines the row order of the heatmap. The gaps_row values in
# section 5 split this list into lineage blocks, so the order here and those
# break positions must stay consistent (see the note in section 5).

genes_of_interest <- c(
  "Nanog","Sox2","Pou5f1","Pdgfa","Zfp42","Klf2","Tfcp2l1","Esrrb","Prdm14",
  "Jam2","Notum","Tcl1",
  "Sox17","Foxa2","Hnf1b","Hnf4a","Gata4","Gata6","Dab2","Amn","Srgn","Runx1",
  "Eomes","Pdgfra",
  "Igf1","Dusp4","Htra1","Serpinh1","P4ha2","Cyp4f14","Gbgt1","Clic3","Ddr2",
  "Sox7",
  "Dkk1","Nodal","Lefty1","Cer1","Hhex","Lhx1","Otx2","Sema6d","Cd8a","Mest",
  "Hes1","Shisa2","Rgs8","Lhx1os","Fzd7","Sfrp5","Irs4","Akr1c12","Akr1c13",
  "Akr1c19",
  "Fgf5","Nog","Fst","Thbd","Cryab","Plat","Drp2","Islr2","Fxyd3","Hs3st1",
  "Trpc5os","Cdh6","Lama1","Lamb1","Lamc1","Col4a1","Col4a2",
  "Nid1","Fabp1","Aldob","Cdhr2","Creb3l3","Ihh","Il22ra1","Slc2a2","Kif12",
  "Mbl2","Rbp4","Mttp","Gjb1","Rdh10","Apoa2","Apoc2","Ttr","Cubn","Bmp6"
)

# Report any panel genes missing from the matrix (so a typo or naming mismatch
# is caught rather than silently dropped).
present <- genes_of_interest %in% rownames(BulkRNA_log)
cat("Marker genes found:", sum(present), "/", length(genes_of_interest), "\n")
if (any(!present))
  cat("Missing:", paste(genes_of_interest[!present], collapse = ", "), "\n")

# Subset to present genes, preserving the panel order
BulkRNA_subset <- as.matrix(BulkRNA_log[genes_of_interest[present], ])

# =============================================================================
# 4. SAMPLE ANNOTATION
# =============================================================================
# Derive a Condition label from each column name by stripping the replicate
# suffix (everything from the first underscore). Adjust the regex to match your
# sample naming scheme.

sample_info <- data.frame(
  Condition = gsub("_.*", "", colnames(BulkRNA_subset)),
  row.names = colnames(BulkRNA_subset)
)

# =============================================================================
# 5. HEATMAP
# =============================================================================
# scale = "row" converts each gene to a z-score across samples, so the colour
# reflects relative expression per gene. Rows and columns are left unclustered
# to preserve the curated lineage order and sample order.
#
# IMPORTANT: gaps_row marks the lineage-block boundaries and is given as row
# positions in the *plotted* matrix. These values assume every panel gene is
# present and in the order above; if any genes are missing (see section 3), the
# downstream boundaries shift and should be adjusted accordingly.

pheatmap(
  BulkRNA_subset,
  annotation_col = sample_info,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  scale          = "row",
  gaps_row       = c(12, 24, 34, 54, 71),  # lineage-block separators
  gaps_col       = c(1),                   # separate the first sample/condition
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  fontsize_row   = 8,
  color          = colorRampPalette(c("dodgerblue4", "white", "red3"))(100),
  main           = "Bulk RNA-seq marker-gene expression",
  filename       = "figures/bulk_marker_heatmap.pdf"
)

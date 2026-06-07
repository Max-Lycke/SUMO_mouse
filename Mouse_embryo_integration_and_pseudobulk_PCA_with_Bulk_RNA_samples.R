# =============================================================================
# Integration of published mouse embryo scRNA-seq datasets and comparison of
# cell-type pseudobulk profiles to in vitro bulk RNA-seq
#
# Overview of the workflow:
#   1. Load five published mouse embryo scRNA-seq datasets
#   2. Harmonise gene identifiers (Ensembl -> MGI symbols where needed)
#   3. Per-dataset QC, normalisation and variable-feature selection
#   4. CCA-based anchor integration (Seurat v4)
#   5. Dimensionality reduction, clustering and lineage annotation
#   6. Cluster marker detection
#   7. Cell-type pseudobulk generation
#   8. Comparison of scRNA-seq pseudobulk to bulk RNA-seq (log-CPM, ComBat,
#      Pearson correlation heatmaps and PCA)
#
# Datasets:
#   Cheng et al. 2019, Cell Rep         (doi:10.1016/j.celrep.2019.02.031)
#   Mohammed et al. 2017, Cell Rep      (doi:10.1016/j.celrep.2017.07.009)
#   Liu et al. 2022, Sci Adv            (doi:10.1126/sciadv.abj3725)
#   Argelaguet et al. 2019, Nature      (doi:10.1038/s41586-019-1825-8)
#   Thowfeequ et al. 2024, Dev Cell     (doi:10.1016/j.devcel.2024.05.014)
# =============================================================================

library(Seurat)        # v4 — CCA anchor integration
library(dplyr)
library(biomaRt)
library(ggplot2)
library(edgeR)
library(sva)
library(pheatmap)

set.seed(42)

# Output directory — edit to your environment
outdir <- "results/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD THE FIVE PUBLISHED DATASETS
# =============================================================================
# Each dataset is loaded as a Seurat object from its published count matrix /
# accession (see DOIs in the header). Insert the appropriate loading code for
# each source below; the rest of the script assumes the objects exist with
# their original author-supplied cell metadata.

# Cheng_data      <- readRDS("data/Cheng.rds")
# Mohammed_data   <- readRDS("data/Mohammed.rds")
# Liu_mTE         <- readRDS("data/Liu.rds")
# Argelaguet_data <- readRDS("data/Argelaguet.rds")
# Thowfeequ       <- readRDS("data/Thowfeequ.rds")

# Tag each object with its dataset of origin (used as the batch label later)
Cheng_data$dataset      <- "Cheng"
Mohammed_data$dataset   <- "Mohammed"
Liu_mTE$dataset         <- "Liu"
Argelaguet_data$dataset <- "Argelaguet"
Thowfeequ$dataset       <- "Thowfeequ"

# =============================================================================
# 2. CONVERT ENSEMBL IDs -> GENE SYMBOLS (datasets that need it)
# =============================================================================
# Cheng, Mohammed and Liu already use gene symbols. Argelaguet and Thowfeequ
# store features as Ensembl IDs, so we map them to MGI symbols with biomaRt.
#
# NOTE: these are MOUSE data, so we query the mouse Ensembl mart and the symbols
# returned are MGI symbols (the reported methods text says "HGNC", which is the
# human nomenclature committee — that wording should read "MGI" for mouse).

convert_ensembl_to_symbols <- function(seurat_obj, species = "mouse") {
  ensembl_ids <- rownames(seurat_obj)

  dataset <- switch(species,
                    human = "hsapiens_gene_ensembl",
                    mouse = "mmusculus_gene_ensembl",
                    stop("species must be 'human' or 'mouse'"))
  mart <- useMart("ensembl", dataset = dataset)

  # Look up symbols for the Ensembl IDs present in the object
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters    = "ensembl_gene_id",
    values     = ensembl_ids,
    mart       = mart
  )
  gene_map <- gene_map[gene_map$external_gene_name != "", ]
  gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]

  mapping <- setNames(gene_map$external_gene_name, gene_map$ensembl_gene_id)
  new_names <- mapping[ensembl_ids]

  # Unmapped IDs keep their Ensembl ID; duplicated symbols get the ID appended
  new_names[is.na(new_names)] <- ensembl_ids[is.na(new_names)]
  dup <- duplicated(new_names) | duplicated(new_names, fromLast = TRUE)
  new_names[dup] <- paste0(new_names[dup], "_", ensembl_ids[dup])

  # Rebuild the object with renamed features (counts only; re-normalise later)
  counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  rownames(counts) <- new_names
  CreateSeuratObject(counts = counts, meta.data = seurat_obj@meta.data)
}

Argelaguet_data <- convert_ensembl_to_symbols(Argelaguet_data, species = "mouse")
Thowfeequ       <- convert_ensembl_to_symbols(Thowfeequ,       species = "mouse")

# CreateSeuratObject does not carry over arbitrary object-level columns, so
# re-apply the dataset tags after conversion
Argelaguet_data$dataset <- "Argelaguet"
Thowfeequ$dataset       <- "Thowfeequ"

# =============================================================================
# 3. PER-DATASET QC, NORMALISATION AND VARIABLE FEATURES
# =============================================================================
# Keep only high-complexity cells (> 5000 detected genes), then log-normalise
# and select the top 2000 variable features per dataset (VST). The strict gene
# threshold reflects the deep, well-covered nature of these embryo datasets.

datasets <- list(Cheng_data, Mohammed_data, Liu_mTE, Argelaguet_data, Thowfeequ)

datasets <- lapply(datasets, function(obj) {
  obj <- subset(obj, subset = nFeature_RNA > 5000)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj
})

# Confirm a usable shared gene space exists across all datasets
common_genes <- Reduce(intersect, lapply(datasets, rownames))
cat("Common genes across all datasets:", length(common_genes), "\n")

# =============================================================================
# 4. CCA-BASED ANCHOR INTEGRATION
# =============================================================================
# Standard Seurat v4 integration: select integration features, find anchors by
# canonical correlation analysis, then integrate. k.weight is lowered to 50
# because some datasets are small and the default (100) can exceed the number
# of anchors available for the smallest dataset.

features <- SelectIntegrationFeatures(object.list = datasets)
anchors  <- FindIntegrationAnchors(object.list = datasets,
                                   anchor.features = features)
integrated_data <- IntegrateData(anchorset = anchors, k.weight = 50)

# =============================================================================
# 5. DIMENSIONALITY REDUCTION AND CLUSTERING
# =============================================================================
# Work on the batch-corrected "integrated" assay for embedding/clustering only.

DefaultAssay(integrated_data) <- "integrated"

integrated_data <- ScaleData(integrated_data, verbose = FALSE)
integrated_data <- RunPCA(integrated_data, npcs = 50, verbose = FALSE)
ElbowPlot(integrated_data, ndims = 50)   # inspect to justify the PCs used below

integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:30)
integrated_data <- FindNeighbors(integrated_data, dims = 1:30)
integrated_data <- FindClusters(integrated_data, resolution = 0.1)

# Integration sanity check — datasets should intermix
DimPlot(integrated_data, group.by = "dataset", label = TRUE, pt.size = 1)

# =============================================================================
# 6. LINEAGE ANNOTATION
# =============================================================================
# Each source dataset carries its own lineage annotation under a differently
# named metadata column. We pull those through onto the integrated object, then
# (see note) manually curate and refine them using marker expression.
#
# Column name per dataset (order matches `datasets`):
#   Cheng -> Lineage1, Mohammed -> Lineage, Liu -> CellType,
#   Argelaguet -> Lineage1, Thowfeequ -> Lineage

all_labels <- c(as.character(datasets[[1]]$Lineage1),
                as.character(datasets[[2]]$Lineage),
                as.character(datasets[[3]]$CellType),
                as.character(datasets[[4]]$Lineage1),
                as.character(datasets[[5]]$Lineage))

integrated_data$Cell_type <- all_labels
Idents(integrated_data)   <- "Cell_type"

# ---- Manual curation -------------------------------------------------------
# Annotations were refined interactively (marker-guided relabelling, e.g. with
# CellSelector) to produce the final cell-type column ("Update"). Because that
# step is manual and not reproducible from code alone, the final curated
# annotations are provided as a supplementary table and reloaded here.

curated <- read.csv(file.path(outdir, "integrated_metadata_curated.csv"),
                    row.names = 1)
curated <- curated[colnames(integrated_data), ]          # match cell order
integrated_data <- AddMetaData(integrated_data, curated) # adds the "Update" column
Idents(integrated_data) <- "Update"

# =============================================================================
# 7. CLUSTER MARKER GENES
# =============================================================================
# Positive markers per annotated cell type, for characterisation / supplement.
# Marker detection uses the RNA assay (not the integrated assay).

DefaultAssay(integrated_data) <- "RNA"

markers <- FindAllMarkers(integrated_data,
                          only.pos = TRUE,
                          logfc.threshold = 0.25,
                          min.pct = 0.1,
                          min.diff.pct = 0.1)
write.csv(markers, file.path(outdir, "integrated_markers.csv"), row.names = TRUE)

top200 <- markers %>% group_by(cluster) %>% slice_head(n = 200)
write.csv(top200, file.path(outdir, "integrated_markers_top200.csv"),
          row.names = TRUE)

# =============================================================================
# 8. VISUALISATION OF THE INTEGRATED OBJECT
# =============================================================================

DimPlot(integrated_data, reduction = "umap", label = TRUE, repel = TRUE,
        pt.size = 1) +
  ggtitle("Integrated mouse embryo datasets")
ggsave(file.path(outdir, "UMAP_integrated_mouse.pdf"), width = 12, height = 10)

# Lineage-ordered marker dot plot. `lineage_markers` is the curated gene panel
# used in the figure; edit to match the populations retained in your object.
lineage_markers <- c(
  "Nanog","Sox2","Pou5f1","Zfp42","Klf2","Tfcp2l1","Esrrb","Prdm14",
  "Sox17","Foxa2","Hnf1b","Hnf4a","Gata4","Gata6","Dab2",
  "Pdgfra","Sox7","Lama1","Lamb1","Lamc1","Col4a1","Col4a2","Nid1",
  "Nodal","Lefty1","Cer1","Hhex","Lhx1","Otx2",
  "Eomes","T","Mixl1","Mesp1","Foxc1","Foxc2","Tbx6"
)
DotPlot(integrated_data, features = lineage_markers,
        col.min = -0.75, cols = c("white", "purple4")) +
  RotatedAxis()
ggsave(file.path(outdir, "DotPlot_lineage_markers.pdf"), width = 18, height = 8)

# =============================================================================
# 9. CELL-TYPE PSEUDOBULK
# =============================================================================
# Sum raw counts across all cells of each annotated cell type to form a
# pseudobulk profile per cell type.

pseudobulk <- AggregateExpression(integrated_data,
                                  assays = "RNA",
                                  group.by = "Update",
                                  slot = "counts",
                                  return.seurat = FALSE)
pseudobulk_mat <- pseudobulk$RNA

# =============================================================================
# 10. COMPARE PSEUDOBULK TO IN VITRO BULK RNA-seq
# =============================================================================
# Combine pseudobulk and bulk on their shared gene set, library-size normalise
# to log-CPM (edgeR), then remove the modality (bulk vs pseudobulk) batch effect
# with ComBat before comparing.

bulk_data <- read.csv("data/BulkRNA_averages.csv", row.names = 1)
bulk_data <- bulk_data[complete.cases(bulk_data), ]      # drop genes with NAs

common_genes <- intersect(rownames(pseudobulk_mat), rownames(bulk_data))
combined <- cbind(pseudobulk_mat[common_genes, ],
                  bulk_data[common_genes, ])
combined <- as.matrix(combined)

combined_cpm <- cpm(combined, log = TRUE)

batch <- c(rep("pseudobulk", ncol(pseudobulk_mat)),
           rep("bulk",       ncol(bulk_data)))
combined_combat <- ComBat(dat = combined_cpm, batch = batch)

# ---- Correlation heatmap (bulk vs pseudobulk) ------------------------------
# Restrict to the top 100 variable genes across samples, excluding cell-cycle
# and housekeeping genes, then correlate (Pearson) bulk against pseudobulk.

gene_vars <- apply(combined_combat, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:100]

cc_hk <- c("Mcm","Pcna","Top2a","Mki67","Ccn","Cdk","Cdc","Ube2c","Birc5",
           "Tpx2","Prc1","Cks","Aurk","Plk",                       # cell cycle
           "Actb","Gapdh","Tubb5","Rpl","Rps","Eef","Tuba","Tubb") # housekeeping
top_genes <- top_genes[!grepl(paste(cc_hk, collapse = "|"), top_genes,
                              ignore.case = TRUE)]

cor_mat <- cor(combined_combat[top_genes, ], method = "pearson")
pb_cols   <- seq_len(ncol(pseudobulk_mat))
bulk_cols <- (ncol(pseudobulk_mat) + 1):ncol(combined_combat)

pheatmap(cor_mat[pb_cols, bulk_cols],
         cluster_rows = TRUE, cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Bulk vs pseudobulk (top 100 variable genes)",
         fontsize = 10,
         filename = file.path(outdir, "heatmap_bulk_vs_pseudobulk.pdf"))

# ---- Global PCA ------------------------------------------------------------
# Top variable genes across all samples to visualise the global relationship.
# (Reported methods state 5000 HVG; set n_hvg accordingly.)

n_hvg <- 5000
top_hvg <- names(sort(gene_vars, decreasing = TRUE))[1:n_hvg]
top_hvg <- top_hvg[!grepl(paste(cc_hk, collapse = "|"), top_hvg,
                          ignore.case = TRUE)]

scaled <- t(scale(t(combined_combat[top_hvg, ]), center = TRUE, scale = TRUE))
pca <- prcomp(t(scaled), center = FALSE, scale. = FALSE)

pca_df <- data.frame(pca$x,
                     Sample = colnames(combined_combat),
                     Type = c(rep("scRNA pseudobulk", ncol(pseudobulk_mat)),
                              rep("Bulk RNA",          ncol(bulk_data))))

ggplot(pca_df, aes(PC1, PC2, colour = Type, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.6, size = 3) +
  theme_bw() +
  labs(title = paste0("PCA: top ", n_hvg, " variable genes"))
ggsave(file.path(outdir, "PCA_bulk_vs_pseudobulk.pdf"), width = 10, height = 8)

cat("Done. Outputs written to", outdir, "\n")

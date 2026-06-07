library(SingleCellExperiment)
library(Seurat)        # v5
library(Matrix)
library(dplyr)

set.seed(42)           

# Paths — edit to your environment ------------------------------------------
atlas_dir <- "data/ExtendedMouseAtlas/"   # where the atlas is unpacked
gast_dir  <- "data/gastruloid/"      

# =============================================================================
# DIRECT DOWNLOAD OF INDIVIDUAL FILES - LOAD THE EXTENDED MOUSE ATLAS REFERENCE
# =============================================================================
# Files are distributed individually from the Göttgens/Marioni lab server
# (landing page: https://marionilab.github.io/ExtendedMouseAtlas/).
# We need the SingleCellExperiment object and the two metadata tables.

base_url <- "https://bioinformatics.stemcells.cam.ac.uk/rlh60/Supplemental/ExtendedMouseAtlas/"

dir.create(atlas_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 3600)   # the SCE is large (~430k cells); allow a long download

for (f in c("embryo_sce.rds", "metadata_cells.csv", "metadata_genes.csv")) {
  download.file(paste0(base_url, f),
                destfile = file.path(atlas_dir, f),
                mode = "wb")
}

# Read the SingleCellExperiment (requires the SingleCellExperiment package)
sce <- readRDS(file.path(atlas_dir, "embryo_sce.rds"))

# Convert to Seurat, keeping raw counts (normalisation is done later, per dataset)
atlas <- as.Seurat(sce, counts = "counts", data = NULL)
atlas$orig.ident <- "ExtendedAtlas"
rm(sce); gc()

# =============================================================================
# Mouse gastruloid × Extended Mouse Atlas integration
# Part 1 — Data download and preprocessing
#
# Reference: Rosshandler et al. (2024) Development — Extended Mouse Atlas
# =============================================================================

# =============================================================================
# 1. BULK DOWNLOAD AND EXTRACTION of EXTENDED MOUSE ATLAS FROM TAR ARCHIVE
# =============================================================================
# The atlas is distributed as a gzipped tar archive containing a
# SingleCellExperiment (.rds) plus cell/gene metadata CSVs. We download the
# archive, extract the files we need, then read the SCE.
#
# NOTE: insert the actual download URL for the atlas archive below.

atlas_dir <- "data/ExtendedMouseAtlas/"
dir.create(atlas_dir, recursive = TRUE, showWarnings = FALSE)

download.file("https://bioinformatics.stemcells.cam.ac.uk/rlh60/Supplemental/ExtendedMouseAtlas/ExtendedMouseAtlas.tar.gz",
              destfile = "ExtendedMouseAtlas.tar.gz",
              mode = "wb")

untar("ExtendedMouseAtlas.tar.gz",
      files = c("ExtendedMouseAtlas/embryo_sce.rds",
                "ExtendedMouseAtlas/metadata_cells.csv",
                "ExtendedMouseAtlas/metadata_genes.csv"),
      exdir = atlas_dir)

sce <- readRDS(file.path(atlas_dir, "embryo_sce.rds"))
atlas <- as.Seurat(sce, counts = "counts", data = NULL)
atlas$orig.ident <- "ExtendedAtlas"
rm(sce); gc()

# =============================================================================
# 2. RENAME GENES: Ensembl IDs -> MGI symbols
# =============================================================================
# The SCE stores rows as Ensembl gene IDs. Swap to MGI symbols using the atlas
# gene metadata, making duplicate symbols unique.

gene_metadata <- read.csv(file.path(atlas_dir, "metadata_genes.csv"),
                          header = TRUE, row.names = 1)
gene_symbols  <- make.unique(gene_metadata$mgi_symbol, sep = ".")

# as.Seurat() names the imported assay "originalexp" by default
counts_matrix <- GetAssayData(atlas, assay = "originalexp", layer = "counts")
rownames(counts_matrix) <- gene_symbols

# Rebuild with gene symbols as feature names
atlas <- CreateSeuratObject(counts = counts_matrix, meta.data = atlas@meta.data)
atlas$orig.ident <- "ExtendedAtlas"
rm(counts_matrix, gene_metadata, gene_symbols); gc()

# =============================================================================
# 3. QC AND DOWNSAMPLING OF THE REFERENCE
# =============================================================================
# Recompute per-cell QC from counts and keep high-complexity cells. The
# reference uses a stricter gene-detection threshold than the query
# (>= 3000 vs >= 2000) owing to its greater per-cell sequencing depth.

counts_mat <- GetAssayData(atlas, assay = "RNA", layer = "counts")
atlas$nFeature_RNA <- colSums(counts_mat > 0)
atlas$nCount_RNA   <- colSums(counts_mat)
rm(counts_mat); gc()

atlas <- subset(atlas, subset = nFeature_RNA >= 3000)

# Downsample to a maximum of 2000 cells per annotated cell type so abundant
# populations do not dominate anchor finding. Only the reference is downsampled;
# all gastruloid cells are retained (Part 2).
cells_keep <- unlist(lapply(
  split(colnames(atlas), atlas$celltype_extended_atlas),
  function(cells) if (length(cells) > 2000) sample(cells, 2000) else cells
))
atlas <- atlas[, cells_keep]
gc()

cat("Reference cells after QC + downsampling:", ncol(atlas), "\n")
print(table(atlas$stage))

# =============================================================================
# 4. LOAD AND QC THE GASTRULOID DATASET
# =============================================================================
# Four culture conditions (CTRL, RACL, NACLB, XAL), each sampled at 48/72/96/120 h.
# Metadata already carries the gastruloid annotation (sc_merged_gast_transfer_label)
# and timepoint (time). Filter each condition to >= 2000 detected genes.
#
# NOTE: NACLB count/metadata files are named "NACL".

gast_files <- list(CTRL  = "601_ctrl_rawcounts.rds",
                   RACL  = "601_RACL_rawcounts.rds",
                   NACLB = "601_NACL_rawcounts.rds",
                   XAL   = "601_XAL_rawcounts.rds")
gast_meta  <- list(CTRL  = "601_ctrl_metadata.csv",
                   RACL  = "601_RACL_metadata.csv",
                   NACLB = "601_NACL_metadata.csv",
                   XAL   = "601_XAL_metadata.csv")

load_condition <- function(cond) {
  counts <- readRDS(file.path(gast_dir, gast_files[[cond]]))
  meta   <- read.csv(file.path(gast_dir, gast_meta[[cond]]),
                     header = TRUE, row.names = 1)
  obj <- CreateSeuratObject(counts = counts, assay = "RNA",
                            min.cells = 0, min.features = 0)
  obj <- AddMetaData(obj, metadata = meta)
  obj$orig.ident <- cond
  obj <- subset(obj, subset = nFeature_RNA >= 2000)
  cat(cond, "cells after QC:", ncol(obj), "\n")
  obj
}

gast_list <- lapply(names(gast_files), load_condition)
names(gast_list) <- names(gast_files)

# Objects ready for integration (Part 2):
#   atlas      — reference, QC-filtered and downsampled
#   gast_list  — list of four QC-filtered gastruloid condition objects

# =============================================================================
# Mouse gastruloid × Extended Mouse Atlas integration
# Part 2 — Reciprocal PCA integration
#
# Inputs from Part 1:
#   atlas      — reference, QC-filtered and downsampled
#   gast_list  — list of four QC-filtered gastruloid condition objects
# =============================================================================

library(Seurat)
library(dplyr)

set.seed(42)

# =============================================================================
# 5. MERGE REFERENCE AND GASTRULOIDS INTO ONE OBJECT
# =============================================================================
# add.cell.ids keeps cell barcodes unique across datasets. The five batches
# (atlas + four conditions) are distinguished by orig.ident.

combined <- merge(atlas,
                  y = gast_list,
                  add.cell.ids = c("Atlas", names(gast_list)),
                  merge.data = TRUE)
combined <- JoinLayers(combined)

cat("Combined object:", ncol(combined), "cells\n")
print(table(combined$orig.ident))

rm(atlas, gast_list); gc()

# =============================================================================
# 6. PER-DATASET NORMALISATION AND INTEGRATION FEATURES
# =============================================================================
# Split by batch, then normalise and find variable features within each batch.
# A shared set of integration features is selected across all five datasets.

obj_list <- SplitObject(combined, split.by = "orig.ident")

obj_list <- lapply(obj_list, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, nfeatures = 3000, verbose = FALSE)
  x
})

features <- SelectIntegrationFeatures(obj_list, nfeatures = 3000)

# =============================================================================
# 7. RECIPROCAL PCA: ANCHORS AND INTEGRATION
# =============================================================================
# Each dataset is scaled and reduced by PCA on the shared features, then
# integration anchors are found in reciprocal-PCA space (more conservative than
# CCA; appropriate for large, only-partially-overlapping datasets). Anchors are
# used to compute the batch-corrected "integrated" assay (Stuart et al., 2019).

obj_list <- lapply(obj_list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, npcs = 40, verbose = FALSE)
  x
})

anchors <- FindIntegrationAnchors(object.list = obj_list,
                                  anchor.features = features,
                                  reduction = "rpca",
                                  dims = 1:40,
                                  verbose = TRUE)
rm(obj_list, features); gc()

combined_integrated <- IntegrateData(anchorset = anchors,
                                     dims = 1:40,
                                     verbose = TRUE)
rm(anchors, combined); gc()

# =============================================================================
# 8. DIMENSIONALITY REDUCTION AND CLUSTERING ON THE INTEGRATED ASSAY
# =============================================================================
# The integrated assay is used for embedding and clustering only; differential
# expression / marker detection later switches back to the RNA assay.

DefaultAssay(combined_integrated) <- "integrated"

combined_integrated <- ScaleData(combined_integrated, verbose = FALSE)
combined_integrated <- RunPCA(combined_integrated, npcs = 50, verbose = FALSE)
combined_integrated <- RunUMAP(combined_integrated, dims = 1:40, verbose = FALSE)
combined_integrated <- FindNeighbors(combined_integrated, dims = 1:40)

# Quick sanity check — batches should be intermixed after integration
DimPlot(combined_integrated, reduction = "umap",
        group.by = "orig.ident", raster = FALSE)

# =============================================================================
# 9. SAVE THE INTEGRATED OBJECT
# =============================================================================

saveRDS(combined_integrated, "results/gastruloid_atlas_RPCA_integrated.rds")
cat("Integration complete.\n")

# =============================================================================
# Mouse gastruloid × Extended Mouse Atlas integration
# Part 3 — Publication UMAP figures
#
# Input from Part 2:
#   combined_integrated — RPCA-integrated Seurat object
#
# Colours follow a lineage-themed palette: hue encodes lineage family,
# lightness/sub-hue encodes sub-identity within a family.
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrastr)       # rasterised points keep vector files small at 600 dpi

figdir <- "figures/"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

combined_integrated <- readRDS("results/gastruloid_atlas_RPCA_integrated.rds")

# Optional: the published figures were mirrored on the x-axis purely for
# orientation. Set to TRUE to reproduce that layout.
FLIP_UMAP_X <- FALSE
if (FLIP_UMAP_X) {
  combined_integrated@reductions$umap@cell.embeddings[, 1] <-
    -combined_integrated@reductions$umap@cell.embeddings[, 1]
}

# =============================================================================
# 10. LINEAGE-THEMED PALETTE  (hue = lineage family; shade = sub-identity)
# =============================================================================

cluster_colors_themed <- c(
  # ---- PLURIPOTENT / EARLY (purples, light -> dark) ----
  "1"  = "#C8A2DB",  # Pluripotent/PGC
  "2"  = "#B87ED5",  # Epiblast
  "3"  = "#A958D2",  # Caudal epiblast
  "4"  = "#9B2FD1",  # Primitive Streak
  "5"  = "#8420B6",  # Anterior PS
  
  # ---- AXIAL MESODERM (oranges -> reds) ----
  "6"  = "#DEB190",  # Nascent mesoderm
  "7"  = "#A84423",  # Node
  "8"  = "#CB624D",  # Notochord
  "9"  = "#DCA37A",  # NMPs
  "10" = "#DB9664",  # NMPs meso-biased
  "11" = "#DA884D",  # Caudal mesoderm
  "12" = "#DB7A35",  # Paraxial mesoderm
  "13" = "#D86D20",  # Presomitic mesoderm
  "14" = "#C86219",  # Somitic mesoderm
  "15" = "#B75712",  # Posterior somitic
  "16" = "#A54C0D",  # Anterior somitic tissues
  "17" = "#CC323F",  # Dermomyotome
  "18" = "#8A2828",  # Sclerotome
  "19" = "#62341C",  # Endotome
  
  # ---- LPM / IM / CARDIAC (pinks -> magentas, FHF warm / SHF cool) ----
  "20" = "#E2AAC6",  # Lateral plate mesoderm
  "21" = "#DC80AE",  # Intermediate mesoderm
  "22" = "#DA5296",  # Kidney primordium
  "23" = "#D9227E",  # Limb mesoderm
  "24" = "#B71466",  # Cranial mesoderm
  "25" = "#AE6980",  # Mesenchyme
  "26" = "#DD5FB3",  # Cardiopharyngeal mesoderm
  "27" = "#D8269D",  # Cardiopharyngeal progenitors
  "28" = "#B02576",  # Anterior cardiopharyngeal progenitors
  "29" = "#D368C1",  # Cardiopharyngeal progenitors FHF
  "30" = "#C525AF",  # Cardiomyocytes FHF 1
  "31" = "#911186",  # Cardiomyocytes FHF 2
  "32" = "#C168D3",  # Cardiopharyngeal progenitors SHF
  "33" = "#A525C5",  # Cardiomyocytes SHF 1
  "34" = "#711191",  # Cardiomyocytes SHF 2
  "35" = "#BF3F74",  # Pharyngeal mesoderm
  "36" = "#5F3181",  # Epicardium
  
  # ---- ENDODERM (greens) ----
  "37" = "#97D7A7",  # Gut tube
  "38" = "#6ACF83",  # Foregut
  "39" = "#3AC95E",  # Midgut
  "40" = "#25A946",  # Hindgut
  "41" = "#168231",  # Pharyngeal endoderm
  "42" = "#288A69",  # Thyroid primordium
  "43" = "#50A63F",  # Visceral endoderm
  "44" = "#537C35",  # Parietal endoderm
  "45" = "#9DBF58",  # ExE endoderm
  
  # ---- ECTODERM / NEURAL (browns + blues) ----
  "46" = "#C09E7C",  # Ectoderm
  "47" = "#B78451",  # Embryonic ectoderm
  "48" = "#9E6B37",  # Surface ectoderm
  "49" = "#7E5124",  # Non-neural ectoderm
  "50" = "#A58C59",  # ExE ectoderm
  "51" = "#D1BA47",  # Placodal ectoderm
  "52" = "#79D2D2",  # Otic placode
  "53" = "#28BDBD",  # Otic vesicle
  "54" = "#7DA5CD",  # Neural tube
  "55" = "#4780D1",  # Spinal cord progenitors
  "56" = "#164E9C",  # Dorsal spinal cord progenitors
  "57" = "#74BED6",  # Hindbrain floor plate
  "58" = "#3EACD1",  # Ventral hindbrain progenitors
  "59" = "#218EB2",  # Hindbrain neural progenitors
  "60" = "#136985",  # Dorsal hindbrain progenitors
  "61" = "#395FAC",  # Midbrain/Hindbrain boundary
  "62" = "#6073D1",  # Midbrain progenitors
  "63" = "#2841BD",  # Dorsal midbrain neurons
  "64" = "#9F8CD8",  # Early dorsal forebrain progenitors
  "65" = "#7050D0",  # Dorsal forebrain progenitors
  "66" = "#4924B8",  # Late dorsal forebrain progenitors
  "67" = "#2F1385",  # Ventral forebrain progenitors
  "68" = "#3BC1DC",  # Optic vesicle
  
  # ---- NEURAL CREST (teals) ----
  "69" = "#4DCBAB",  # Migratory neural crest
  "70" = "#25B28F",  # Branchial arch neural crest
  "71" = "#138569",  # Frontonasal mesenchyme
  
  # ---- VASCULAR / ENDOTHELIAL (pink-reds) ----
  "72" = "#D49B9F",  # Haematoendothelial progenitors
  "73" = "#CB6E75",  # Embryo proper endothelium
  "74" = "#C63D48",  # Venous endothelium
  "75" = "#CB4D6C",  # YS endothelium
  "76" = "#B13368",  # Allantois endothelium
  "77" = "#82161F",  # Endocardium
  
  # ---- BLOOD (deep reds -> burgundy) ----
  "78" = "#CB4D4D",  # Blood progenitors
  "79" = "#C23232",  # MEP
  "80" = "#AA2626",  # EMP
  "81" = "#901C1C",  # Megakaryocyte progenitors
  "82" = "#751313",  # Erythroid
  "83" = "#590C0C",  # Chorioallantoic-derived erythroid progenitors
  
  # ---- EXTRAEMBRYONIC (yellows / olives) ----
  "84" = "#D2CA79",  # Allantois
  "85" = "#CEC23D",  # YS mesothelium
  "86" = "#ADA11E",  # YS mesothelium-derived endothelial progenitors
  "87" = "#B89C46",  # Amniotic ectoderm
  
  # ---- FALLBACK (greys) ----
  "88" = "#BDBDBD",  # NA
  "89" = "#E0E0E0"   # Unassigned
)


# Cell-type label -> palette number. Synonyms and alternate dashes map to the
# same number; lineage sub-types are collapsed where noted. Unmapped labels
# fall through to grey (#89).

label_to_number <- c(
  # ---- PLURIPOTENT / EARLY ----
  "PGC"                                           = "1",
  "Pluripotent/PGC"                               = "1",
  "Epiblast"                                      = "2",
  "Caudal epiblast"                               = "3",
  "Primitive Streak"                              = "4",
  "Primitive Streak (PS)"                         = "4",
  "Anterior Primitive Streak"                     = "5",
  "Anterior PS"                                   = "5",
  
  # ---- AXIAL MESODERM ----
  "Nascent mesoderm"                              = "6",
  "Node"                                          = "7",
  "Notochord"                                     = "8",
  "NMPs"                                          = "9",
  "NMPs/Mesoderm-biased"                          = "10",
  "NMPs meso-biased"                              = "10",
  "Caudal mesoderm"                               = "11",
  "Paraxial mesoderm"                             = "12",
  "Presomitic mesoderm"                           = "13",
  "Somitic mesoderm"                              = "14",
  "Posterior somitic tissues"                     = "15",
  "Posterior somitic"                             = "15",
  "Anterior somitic tissues"                      = "16",
  "Dermomyotome"                                  = "17",
  "Sclerotome"                                    = "18",
  "Endotome"                                      = "19",
  
  # ---- LPM / IM / CARDIAC ----
  "Lateral plate mesoderm"                        = "20",
  "Intermediate mesoderm"                         = "21",
  "Kidney primordium"                             = "22",
  "Limb mesoderm"                                 = "23",
  "Forelimb"                                      = "23",   # -> Limb mesoderm
  "Cranial mesoderm"                              = "24",
  "Mesenchyme"                                    = "25",
  "Embryo proper mesothelium"                     = "25",   # -> Mesenchyme
  "Cardiopharyngeal mesoderm"                     = "26",
  "Cardiopharyngeal progenitors"                  = "27",
  "Anterior cardiopharyngeal progenitors"         = "28",
  "Cardiopharyngeal progenitors FHF"              = "29",
  "Cardiomyocytes FHF 1"                          = "30",
  "Cardiomyocytes FHF 2"                          = "31",
  "Cardiopharyngeal progenitors SHF"              = "32",
  "Cardiomyocytes SHF 1"                          = "33",
  "Cardiomyocytes SHF 2"                          = "34",
  "Pharyngeal mesoderm"                           = "35",
  "Epicardium"                                    = "36",
  
  # ---- ENDODERM ----
  "Gut tube"                                      = "37",
  "Foregut"                                       = "38",
  "Midgut"                                        = "39",
  "Hindgut"                                       = "40",
  "Pharyngeal endoderm"                           = "41",
  "Thyroid primordium"                            = "42",
  "Visceral endoderm"                             = "43",
  "Parietal endoderm"                             = "44",
  "ExE endoderm"                                  = "45",
  
  # ---- ECTODERM / NEURAL ----
  "Ectoderm"                                      = "46",
  "Embryonic ectoderm"                            = "47",
  "Surface ectoderm"                              = "48",
  "Epidermis"                                     = "48",   # -> Surface ectoderm
  "Limb ectoderm"                                 = "48",   # -> Surface ectoderm
  "Non-neural ectoderm"                           = "49",
  "Non\u2212neural ectoderm"                      = "49",   # alt unicode dash
  "ExE ectoderm"                                  = "50",
  "Placodal ectoderm"                             = "51",
  "Otic placode"                                  = "52",
  "Otic vesicle"                                  = "53",
  "Otic neural progenitors"                       = "53",   # -> Otic vesicle
  "Neural tube"                                   = "54",
  "Spinal cord progenitors"                       = "55",
  "Dorsal spinal cord progenitors"                = "56",
  "Hindbrain floor plate"                         = "57",
  "Ventral hindbrain progenitors"                 = "58",
  "Hindbrain neural progenitors"                  = "59",
  "Dorsal hindbrain progenitors"                  = "60",
  "Midbrain/Hindbrain boundary"                   = "61",
  "Midbrain progenitors"                          = "62",
  "Dorsal midbrain neurons"                       = "63",
  "Early dorsal forebrain progenitors"            = "64",
  "Dorsal forebrain progenitors"                  = "65",
  "Late dorsal forebrain progenitors"             = "66",
  "Ventral forebrain progenitors"                 = "67",
  "Optic vesicle"                                 = "68",
  
  # ---- NEURAL CREST ----
  "Migratory neural crest"                        = "69",
  "Branchial arch neural crest"                   = "70",
  "Frontonasal mesenchyme"                        = "71",
  
  # ---- VASCULAR / ENDOTHELIAL ----
  "Haematoendothelial progenitors"                = "72",
  "Embryo proper endothelium"                     = "73",
  "Venous endothelium"                            = "74",
  "YS endothelium"                                = "75",
  "Allantois endothelium"                         = "76",
  "Endocardium"                                   = "77",
  
  # ---- BLOOD ----
  "Blood progenitors"                             = "78",
  "MEP"                                           = "79",
  "EMP"                                           = "80",
  "Megakaryocyte progenitors"                     = "81",
  "Erythroid"                                     = "82",
  "Chorioallantoic-derived erythroid progenitors" = "83",
  "Chorioallantoic\u2212derived erythroid progenitors" = "83",  # alt unicode
  
  # ---- EXTRAEMBRYONIC ----
  "Allantois"                                     = "84",
  "YS mesothelium"                                = "85",
  "YS mesothelium-derived endothelial progenitors" = "86",
  "YS mesothelium\u2212derived endothelial progenitors" = "86",
  "Amniotic ectoderm"                             = "87"
)

# Helpers: number -> canonical label, and "N — Label" legend strings
number_to_label  <- tapply(names(label_to_number), label_to_number, `[`, 1)
make_legend_label <- function(num) paste0(num, " \u2014 ", number_to_label[as.character(num)])

# =============================================================================
# 11. BUILD A PLOTTING DATA FRAME + CHECK FOR UNMAPPED LABELS
# =============================================================================

umap_full <- as.data.frame(combined_integrated@reductions$umap@cell.embeddings)
colnames(umap_full)[1:2] <- c("umap_1", "umap_2")
umap_full$orig.ident <- combined_integrated$orig.ident
umap_full$atlas_label <- combined_integrated$celltype_extended_atlas
umap_full$gast_label  <- combined_integrated$sc_merged_gast_transfer_label
umap_full$stage       <- combined_integrated$stage
umap_full$time        <- combined_integrated$time

atlas_bg   <- umap_full[umap_full$orig.ident == "ExtendedAtlas", ]
conditions <- c("CTRL", "RACL", "NACLB", "XAL")

# Any label not in the mapping will be drawn grey — check before plotting
report_unmapped <- function(labels, what) {
  miss <- setdiff(unique(na.omit(labels)), names(label_to_number))
  cat("Unmapped", what, "labels (drawn grey):",
      if (length(miss)) paste(miss, collapse = ", ") else "none", "\n")
}
report_unmapped(umap_full$atlas_label, "atlas")
report_unmapped(umap_full$gast_label,  "gastruloid")

# =============================================================================
# 12. PLOT STYLE + REUSABLE CELL-TYPE UMAP FUNCTION
# =============================================================================

PT_SIZE    <- 0.4    # foreground points
BG_SIZE    <- 0.2    # grey background points
RASTER_DPI <- 600
FIG_W <- 6; FIG_H <- 6   # inches (adjust per panel as needed)

base_theme <- theme_classic(base_size = 9) +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 10, face = "bold"),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

# Colour a set of cells by a cell-type label column, over an optional grey
# background. Works for the atlas (atlas_label) or gastruloids (gast_label).
plot_celltype_umap <- function(fg, bg = NULL, label_col,
                               title = NULL, show_legend = TRUE) {
  fg <- fg[!is.na(fg[[label_col]]), ]
  fg$number <- label_to_number[as.character(fg[[label_col]])]
  fg$number[is.na(fg$number)] <- "89"
  present <- sort(unique(as.integer(fg$number)))
  fg$number <- factor(fg$number, levels = as.character(present))
  set.seed(42); fg <- fg[sample(nrow(fg)), ]   # shuffle so no type sits on top
  
  cols <- cluster_colors_themed[as.character(present)]
  labs <- make_legend_label(as.character(present)); names(labs) <- names(cols)
  
  p <- ggplot()
  if (!is.null(bg))
    p <- p + geom_point_rast(data = bg, aes(umap_1, umap_2),
                             colour = "grey85", size = BG_SIZE, alpha = 0.4,
                             raster.dpi = RASTER_DPI)
  p <- p +
    geom_point_rast(data = fg, aes(umap_1, umap_2, colour = number),
                    size = PT_SIZE, alpha = 0.85, raster.dpi = RASTER_DPI) +
    scale_colour_manual(values = cols, labels = labs, name = NULL, drop = FALSE) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    base_theme
  
  if (show_legend)
    p + guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 2))
  else
    p + guides(colour = "none")
}

# =============================================================================
# 13. FIGURE: ATLAS REFERENCE COLOURED BY CELL TYPE
# =============================================================================

p_atlas <- plot_celltype_umap(
  fg = atlas_bg, bg = NULL, label_col = "atlas_label",
  title = "Extended Mouse Atlas reference", show_legend = TRUE
)
ggsave(file.path(figdir, "atlas_reference_celltypes.pdf"),
       p_atlas, width = 9, height = 6)

# =============================================================================
# 14. FIGURE: GASTRULOID CELL TYPES ON GREY ATLAS BACKGROUND
# =============================================================================

# Per condition (shared legend, 2x2 grid)
cond_panels <- lapply(conditions, function(cond) {
  plot_celltype_umap(
    fg = umap_full[umap_full$orig.ident == cond, ],
    bg = atlas_bg, label_col = "gast_label",
    title = cond, show_legend = FALSE
  )
})
names(cond_panels) <- conditions

cond_grid <- (cond_panels[["CTRL"]]  + cond_panels[["RACL"]]) /
  (cond_panels[["NACLB"]] + cond_panels[["XAL"]])
ggsave(file.path(figdir, "gastruloid_celltypes_by_condition.pdf"),
       cond_grid, width = 10, height = 10)

# All gastruloids pooled (with legend)
p_gast_all <- plot_celltype_umap(
  fg = umap_full[umap_full$orig.ident %in% conditions, ],
  bg = atlas_bg, label_col = "gast_label",
  title = "All gastruloids on atlas", show_legend = TRUE
)
ggsave(file.path(figdir, "gastruloid_celltypes_all.pdf"),
       p_gast_all, width = 9, height = 6)

# =============================================================================
# 15. FIGURE: EMBRYONIC STAGES (atlas coloured by stage, gastruloids grey)
# =============================================================================

stage_levels <- c("E6.5","E6.75","E7.0","E7.25","E7.5","E7.75","E8.0","E8.25",
                  "E8.5","E8.75","E9.0","E9.25","E9.5","Mixed gastrulation")
fg_stage <- atlas_bg[!is.na(atlas_bg$stage), ]
fg_stage$stage <- factor(fg_stage$stage, levels = stage_levels)
set.seed(42); fg_stage <- fg_stage[sample(nrow(fg_stage)), ]

p_stage <- ggplot() +
  geom_point_rast(data = umap_full[umap_full$orig.ident %in% conditions, ],
                  aes(umap_1, umap_2), colour = "grey85",
                  size = BG_SIZE, alpha = 0.4, raster.dpi = RASTER_DPI) +
  geom_point_rast(data = fg_stage, aes(umap_1, umap_2, colour = stage),
                  size = PT_SIZE, alpha = 0.85, raster.dpi = RASTER_DPI) +
  scale_colour_viridis_d(option = "plasma", na.value = "grey85", name = "Stage") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(title = "Embryonic stage (gastruloids in grey)",
       x = "UMAP 1", y = "UMAP 2") +
  base_theme
ggsave(file.path(figdir, "atlas_embryonic_stages.pdf"),
       p_stage, width = 8, height = 6)

# =============================================================================
# 16. FIGURE: GASTRULOID TIMEPOINTS (per condition, on grey atlas)
# =============================================================================

timepoint_colours <- c("48h"="#FDE725","72h"="#5DC863","96h"="#21908C","120h"="#440154")
tp_levels <- names(timepoint_colours)

plot_timepoints <- function(cond) {
  fg <- umap_full[umap_full$orig.ident == cond & !is.na(umap_full$time), ]
  fg$time <- factor(fg$time, levels = tp_levels)
  fg <- fg[order(fg$time), ]   # earliest underneath, latest on top
  ggplot() +
    geom_point_rast(data = atlas_bg, aes(umap_1, umap_2), colour = "grey85",
                    size = BG_SIZE, alpha = 0.4, raster.dpi = RASTER_DPI) +
    geom_point_rast(data = fg, aes(umap_1, umap_2, colour = time),
                    size = PT_SIZE, alpha = 0.85, raster.dpi = RASTER_DPI) +
    scale_colour_manual(values = timepoint_colours, name = "Timepoint", drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(title = cond, x = "UMAP 1", y = "UMAP 2") +
    base_theme
}

tp_grid <- (plot_timepoints("CTRL")  + plot_timepoints("RACL")) /
  (plot_timepoints("NACLB") + plot_timepoints("XAL")) +
  plot_layout(guides = "collect") & theme(legend.position = "right")
ggsave(file.path(figdir, "gastruloid_timepoints_by_condition.pdf"),
       tp_grid, width = 11, height = 11)

# =============================================================================
# 17. FIGURE: ALL GASTRULOIDS COLOURED BY CONDITION (on grey atlas)
# =============================================================================

condition_colours <- c("CTRL"="#E41A1C","RACL"="#377EB8",
                       "NACLB"="#4DAF4A","XAL"="#984EA3")
fg_cond <- umap_full[umap_full$orig.ident %in% conditions, ]
fg_cond$orig.ident <- factor(fg_cond$orig.ident, levels = conditions)
set.seed(42); fg_cond <- fg_cond[sample(nrow(fg_cond)), ]

p_cond <- ggplot() +
  geom_point_rast(data = atlas_bg, aes(umap_1, umap_2), colour = "grey85",
                  size = BG_SIZE, alpha = 0.4, raster.dpi = RASTER_DPI) +
  geom_point_rast(data = fg_cond, aes(umap_1, umap_2, colour = orig.ident),
                  size = PT_SIZE, alpha = 0.85, raster.dpi = RASTER_DPI) +
  scale_colour_manual(values = condition_colours, name = "Condition") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(title = "Gastruloids by condition", x = "UMAP 1", y = "UMAP 2") +
  base_theme
ggsave(file.path(figdir, "gastruloid_by_condition.pdf"),
       p_cond, width = 8, height = 6)

cat("\nFigures written to", figdir, "\n")
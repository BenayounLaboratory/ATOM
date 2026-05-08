###########################################################################################
# Aging transcriptomics of mouse macrophages 
# Process publicly available single cell feature matrix and pseudobulk the single cell data
###########################################################################################

setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/new_data_08_01_2025/GSE205395_skeletal")

################################################################################################
######################## 1. Load in necessary data and packages

library(Seurat)
library(hdf5r)
library(DelayedArray)
library(muscat)
library(glmGamPoi)


# 1. Read and label each sample ---------------------------------------

# Young Injured Rep 1
yi1 <- Read10X(data.dir = "young_injured_rep1_feature_bc_matrix")
yi1_obj <- CreateSeuratObject(yi1)
yi1_obj$sample <- "young_injured_rep1"
yi1_obj$age <- "young"
yi1_obj$condition <- "injured"
yi1_obj$replicate <- "1"

# Young Injured Rep 2
yi2 <- Read10X(data.dir = "young_injured_rep2_feature_bc_matrix")
yi2_obj <- CreateSeuratObject(yi2)
yi2_obj$sample <- "young_injured_rep2"
yi2_obj$age <- "young"
yi2_obj$condition <- "injured"
yi2_obj$replicate <- "2"

# Aged Injured Rep 1
ai1 <- Read10X(data.dir = "aged_injured_rep1_feature_bc_matrix")
ai1_obj <- CreateSeuratObject(ai1)
ai1_obj$sample <- "aged_injured_rep1"
ai1_obj$age <- "aged"
ai1_obj$condition <- "injured"
ai1_obj$replicate <- "1"

# Aged Injured Rep 2
ai2 <- Read10X(data.dir = "aged_injured_rep2_feature_bc_matrix")
ai2_obj <- CreateSeuratObject(ai2)
ai2_obj$sample <- "aged_injured_rep2"
ai2_obj$age <- "aged"
ai2_obj$condition <- "injured"
ai2_obj$replicate <- "2"


# 2. Merge all into one Seurat object ---------------------------------

combined <- merge(
  yi1_obj,
  y = list(yi2_obj, ai1_obj, ai2_obj),
  add.cell.ids = c("yi1", "yi2", "ai1", "ai2"),
  project = "GSE205395_injured",
  merge.data = TRUE
)

# 3. Join Layers ------------------------------------------------------

combined <- JoinLayers(combined)

# Verify cell counts
expected_total <- ncol(yi1_obj) + ncol(yi2_obj) + ncol(ai1_obj) + ncol(ai2_obj)
print(paste("Expected:", expected_total, "| Got:", ncol(combined)))


# 4. QC and filtering -------------------------------------------------

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")

VlnPlot(combined, group.by = "sample", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

combined <- subset(
  combined,
  subset = percent.mt < 25 &
    nFeature_RNA > 250 &
    nCount_RNA < 50000
)

VlnPlot(combined, group.by = "sample", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))


# 5. Normalize with SCTransform ---------------------------------------

combined <- SCTransform(
  combined,
  vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  method = "glmGamPoi",
  conserve.memory = TRUE,
  verbose = TRUE
)

combined

# 6. Select for cells with macrophage markers -------------------------

combined <- RunPCA(combined)

FeaturePlot(combined, features = c("Itgam", "Adgre1"))

combined_macrophage_filter <- subset(
  combined,
  Itgam > 0 &
    Adgre1 > 0
)

combined_macrophage_filter


# 7. Pseudobulking ----------------------------------------------------

pseudobulk_combined <- AggregateExpression(
  combined_macrophage_filter,
  group.by = "sample",
  assay = "RNA",
  slot = "counts"
)

count_data <- as.data.frame(as.matrix(pseudobulk_combined$RNA))

colnames(count_data) <- c(
  "AgedInjured1",
  "AgedInjured2",
  "YoungInjured1",
  "YoungInjured2"
)

filename <- paste0(Sys.Date(), "_GSE205395_SkM_macrophages_injured_counts.txt")
write.table(count_data, file = filename, sep = "\t", quote = FALSE, col.names = NA)

# Session info
sink(file = paste0(Sys.Date(), "_GSE205395_SkM_scRNAseq_pseudoBulk_session_info.txt"))
sessionInfo()
sink()
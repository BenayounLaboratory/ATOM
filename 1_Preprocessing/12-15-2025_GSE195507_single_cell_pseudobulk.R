###########################################################################################
# Aging transcriptomics of mouse macrophages 
# Process publicly available single cell feature matrix and pseudobulk the single cell data
###########################################################################################

setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/new_data_08_01_2025/PRJNA800823_SkM/GSE195507_RAW/")

################################################################################################
######################## 1. Load in necessary data and packages

library(Seurat)
library(hdf5r)
library(DelayedArray)
library(muscat)
library(glmGamPoi)


# Set paths to H5 files
Macrophage_1_Old_h5 <- "GSM5839574_Macrophage_1_Old_filtered_feature_bc_matrix.h5"
Macrophage_1_Young_h5 <- "GSM5839575_Macrophage_1_Young_filtered_feature_bc_matrix.h5"

Macrophage_2_Old_h5 <- "GSM5839576_Macrophage_2_Old_filtered_feature_bc_matrix.h5"
Macrophage_2_Young_h5 <- "GSM5839577_Macrophage_2_Young_filtered_feature_bc_matrix.h5"

Macrophage_3_Old_h5 <- "GSM5839578_Macrophage_3_Old_filtered_feature_bc_matrix.h5"
Macrophage_3_Young_h5 <- "GSM5839579_Macrophage_3_Young_filtered_feature_bc_matrix.h5"

# 1. Read and label each sample ---------------------------------------

# Macrophage 1 - Old
m1_old <- Read10X_h5("GSM5839574_Macrophage_1_Old_filtered_feature_bc_matrix.h5")
m1_old_obj <- CreateSeuratObject(m1_old)
m1_old_obj$sample <- "Macrophage_1_Old"
m1_old_obj$age <- "old"
m1_old_obj$replicate <- "1"

# Macrophage 1 - Young
m1_young <- Read10X_h5("GSM5839575_Macrophage_1_Young_filtered_feature_bc_matrix.h5")
m1_young_obj <- CreateSeuratObject(m1_young)
m1_young_obj$sample <- "Macrophage_1_Young"
m1_young_obj$age <- "young"
m1_young_obj$replicate <- "1"

# Macrophage 2 - Old
m2_old <- Read10X_h5("GSM5839576_Macrophage_2_Old_filtered_feature_bc_matrix.h5")
m2_old_obj <- CreateSeuratObject(m2_old)
m2_old_obj$sample <- "Macrophage_2_Old"
m2_old_obj$age <- "old"
m2_old_obj$replicate <- "2"

# Macrophage 2 - Young
m2_young <- Read10X_h5("GSM5839577_Macrophage_2_Young_filtered_feature_bc_matrix.h5")
m2_young_obj <- CreateSeuratObject(m2_young)
m2_young_obj$sample <- "Macrophage_2_Young"
m2_young_obj$age <- "young"
m2_young_obj$replicate <- "2"

# Macrophage 3 - Old
m3_old <- Read10X_h5("GSM5839578_Macrophage_3_Old_filtered_feature_bc_matrix.h5")
m3_old_obj <- CreateSeuratObject(m3_old)
m3_old_obj$sample <- "Macrophage_3_Old"
m3_old_obj$age <- "old"
m3_old_obj$replicate <- "3"

# Macrophage 3 - Young
m3_young <- Read10X_h5("GSM5839579_Macrophage_3_Young_filtered_feature_bc_matrix.h5")
m3_young_obj <- CreateSeuratObject(m3_young)
m3_young_obj$sample <- "Macrophage_3_Young"
m3_young_obj$age <- "young"
m3_young_obj$replicate <- "3"


# 2. Merge all into one Seurat object ---------------------------------

combined <- merge(
  m1_old_obj,
  y = list(m1_young_obj, m2_old_obj, m2_young_obj, m3_old_obj, m3_young_obj),
  add.cell.ids = c(
    "m1_old", "m1_young", "m2_old", "m2_young", "m3_old", "m3_young"
  ),
  project = "Macrophage_Aging",
  merge.data = TRUE
)

names(combined[["RNA"]]@layers)
combined
# An object of class Seurat 
# 31053 features across 21914 samples within 1 assay 
# Active assay: RNA (31053 features, 0 variable features)
# 6 layers present: counts.1, counts.2, counts.3, counts.4, counts.5, counts.6

# 3. Join Layers ---------------------------------

combined <- JoinLayers(combined)

# count to make sure the number of cell are still the same 
original_cells_m1_old <- colnames(m1_old_obj)


# Get the total number of cells in the final object
total_cells_combined <- ncol(combined)
print(total_cells_combined)

# Get the number of cells for each object
n_m1_old <- ncol(m1_old_obj)
n_m1_young <- ncol(m1_young_obj)
n_m2_old <- ncol(m2_old_obj)
n_m2_young <- ncol(m2_young_obj)
n_m3_old <- ncol(m3_old_obj)
n_m3_young <- ncol(m3_young_obj)

# Verify that the total count is correct
expected_total <- n_m1_old + n_m1_young + n_m2_old + n_m2_young + n_m3_old + n_m3_young
print(expected_total)



# 4. QC and filtering  ---------------------------------

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")

VlnPlot(combined, group.by = "sample", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

# Add macrophage marker expression filtering
combined <- subset(
  combined,
  subset = percent.mt < 25 & 
    nFeature_RNA > 250 & 
    nCount_RNA < 50000 
)

VlnPlot(combined, group.by = "sample", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

# 5. Normalize with SCTransform  ---------------------------------


combined <- SCTransform(
  combined,
  vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  method = "glmGamPoi",             
  conserve.memory = TRUE,           # Reduce RAM usage
  verbose = TRUE
)

combined
#48972 features across 20286 samples within 2 assays 
#Active assay: SCT (17919 features, 3000 variable features)
#3 layers present: counts, data, scale.data
#1 other assay present: RNA
# 

# 6. Select for cells with macrophage markers   ---------------------------------

# Run PCA 
 combined <- RunPCA(combined)
 
# plot the feature expression
 FeaturePlot(combined, features = c("Itgam", "Adgre1"))
 
 combined_macrophage_filter <- subset(
   combined,
     Itgam > 0 &    # Filter for Itgam expression (Cd11b)
     Adgre1 > 0     # Filter for Adgre1 expression (F4/80)
 )

 combined_macrophage_filter

 # 7. Pseudobulking   ---------------------------------
setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/new_data_08_01_2025/PRJNA800823_SkM/")

pseudobulk_combined <- AggregateExpression(combined_macrophage_filter, group.by = "sample", assay = "RNA", slot = "counts")

count_data <- as.data.frame(as.matrix(pseudobulk_combined$RNA))

# rename the samples
colnames(count_data) <- c(
  "Macrophage1_24m_male",
  "Macrophage1_3m_male",
  "Macrophage2_24m_male",
  "Macrophage2_3m_male",
  "Macrophage3_24m_male",
  "Macrophage3_3m_male"
)

# Construct the filename
filename <- paste0(Sys.Date(), "_PRJNA800823_Skeletal_Muscle_macrophages_aging_counts.txt")

# Save the count_data to a tab-delimited text file
write.table(count_data, file = filename, sep = "\t", quote = FALSE, col.names = NA)

#######################
setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/new_data_08_01_2025/PRJNA800823_SkM/")
sink(file = paste(Sys.Date(),"_PRJNA800823_SkM_scRNAseq_pseudoBulk_analysis_session_info_.txt", sep =""))
sessionInfo()
sink()











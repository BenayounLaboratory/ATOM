####################################################################################
# Aging transcriptomics of mouse macrophages 
# Mitch multivariate gene set enrichment analysis - common pathway analysis
####################################################################################

set.seed(1234) # set seed for reproducibility

################################################################################
# 1. Read in the data
################################################################################


library("mitch")
library("msigdbr")
library(org.Mm.eg.db)  
library(ggplot2)
library(reshape2)  
library(stringr)
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)



#put into list
setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/2026-01-28_continuous_age_output")
deseq_files <- list.files(pattern = "all_genes_statistics\\.txt$")

my_deseq_list <- list()
# Loop through each file in the directory
for (file in deseq_files) {
  clean_name <- str_replace(file, "^\\d{4}-\\d{2}-\\d{2}_", "")
  dataset_name <- str_replace(clean_name, "_AGE_DIM_all_genes_statistics\\.txt$", "")
  dataset_name <- str_remove(dataset_name, "_$")  # Remove trailing underscore
  dataset <- read.csv(file, header = TRUE, row.names = NULL, sep = "\t")
  my_deseq_list[[dataset_name]] <- dataset
}

# Filter for datasets that passed QC
pass_qc <- readLines("2026-01-28_datasets_pass_qc.txt")
my_deseq_list <- my_deseq_list[names(my_deseq_list) %in% pass_qc]


padj_cutoff <- 0.05



################################################################################
# 2. Reorder based on niche
################################################################################

setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/MitchOutput")

niche_map <- list(
  Microglia = c(
    "GSE98401_Microglia",  
    "GSE131869_Microglia_M",
    "GSE131869_Microglia_F",
    "GSE137028_Microglia",
    "GSE156762_Microglia",
    "GSE267529_Microglia_M",
    "GSE267529_Microglia_F",
    "PRJNA524906_Microglia"
  ),
  
  Alveolar = c(
    "GSE134397_Alveolar",
    "GSE134397_Alveolar_CTL",
    "GSE145295_Alveolar",
    "GSE190689_Alveolar"
  ),
  
  Spleen = c(
    "GSE93202_Spleen",
    "GSE199879_Spleen_Red_Pulp"
  ),
  
  Adipose = c(
    "GSE154832_eWAT",
    "GSE93202_VAT"
  ),
  
  Nerve = c(
    "GSE132882_Nerve"
  ),
  
  Peritoneal = c(
    "GSE128830_Peritoneal"
  ),
  
  Callus = c(
    "PRJNA682234_Callus",
    "PRJNA816431_Callus"
  ),
  
  Skin = c(
    "GSE199763_SkinWound"
  ),
  
  Skeletal = c(
    "GSE205395_Skeletal",
    "GSE142580_SkM",
    "PRJNA800823_SkM"
  ),
  
  BMDM = c(
    "BMDM_NIA_F",
    "PRJNA1173774_BMDM"
  )
)



sex_map <- c(
  # Males
  "GSE93202_Spleen" = "Male",
  "GSE93202_VAT" = "Male",
  "GSE98401_Microglia" = "Male",
  "GSE134397_Alveolar_CTL" = "Male",
  "GSE131869_Microglia_M" = "Male",
  "GSE128830_Peritoneal" = "Male",
  "GSE137028_Microglia" = "Male",
  "GSE154832_eWAT" = "Male",
  "PRJNA682234_Callus" = "Male",
  "GSE267529_Microglia_M" = "Male",
  "GSE145295_Alveolar" = "Male",
  "GSE142580_SkM" = "Male",
  "GSE190689_Alveolar" = "Male",
  "PRJNA800823_SkM" = "Male",
  "PRJNA1173774_BMDM" = "Male",
  "PRJNA524906_Microglia" = "Male",
  "PRJNA816431_Callus" = "Male",
  
  # Females
  "GSE199763_SkinWound" = "Female",
  "GSE199879_Spleen_Red_Pulp" = "Female",
  "GSE156762_Microglia" = "Female",
  "GSE267529_Microglia_F" = "Female",
  "GSE132882_Nerve" = "Female",
  "GSE131869_Microglia_F" = "Female",
  "BMDM_NIA_F" = "Female",
  "GSE134397_Alveolar" = "Female",
  "GSE205395_Skeletal" = "Female"
)

# temporary aliases to avoid name duplication with Mitch
dataset.aliases <- list(
  # Male datasets
  "GSE93202_Spleen"           = "Spleen1",
  "GSE93202_VAT"              = "Adipose1",
  "GSE98401_Microglia"        = "Microglia1",
  "GSE134397_Alveolar_CTL"    = "Alveolar1",
  "GSE131869_Microglia_M"     = "Microglia2",
  "GSE128830_Peritoneal"      = "Peritoneal1",
  "GSE137028_Microglia"       = "Microglia3",
  "GSE154832_eWAT"            = "Adipose2",
  "PRJNA682234_Callus"        = "Callus1",
  "GSE267529_Microglia_M"     = "Microglia4",
  "GSE145295_Alveolar"        = "Alveolar2",
  "GSE142580_SkM"             = "SkM1",
  "GSE190689_Alveolar"        = "Alveolar3",
  "PRJNA800823_SkM"           = "SkM2",
  "PRJNA1173774_BMDM"         = "BMDM1",
  "PRJNA524906_Microglia"     = "Microglia5",
  "PRJNA816431_Callus"        = "Callus2",
  
  # Female datasets
  "GSE199763_SkinWound"       = "Skin1",
  "GSE199879_Spleen_Red_Pulp" = "Spleen2",
  "GSE156762_Microglia"       = "Microglia6",
  "GSE267529_Microglia_F"     = "Microglia7",
  "GSE131869_Microglia_F"     = "Microglia8",
  "BMDM_NIA_F"                = "BMDM2",
  "GSE134397_Alveolar"        = "Alveolar4"
)


alias_to_dataset <- c(
  # Male
  "Spleen1"      = "GSE93202_Spleen",
  "Adipose1"     = "GSE93202_VAT",
  "Microglia1"   = "GSE98401_Microglia",
  "Alveolar1"    = "GSE134397_Alveolar_CTL",
  "Microglia2"   = "GSE131869_Microglia_M",
  "Peritoneal1"  = "GSE128830_Peritoneal",
  "Microglia3"   = "GSE137028_Microglia",
  "Adipose2"     = "GSE154832_eWAT",
  "Callus1"      = "PRJNA682234_Callus",
  "Microglia4"   = "GSE267529_Microglia_M",
  "Alveolar2"    = "GSE145295_Alveolar",
  "SkM1"         = "GSE142580_SkM",
  "Alveolar3"    = "GSE190689_Alveolar",
  "SkM2"         = "GSE195507_SkM",
  "BMDM1"        = "GSE279654_BMDM",
  "Microglia5"   = "GSE127542_Microglia",
  "Callus2"      = "GSE198666_Callus",
  
  # Female
  "Skin1"        = "GSE199763_SkinWound",
  "Spleen2"      = "GSE199879_Spleen_Red_Pulp",
  "Microglia6"   = "GSE156762_Microglia",
  "Microglia7"   = "GSE267529_Microglia_F",
  "Microglia8"   = "GSE131869_Microglia_F",
  "BMDM2"        = "PRJNA1029936",
  "Alveolar4"    = "GSE134397_Alveolar"
)

desired_alias_order <- c(
  # Adipose
  "GSE154832_eWAT",
  "GSE93202_VAT",
  
  # Alveolar
  "GSE134397_Alveolar",
  "GSE134397_Alveolar_CTL",
  "GSE145295_Alveolar",
  "GSE190689_Alveolar",
  
  # BMDM
  "PRJNA1029936",
  "GSE279654_BMDM",
  
  # Callus
  "PRJNA682234_Callus",
  "GSE198666_Callus",
  
  # Microglia
  "GSE98401_Microglia",
  "GSE131869_Microglia_M",
  "GSE131869_Microglia_F",
  "GSE137028_Microglia",
  "GSE156762_Microglia",
  "GSE267529_Microglia_M",
  "GSE267529_Microglia_F",
  "GSE127542_Microglia",
  
  
  # Peritoneal
  "GSE128830_Peritoneal",
  
  # Skeletal
  "GSE142580_SkM",
  "GSE195507_SkM",
  
  # Skin
  "GSE199763_SkinWound",
  
  # Spleen
  "GSE93202_Spleen",
  "GSE199879_Spleen_Red_Pulp"
)



################################################################################
# 3. Configure the gene lists and run Mitch
################################################################################

# 1. correct the formatting of DESeq2 output
my_deseq_list <- lapply(my_deseq_list, function(df) {
  rownames(df) <- df$row.names  # Set gene names as row names
  df$row.names <- NULL         
  df
})

# 2. initialize mitch input
mitch_input <- mitch_import(my_deseq_list, DEtype = "deseq2")


# Note: Mean no. genes in input = 14493.2083333333
# Note: no. genes in output = 9225
# Note: estimated proportion of input genes in output = 0.637


# 3. change the column names to the alias names so they don't register as duplicate
colnames(mitch_input) <- dataset.aliases[colnames(mitch_input)]


# 4. Get Reactome, KEGG, and GOBP gene sets for mouse 
mouse_reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
mouse_kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG_LEGACY")
mouse_GOBP <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")


# 5. Convert to a list format (named list of gene vectors)
genesets_mouse_reactome <- split(mouse_reactome$gene_symbol, mouse_reactome$gs_name)
genesets_mouse_kegg <- split(mouse_kegg$gene_symbol, mouse_kegg$gs_description)
genesets_mouse_GOBP <- split(mouse_GOBP$gene_symbol, mouse_GOBP$gs_name)



# 6a. run Mitch for each database - prioritize significance
## reactome pathway analysis 
res_reactome_sig <- mitch_calc(mitch_input,genesets_mouse_reactome, priority="significance",cores=2)

## kegg pathway analysis 
res_kegg_sig <- mitch_calc(mitch_input,genesets_mouse_kegg, priority="significance",cores=2)

## gene ontology pathway analysis 
res_GOBP_sig <- mitch_calc(mitch_input,genesets_mouse_GOBP, priority="significance",cores=2)



# 6b. run Mitch for each database - prioritize effect size
## reactome pathway analysis 
res_reactome_effect <- mitch_calc(mitch_input,genesets_mouse_reactome, priority="effect",cores=2)

## kegg pathway analysis 
res_kegg_effect <- mitch_calc(mitch_input,genesets_mouse_kegg, priority="effect",cores=2)

## gene ontology pathway analysis 
res_GOBP_effect <- mitch_calc(mitch_input,genesets_mouse_GOBP, priority="effect",cores=2)



# 7. Save all intermediate Mitch results (significance + effect)

# Function to rename columns in the enrichment_result table
rename_sp_columns <- function(df, alias_to_dataset) {
  
  old_names <- names(df)
  
  new_names <- sapply(old_names, function(col) {
    # Only rename columns starting with s. or p.
    if (grepl("^(s|p)\\.", col)) {
      prefix <- substr(col, 1, 2)           # "s." or "p."
      alias  <- sub("^(s|p)\\.", "", col)   # e.g., "Microglia4"
      
      if (alias %in% names(alias_to_dataset)) {
        # Replace alias with actual dataset name
        return(paste0(prefix, alias_to_dataset[alias]))
      }
    }
    # Otherwise return unchanged
    return(col)
  })
  
  names(df) <- new_names
  return(df)
}



## --- Reactome ---
# Apply renaming before writing CSVs
reactome_sig_df <- rename_sp_columns(res_reactome_sig$enrichment_result, alias_to_dataset)
reactome_effect_df <- rename_sp_columns(res_reactome_effect$enrichment_result, alias_to_dataset)


# Significance
write.csv(reactome_sig_df, 
          file = "reactome_mitch_results_significance.csv", row.names = FALSE)

saveRDS(res_reactome_sig, 
        file = "res_reactome_significance.rds")

# Effect
write.csv(res_reactome_effect$enrichment_result, 
          file = "reactome_mitch_results_effect.csv", row.names = FALSE)

saveRDS(res_reactome_effect, 
        file = "res_reactome_effect.rds")


## --- KEGG ---
# Apply renaming before writing CSVs
kegg_sig_df    <- rename_sp_columns(res_kegg_sig$enrichment_result, alias_to_dataset)
kegg_effect_df <- rename_sp_columns(res_kegg_effect$enrichment_result, alias_to_dataset)

# Significance
write.csv(kegg_sig_df,
          file = "kegg_mitch_results_significance.csv",
          row.names = FALSE)

saveRDS(res_kegg_sig,
        file = "res_kegg_significance.rds")

# Effect
write.csv(kegg_effect_df,
          file = "kegg_mitch_results_effect.csv",
          row.names = FALSE)

saveRDS(res_kegg_effect,
        file = "res_kegg_effect.rds")


## --- GO Biological Process (GOBP) ---
# Apply renaming before writing CSVs
GOBP_sig_df    <- rename_sp_columns(res_GOBP_sig$enrichment_result, alias_to_dataset)
GOBP_effect_df <- rename_sp_columns(res_GOBP_effect$enrichment_result, alias_to_dataset)

# Significance
write.csv(GOBP_sig_df,
          file = "GOBP_mitch_results_significance.csv",
          row.names = FALSE)

saveRDS(res_GOBP_sig,
        file = "res_GOBP_significance.rds")

# Effect
write.csv(GOBP_effect_df,
          file = "GOBP_mitch_results_effect.csv",
          row.names = FALSE)

saveRDS(res_GOBP_effect,
        file = "res_GOBP_effect.rds")




# 8. Load all intermediate Mitch results both significance and effect 

## --- Reactome ---
res_reactome_sig <- readRDS("res_reactome_significance.rds")
res_reactome_effect <- readRDS("res_reactome_effect.rds")

## --- KEGG ---
res_kegg_sig <- readRDS("res_kegg_significance.rds")
res_kegg_effect <- readRDS("res_kegg_effect.rds")

## --- GO Biological Process (GOBP) ---
res_GOBP_sig <- readRDS("res_GOBP_significance.rds")
res_GOBP_effect <- readRDS("res_GOBP_effect.rds")


################################################################################
# 4. Plot heatmap of top 10 pathways either by significance or effect
################################################################################

### No consistency, just top 10

plot_heatmap_with_fdr <- function(enrichment_result, score_range = c(-1, 1), top_n = 10, title = "Pathways Heatmap") {
  
  # 1. Extract enrichment score columns
  score_cols <- grep("^s\\.", colnames(enrichment_result), value = TRUE)
  score_cols <- setdiff(score_cols, "s.dist")
  
  
  # 2. Restrict to top N pathways
  top_paths <- enrichment_result[1:min(top_n, nrow(enrichment_result)), ]
  
  # 3. Create matrix for heatmap
  mat <- as.matrix(top_paths[, score_cols])
  rownames(mat) <- top_paths$set
  colnames(mat) <- gsub("^s\\.", "", colnames(mat))
  
  # Convert pathway names to lowercase
  rownames(mat) <- tolower(top_paths$set)
  
  # 3.5 Change col names back to dataset name and sort
    # Rename aliases to real dataset names
    colnames(mat) <- alias_to_dataset[colnames(mat)]
  
    # Reorder columns to match desired_alias_order
    common_cols <- intersect(desired_alias_order, colnames(mat))
    mat <- mat[, common_cols, drop = FALSE]
  
  # 4. Create barplot for -log10(FDR)
  fdr_values <- -log10(top_paths$p.adjustMANOVA)
  
  # impute a nonzero/not infinite value for very small p-values
  cap <- 300
  fdr_values[is.infinite(fdr_values)] <- cap
  
  names(fdr_values) <- top_paths$set
  row_ha <- rowAnnotation(`-log10(FDR)` = anno_barplot(fdr_values, border = FALSE, gp = gpar(fill = "black")))
  
  # 5. Draw heatmap with row annotation
  Heatmap(
    mat,
    name = "Enrichment\nScore",
    col = colorRamp2(c(-1, 0, 1), c("darkslateblue", "white", "firebrick3")),
    top_annotation = NULL,
    right_annotation = row_ha,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(title = "Enrichment\nScore"),
    row_title = NULL,
    row_names_side = "left",
    rect_gp = gpar(col = "black", lwd = 0.5),  # add gridlines between cells
    cluster_rows = T,
    cluster_columns = FALSE
  )
}


################################################################################################
# Reactome
pdf(paste0(Sys.Date(),"_Reactome_top10_significance.pdf"), width = 20, height = 11)
plot_heatmap_with_fdr(res_reactome_sig$enrichment_result, top_n = 10)
dev.off()

pdf(paste0(Sys.Date(),"_Reactome_top10_effect.pdf"), width = 20, height = 11)
plot_heatmap_with_fdr(res_reactome_effect$enrichment_result, top_n = 10)
dev.off()

################################################################################################
# KEGG
pdf(paste0(Sys.Date(),"_KEGG_top10_significance.pdf"), width = 20, height = 11)
plot_heatmap_with_fdr(res_kegg_sig$enrichment_result, top_n = 10)
dev.off()

pdf(paste0(Sys.Date(),"_KEGG_top10_effect.pdf"), width = 20, height = 11)
plot_heatmap_with_fdr(res_kegg_effect$enrichment_result, top_n = 10)
dev.off()

################################################################################################
# GO Biological Process
pdf(paste0(Sys.Date(),"_GOBP_top10_significance.pdf"), width = 20, height = 11)
plot_heatmap_with_fdr(res_GOBP_sig$enrichment_result, top_n = 10)
dev.off()

pdf(paste0(Sys.Date(),"_GOBP_top10_effect.pdf"), width = 20, height = 11)
plot_heatmap_with_fdr(res_GOBP_effect$enrichment_result, top_n = 10)
dev.off()

################################################################################################

########################################################################################
# 5. Plot heatmap of top 10 pathways consistently upregulated or downregulated with age
########################################################################################


get_consistent_up_down_heatmaps <- function(
    enrichment_result,
    score_range = c(-1, 1),
    top_n = 10,
    consistency_count = 18
) {
  # 1. Extract enrichment score columns
  score_cols <- grep("^s\\.", colnames(enrichment_result), value = TRUE)
  score_cols <- setdiff(score_cols, "s.dist")
  if (length(score_cols) == 0)
    stop("No columns starting with 's.' found in enrichment_result.")
  
  # Convert to numeric 
  enrichment_result[, score_cols] <- lapply(enrichment_result[, score_cols, drop = FALSE], function(x) as.numeric(as.character(x)))
  
  # 2. Determine consistency (≥ consistency_count datasets with same sign)
  consistency_sign <- apply(enrichment_result[, score_cols, drop = FALSE], 1, function(x) {
    pos <- sum(x > 0, na.rm = TRUE)
    neg <- sum(x < 0, na.rm = TRUE)
    if (pos >= consistency_count) return(1)
    if (neg >= consistency_count) return(-1)
    return(NA)
  })
  enrichment_result$consistency_sign <- consistency_sign
  
  # 3. Subset consistent ones in original order
  consistent_up   <- enrichment_result[enrichment_result$consistency_sign ==  1, ]
  consistent_down <- enrichment_result[enrichment_result$consistency_sign == -1, ]
  
  # Remove rows where 'set' is NA
  consistent_up <- consistent_up[!is.na(consistent_up$set), ]
  consistent_down <- consistent_down[!is.na(consistent_down$set), ]
  
  # 4. Take top N from each group, preserving input order
  top_up   <- head(consistent_up,   n = min(top_n, nrow(consistent_up)))
  top_down <- head(consistent_down, n = min(top_n, nrow(consistent_down)))
  
  # 5. Helper to make annotated heatmap
  make_heatmap <- function(top_paths, title) {
    if (nrow(top_paths) == 0) return(NULL)
    
    # Extract numeric matrix
    mat <- as.matrix(top_paths[, score_cols, drop = FALSE])
    storage.mode(mat) <- "numeric"
    
    # Drop all-NA rows first
    keep <- rowSums(!is.na(mat)) > 0
    mat <- mat[keep, , drop = FALSE]
    top_paths <- top_paths[keep, , drop = FALSE]
    
    # Assign dimnames after filtering
    rownames(mat) <- top_paths$set
    colnames(mat) <- gsub("^s\\.", "", colnames(mat))
    
    # Change col names back to dataset name and sort
    # Rename aliases to real dataset names
    colnames(mat) <- alias_to_dataset[colnames(mat)]
    
    # Reorder columns to match desired_alias_order
    common_cols <- intersect(desired_alias_order, colnames(mat))
    mat <- mat[, common_cols, drop = FALSE]
    
    rownames(mat) <- tolower(top_paths$set)              # lowercase pathway names
    
    # -log10(FDR) annotation
    fdr_values <- -log10(top_paths$p.adjustMANOVA)
    
    # impute a non-zero/not infinite value for very small p-values
    cap <- 300
    fdr_values[is.infinite(fdr_values)] <- cap
    
    names(fdr_values) <- top_paths$set
    row_ha <- rowAnnotation(`-log10(FDR)` = anno_barplot(
      fdr_values,
      border = FALSE,
      width = unit(4, "cm"),
      gp = gpar(fill = "black"),
      ylim = c(0, 12),  # fixed y-axis limits
      axis = TRUE,     # show axis
      axis_param = list(at = seq(0, 12, by = 2), labels = seq(0, 12, by = 2), gp = gpar(fontsize = 6))
    ))
    
    # Heatmap
    Heatmap(
      mat,
      name = "Enrichment\nScore",
      col = colorRamp2(c(-1, 0, 1), c("darkslateblue", "white", "firebrick3")),
      right_annotation = row_ha,
      show_row_names = TRUE,
      row_names_side = "left",
      show_column_names = TRUE,
      column_names_rot = 45,
      row_names_gp = gpar(fontsize = 10),
      heatmap_legend_param = list(title = "Enrichment\nScore"),
      rect_gp = gpar(col = "black", lwd = 0.5),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title = NULL
    )
  }
  
  # 6. Generate plots
  hm_up <- make_heatmap(top_up, title_up)
  hm_down <- make_heatmap(top_down, title_down)
  
  message(sprintf(
    "Found %d consistent up and %d consistent down pathways (showing top %d each).",
    nrow(consistent_up), nrow(consistent_down), top_n
  ))
  
  return(list(up = hm_up, down = hm_down))
}



################################################################################################
# GO Biological Process
################################################################################################

pdf(paste0(Sys.Date(), "_GOBP_top10_significance_UP.pdf"), width = 10, height = 4)
heatmaps <- get_consistent_up_down_heatmaps(res_GOBP_sig$enrichment_result, top_n = 10)
draw(heatmaps$up)

dev.off()

pdf(paste0(Sys.Date(), "_GOBP_top10_significance_DOWN.pdf"), width = 10, height = 4)
draw(heatmaps$down)
dev.off()

pdf(paste0(Sys.Date(), "_GOBP_top10_effect_UP.pdf"), width = 10, height = 4)
heatmaps <- get_consistent_up_down_heatmaps(res_GOBP_effect$enrichment_result, top_n = 10)
draw(heatmaps$up)
dev.off()

pdf(paste0(Sys.Date(), "_GOBP_top10_effect_DOWN.pdf"), width = 10, height = 4)
draw(heatmaps$down)
dev.off()


################################################################################################
# KEGG
################################################################################################

pdf(paste0(Sys.Date(), "_KEGG_top10_significance_UP.pdf"), width = 10, height = 4)
heatmaps <- get_consistent_up_down_heatmaps(res_kegg_sig$enrichment_result, top_n = 10)
draw(heatmaps$up)
dev.off()

pdf(paste0(Sys.Date(), "_KEGG_top10_significance_DOWN.pdf"), width = 10, height = 4)
draw(heatmaps$down)
dev.off()

pdf(paste0(Sys.Date(), "_KEGG_top10_effect_UP.pdf"), width = 10, height = 4)
heatmaps <- get_consistent_up_down_heatmaps(res_kegg_effect$enrichment_result, top_n = 10)
draw(heatmaps$up)
dev.off()

pdf(paste0(Sys.Date(), "_KEGG_top10_effect_DOWN.pdf"), width = 10, height = 4)
draw(heatmaps$down)
dev.off()


################################################################################################
# Reactome
################################################################################################

pdf(paste0(Sys.Date(), "_Reactome_top10_significance_UP.pdf"), width = 10, height = 4)
heatmaps <- get_consistent_up_down_heatmaps(res_reactome_sig$enrichment_result, top_n = 10)
draw(heatmaps$up)
dev.off()

pdf(paste0(Sys.Date(), "_Reactome_top10_significance_DOWN.pdf"), width = 10, height = 4)
draw(heatmaps$down)
dev.off()

pdf(paste0(Sys.Date(), "_Reactome_top10_effect_UP.pdf"), width = 10, height = 4)
heatmaps <- get_consistent_up_down_heatmaps(res_reactome_effect$enrichment_result, top_n = 10)
draw(heatmaps$up)
dev.off()

pdf(paste0(Sys.Date(), "_Reactome_top10_effect_DOWN.pdf"), width = 10, height = 4)
draw(heatmaps$down)
dev.off()



################################################################################
# Save package versions
sink(file = paste(Sys.Date(),"Mitch_Session_Info.txt", sep =""))
sessionInfo()
sink()




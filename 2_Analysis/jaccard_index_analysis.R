####################################################################################
# Aging transcriptomics of mouse macrophages 
# Jaccard index plots, gene venn diagrams for alveolar and microglia, ORA plots
####################################################################################

set.seed(1234) # set seed for reproducibility

################################################################################
# 1. Read in the data
################################################################################

library(pheatmap)
library(stringr)
library(VennDiagram)
library(scales)
library(clusterProfiler)
library(org.Mm.eg.db)  
library(ggplot2)


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


################################################################################
# 2. Reorder based on niche
################################################################################

setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/jaccard_index")

padj_cutoff <- 0.05
log2fc_cutoff <- 0  


male_datasets <- c(
  "GSE93202_Spleen", "GSE93202_VAT", "GSE98401_Microglia", "GSE134397_Alveolar_CTL",
  "GSE131869_Microglia_M", "GSE128830_Peritoneal", "GSE137028_Microglia",
  "GSE154832_eWAT", "PRJNA682234_Callus", "GSE267529_Microglia_M",
  "GSE145295_Alveolar", "GSE142580_SkM", "GSE190689_Alveolar", 
  "PRJNA800823_SkM", "PRJNA1173774_BMDM", "PRJNA524906_Microglia", 
  "PRJNA816431_Callus"
)

female_datasets <- c(
  "GSE199763_SkinWound", "GSE199879_Spleen_Red_Pulp", "GSE156762_Microglia",
  "GSE267529_Microglia_F", "GSE131869_Microglia_F",
  "BMDM_NIA_F", "GSE134397_Alveolar" 
)


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


dataset.aliases <- list(
  # Male datasets
  "GSE93202_Spleen"           = "GSE93202_Spleen",
  "GSE93202_VAT"              = "GSE93202_VAT",
  "GSE98401_Microglia"        = "GSE98401_Microglia",
  "GSE134397_Alveolar_CTL"    = "GSE134397_Alveolar_CTL",
  "GSE131869_Microglia_M"     = "GSE131869_Microglia_M",
  "GSE128830_Peritoneal"      = "GSE128830_Peritoneal",
  "GSE137028_Microglia"       = "GSE137028_Microglia",
  "GSE154832_eWAT"            = "GSE154832_eWAT",
  "PRJNA682234_Callus"        = "PRJNA682234_Callus", 
  "GSE267529_Microglia_M"     = "GSE267529_Microglia_M",
  "GSE145295_Alveolar"        = "GSE145295_Alveolar",
  "GSE142580_SkM"             = "GSE142580_SkM",
  "GSE190689_Alveolar"        = "GSE190689_Alveolar",
  "PRJNA800823_SkM"           = "GSE195507_SkM",
  "PRJNA1173774_BMDM"         = "GSE279654_BMDM",
  "PRJNA524906_Microglia"     = "GSE127542_Microglia",
  "PRJNA816431_Callus"        = "GSE198666_Callus",
  
  # Female datasets
  "GSE199763_SkinWound"       = "GSE199763_SkinWound",
  "GSE199879_Spleen_Red_Pulp" = "GSE199879_Spleen_Red_Pulp",
  "GSE156762_Microglia"       = "GSE156762_Microglia",
  "GSE267529_Microglia_F"     = "GSE267529_Microglia_F",
  "GSE131869_Microglia_F"     = "GSE131869_Microglia_F",
  "BMDM_NIA_F"                = "PRJNA1029936",
  "GSE134397_Alveolar"        = "GSE134397_Alveolar"
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
# 3. Make Jaccard Index Plots - by FDR < 5%
################################################################################

# Split DEG sets by direction
young_sets <- list()
old_sets <- list()


# Loop over my_deseq_list with names
for (dataset_name in names(my_deseq_list)) {
  df <- my_deseq_list[[dataset_name]]
  
  # Filter for significant DEGs
  df <- df[!is.na(df$padj) & df$padj < padj_cutoff, ]
  df <- df[!is.na(df$log2FoldChange), ]
  
  # Split by log2 fold change direction
  old_sets[[dataset_name]] <- df[df$log2FoldChange >  log2fc_cutoff, "row.names"]
  young_sets[[dataset_name]] <- df[df$log2FoldChange < -log2fc_cutoff, "row.names"]
}

### Save young and old sets (FDR cutoff)
make_padded_df <- function(dataset_list) {
  max_len <- max(sapply(dataset_list, length))
  padded <- lapply(dataset_list, function(x) {
    length(x) <- max_len      
    x[is.na(x)] <- ""         
    return(x)
  })
  df <- as.data.frame(padded, stringsAsFactors = FALSE)
  return(df)
}


young_df <- make_padded_df(young_sets)
old_df   <- make_padded_df(old_sets)

date_stamp <- Sys.Date()  

# Save CSVs with date in filename
write.csv(young_df, paste0(date_stamp,"_young_sets.csv"), row.names = FALSE)
write.csv(old_df,   paste0(date_stamp,"_old_sets.csv"), row.names = FALSE)

####


# Intersect datasets with valid data in both sets
dataset_names <- intersect(names(old_sets), names(young_sets))

# Sex annotation dataframe
sex_vector <- sapply(dataset_names, function(x) {
  if (x %in% male_datasets) {
    return("Male")
  } else if (x %in% female_datasets) {
    return("Female")
  } else {
    return(NA)
  }
})

# Create annotation dataframe
annotation_col <- data.frame(Sex = sex_vector)
rownames(annotation_col) <- dataset_names

# Define color mapping
sex_colors <- list(Sex = c(Male = "deepskyblue", Female = "deeppink"))

# Jaccard calculation
compute_jaccard_matrix <- function(dataset_names, gene_sets) {
  mat <- matrix(0, nrow = length(dataset_names), ncol = length(dataset_names),
                dimnames = list(dataset_names, dataset_names))
  
  for (i in dataset_names) {
    for (j in dataset_names) {
      set1 <- gene_sets[[i]]
      set2 <- gene_sets[[j]]
      intersection <- length(intersect(set1, set2))
      union <- length(union(set1, set2))
      mat[i, j] <- if (union > 0) intersection / union else 0
    }
  }
  return(mat)
}

# Compute matrices
old_jaccard <- compute_jaccard_matrix(dataset_names, old_sets)
young_jaccard <- compute_jaccard_matrix(dataset_names, young_sets)

# Convert alias list to named vector
dataset_aliases_vec <- unlist(dataset.aliases)

# Rename axes
rownames(old_jaccard)   <- dataset_aliases_vec[rownames(old_jaccard)]
colnames(old_jaccard)   <- dataset_aliases_vec[colnames(old_jaccard)]
rownames(young_jaccard) <- dataset_aliases_vec[rownames(young_jaccard)]
colnames(young_jaccard) <- dataset_aliases_vec[colnames(young_jaccard)]

# Rename annotation rownames
rownames(annotation_col) <- dataset_aliases_vec[rownames(annotation_col)]

# Plot heatmaps with sex annotations
down_palette <- colorRampPalette(c("white", "darkslateblue"))(100)
up_palette <- colorRampPalette(c("white", "firebrick3"))(100)


breaks <- c(
  seq(0, 0.1, length.out = 20),         
  seq(0.100001, 0.4, length.out = 60),   
  seq(0.400001, 1, length.out = 21)      
)

length(breaks)  

# --- Ensure desired order matches datasets present ---
ordered_aliases <- desired_alias_order[
  desired_alias_order %in% rownames(old_jaccard)
]

# Reorder matrices
old_jaccard   <- old_jaccard[ordered_aliases, ordered_aliases]
young_jaccard <- young_jaccard[ordered_aliases, ordered_aliases]

# Reorder annotation dataframe to match
annotation_col <- annotation_col[ordered_aliases, , drop = FALSE]

# --- Plot heatmaps without clustering ---
pdf(paste0(Sys.Date(), "_Jaccard_DEG_Overlap_OldUp_Annotated_ordered.pdf"), width = 17, height = 14)
pheatmap(old_jaccard,
         main = "Upregulated Genes With Age (padj < 0.05)",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         annotation_colors = sex_colors,
         color = up_palette,
         breaks = breaks,
         display_numbers = FALSE,
         number_format = "%.2f")
dev.off()

pdf(paste0(Sys.Date(), "_Jaccard_DEG_Overlap_YoungUp_Annotated_ordered.pdf"), width = 17, height = 14)
pheatmap(young_jaccard,
         main = "Downregulated Genes With Age (padj < 0.05)",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         annotation_colors = sex_colors,
         color = down_palette,
         breaks = breaks,
         display_numbers = FALSE,
         number_format = "%.2f")
dev.off()



################################################################################
# 4. Make Jaccard Index Plots - by 1500 upregulated and 1500 downregulated genes 
################################################################################

# Initialize lists to store top genes
top_young_sets <- list()
top_old_sets <- list()

# Loop over each dataset in your list
for (dataset_name in names(my_deseq_list)) {
  df <- my_deseq_list[[dataset_name]]
  
  df <- df[!is.na(df$log2FoldChange), ]
  
  # Sort by log2FC for top 1500 up and down
  top_old <- df[order(df$log2FoldChange, decreasing = TRUE), ][1:min(1500, nrow(df)), ]
  top_young <- df[order(df$log2FoldChange, decreasing = FALSE), ][1:min(1500, nrow(df)), ]
  
  # Store row names (gene IDs)
  top_old_sets[[dataset_name]] <- top_old[["row.names"]]
  top_young_sets[[dataset_name]] <- top_young[["row.names"]]
  
}

### Save young and old sets (FDR cutoff)
make_padded_df <- function(dataset_list) {
  max_len <- max(sapply(dataset_list, length))
  padded <- lapply(dataset_list, function(x) {
    length(x) <- max_len      
    x[is.na(x)] <- ""         
    return(x)
  })
  df <- as.data.frame(padded, stringsAsFactors = FALSE)
  return(df)
}

# Create data frames
young_df <- make_padded_df(top_young_sets)
old_df   <- make_padded_df(top_old_sets)

# Save as CSV files
date_stamp <- Sys.Date()  

# Save CSVs with date in filename
write.csv(young_df, paste0(date_stamp,"_top_young_sets.csv"), row.names = FALSE)
write.csv(old_df,   paste0(date_stamp,"_top_old_sets.csv"), row.names = FALSE)

####

# Find intersecting dataset names that exist in both
dataset_names <- intersect(names(top_old_sets), names(top_young_sets))

# Build sex annotation vector based on dataset name
sex_vector <- sapply(dataset_names, function(x) {
  if (x %in% male_datasets) {
    return("Male")
  } else if (x %in% female_datasets) {
    return("Female")
  } else {
    return(NA)
  }
})

# Create annotation dataframe
annotation_col <- data.frame(Sex = sex_vector)
rownames(annotation_col) <- dataset_names

# Define sex color mapping
sex_colors <- list(Sex = c(Male = "deepskyblue", Female = "deeppink"))

# Reusable Jaccard matrix computation function
compute_jaccard_matrix <- function(dataset_names, gene_sets) {
  mat <- matrix(0, nrow = length(dataset_names), ncol = length(dataset_names),
                dimnames = list(dataset_names, dataset_names))
  
  for (i in dataset_names) {
    for (j in dataset_names) {
      set1 <- gene_sets[[i]]
      set2 <- gene_sets[[j]]
      intersection <- length(intersect(set1, set2))
      union <- length(union(set1, set2))
      mat[i, j] <- if (union > 0) intersection / union else 0
    }
  }
  return(mat)
}

# Compute Jaccard matrices using top 1500 genes
old_jaccard <- compute_jaccard_matrix(dataset_names, top_old_sets)
young_jaccard <- compute_jaccard_matrix(dataset_names, top_young_sets)

# Convert alias list to named vector
dataset_aliases_vec <- unlist(dataset.aliases)

# Rename axes
rownames(old_jaccard)   <- dataset_aliases_vec[rownames(old_jaccard)]
colnames(old_jaccard)   <- dataset_aliases_vec[colnames(old_jaccard)]
rownames(young_jaccard) <- dataset_aliases_vec[rownames(young_jaccard)]
colnames(young_jaccard) <- dataset_aliases_vec[colnames(young_jaccard)]

# Rename annotation rownames
rownames(annotation_col) <- dataset_aliases_vec[rownames(annotation_col)]

breaks <- c(
  seq(0, 0.1, length.out = 20),         
  seq(0.100001, 0.4, length.out = 60),   
  seq(0.400001, 1, length.out = 21)      
)

# Plot heatmaps with sex annotations
down_palette <- colorRampPalette(c("white", "darkslateblue"))(100)
up_palette <- colorRampPalette(c("white", "firebrick3"))(100)


length(breaks)  

# --- Ensure desired order matches datasets present ---
ordered_aliases <- desired_alias_order[
  desired_alias_order %in% rownames(old_jaccard)
]

# Reorder matrices
old_jaccard   <- old_jaccard[ordered_aliases, ordered_aliases]
young_jaccard <- young_jaccard[ordered_aliases, ordered_aliases]

# Reorder annotation dataframe to match
annotation_col <- annotation_col[ordered_aliases, , drop = FALSE]


# Old-upregulated heatmap
pdf(paste0(Sys.Date(),"_Jaccard_TOP1500_DEG_Overlap_OldUp_Annotated.pdf"), width = 17, height = 14)
pheatmap(old_jaccard,
         main = "Top 1500 Upregulated Genes With Age",
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = annotation_col,
         annotation_colors = sex_colors,
         color = up_palette,
         breaks = breaks,
         display_numbers = F,
         number_format = "%.2f")
dev.off()

# Young-upregulated heatmap
pdf(paste0(Sys.Date(),"_Jaccard_TOP1500_DEG_Overlap_YoungUp_Annotated_Updated.pdf"), width = 17, height = 14)
pheatmap(young_jaccard,
         main = "Top 1500 Downregulated Genes With Age",
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = annotation_col,
         annotation_colors = sex_colors,
         color = down_palette,
         breaks = breaks,
         display_numbers = F,
         number_format = "%.2f")
dev.off()



################################################################################
# 5. Make Venn Diagram Plots for Alveolar and Microglia sets (FDR < 5%)
################################################################################

# Get male/female alveolar datasets
alveolar_male <- intersect(niche_map$Alveolar, male_datasets)
alveolar_female <- intersect(niche_map$Alveolar, female_datasets)

# Get alveolar dataset names from both sexes
alveolar_datasets <- niche_map$Alveolar

# Pull all genes with non-NA p-values from each dataset
alveolar_universe <- unique(unlist(
  lapply(my_deseq_list[alveolar_datasets], function(df) {
    df$row.names[!is.na(df$padj)]
  })
))


# Get male/female microglia datasets
microglia_male <- intersect(niche_map$Microglia, male_datasets)
microglia_female <- intersect(niche_map$Microglia, female_datasets)

# Get microglia dataset names from both sexes
microglia_datasets <- niche_map$Microglia

# Pull all genes with non-NA p-values from each dataset
microglia_universe <- unique(unlist(
  lapply(my_deseq_list[microglia_datasets], function(df) {
    df$row.names[!is.na(df$padj)]
  })
))



### p value adjusted 
# Alveolar
alveolar_up_male <- unique(unlist(old_sets[names(old_sets) %in% alveolar_male]))
alveolar_up_female <- unique(unlist(old_sets[names(old_sets) %in% alveolar_female]))

alveolar_down_male <- unique(unlist(young_sets[names(young_sets) %in% alveolar_male]))
alveolar_down_female <- unique(unlist(young_sets[names(young_sets) %in% alveolar_female]))

# Microglia
microglia_up_male <- unique(unlist(old_sets[names(old_sets) %in% microglia_male]))
microglia_up_female <- unique(unlist(old_sets[names(old_sets) %in% microglia_female]))

microglia_down_male <- unique(unlist(young_sets[names(young_sets) %in% microglia_male]))
microglia_down_female <- unique(unlist(young_sets[names(young_sets) %in% microglia_female]))




myCol <- c("deepskyblue", "deeppink")


#### function to save ora plots and csv

# Function to create Venn diagram and save gene sets to one CSV
create_venn_and_save_csv <- function(set_male, set_female, title, filename_prefix) {
  
  # Generate the Venn diagram
  venn.diagram(
    x = list(
      Male = set_male,
      Female = set_female
    ),
    filename = paste0(filename_prefix, "_venn.png"),
    imagetype = "png",
    fill = myCol,
    alpha = 0.2,
    lwd = 0,
    col = "transparent",
    cex = 2.75,
    cat.cex = 1.5,
    cat.col = "black",
    cat.fontface = "bold",
    fontfamily = "Arial",
    cat.fontfamily = "Arial",
    main.fontfamily = "Arial",
    cat.pos = c(1.5, 1.5),
    cat.dist = c(0.02, 0.02),
    cat.default.pos = "outer",
    margin = 0.05,
    main = title,
    main.pos = c(0.5, 0.95),
    main.cex = 2,
    height = 3000,
    width = 3000,
    resolution = 300,
    disable.logging = TRUE
  )
  
  # Find shared and unique genes
  male_only <- setdiff(set_male, set_female)
  female_only <- setdiff(set_female, set_male)
  shared_genes <- intersect(set_male, set_female)
  
  # Combine into one data frame (pad shorter vectors with NA)
  max_len <- max(length(male_only), length(female_only), length(shared_genes))
  df <- data.frame(
    male_only = c(male_only, rep(NA, max_len - length(male_only))),
    female_only = c(female_only, rep(NA, max_len - length(female_only))),
    shared_genes = c(shared_genes, rep(NA, max_len - length(shared_genes)))
  )
  
  # Save one CSV
  write.csv(df, paste0(filename_prefix, "_gene_sets.csv"), row.names = FALSE)
}


# Create all four Venn diagrams and save gene sets
create_venn_and_save_csv(microglia_up_male, microglia_up_female, 
                         "Microglia: Male vs Female Upregulated Genes", 
                         "microglia_up")

create_venn_and_save_csv(microglia_down_male, microglia_down_female, 
                         "Microglia: Male vs Female Downregulated Genes", 
                         "microglia_down")

create_venn_and_save_csv(alveolar_up_male, alveolar_up_female, 
                         "Alveolar: Male vs Female Upregulated Genes", 
                         "alveolar_up")

create_venn_and_save_csv(alveolar_down_male, alveolar_down_female, 
                         "Alveolar: Male vs Female Downregulated Genes", 
                         "alveolar_down")



################################################################################
# 6. Make ORA plots for Alveolar and Microglia sets (FDR < 5%)
################################################################################


# 1. Extract sets
# Alveolar Up
shared_alveolar_up <- intersect(alveolar_up_male, alveolar_up_female)
male_only_alveolar_up <- setdiff(alveolar_up_male, alveolar_up_female)
female_only_alveolar_up <- setdiff(alveolar_up_female, alveolar_up_male)

# Alveolar Down
shared_alveolar_down <- intersect(alveolar_down_male, alveolar_down_female)
male_only_alveolar_down <- setdiff(alveolar_down_male, alveolar_down_female)
female_only_alveolar_down <- setdiff(alveolar_down_female, alveolar_down_male)

# Microglia Up
shared_microglia_up <- intersect(microglia_up_male, microglia_up_female)
male_only_microglia_up <- setdiff(microglia_up_male, microglia_up_female)
female_only_microglia_up <- setdiff(microglia_up_female, microglia_up_male)

# Microglia Down
shared_microglia_down <- intersect(microglia_down_male, microglia_down_female)
male_only_microglia_down <- setdiff(microglia_down_male, microglia_down_female)
female_only_microglia_down <- setdiff(microglia_down_female, microglia_down_male)



# 2. Run ORA
run_ora <- function(genes, universe, title = "ORA", color) {
  
  # Run enrichment
  enrich_result <- enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",  
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
 
  df <- as.data.frame(enrich_result)
  
  # ---- Save results as CSV ----
  if (nrow(df) > 0) {
    safe_title <- gsub("[^A-Za-z0-9_\\-]", "_", title)  
    date_str <- format(Sys.Date(), "%Y%m%d")
    filename <- paste0(date_str, "_", safe_title, "_ORA_results.csv")
    write.csv(df, filename, row.names = FALSE)
  } else {
    message(paste("No significant terms for:", title))
  }
  # ---------------------------------
  
  # Take top 15 categories
  df <- df[order(df$p.adjust), ][1:min(15, nrow(df)), ]
  
  # Gene ratio to numeric 
  if (is.character(df$GeneRatio)) {
    df$GeneRatio <- sapply(df$GeneRatio, function(x) eval(parse(text = x)))
  }
  
  # Take top 15 by adjusted p-value
  df <- df[order(df$p.adjust), ][1:min(15, nrow(df)), ]
  
  # Reorder factor levels by GeneRatio (ascending)
  df$Description <- factor(df$Description, levels = df$Description[order(df$GeneRatio)])
  
  
  max_log10p <- ceiling(max(-log10(df$p.adjust), na.rm = TRUE))
  
  df$Description <- stringr::str_wrap(df$Description, width = 40) # wrap long names
  
  # Create a light version of the input color
  light_color <- colorspace::lighten(color, amount = 0.7)
  
  # Compute cutoff and max
  cutoff_log10 <- -log10(0.05)
  max_log10p <- ceiling(max(-log10(df$p.adjust), na.rm = TRUE))
  
  ggplot(df, aes(x = 1, y = Description)) +   # x is constant = 1, so all points line up vertically
    geom_point(aes(size = Count, color = -log10(p.adjust))) +
    scale_color_gradientn(
      colours = c("grey", light_color, color),
      values = scales::rescale(c(0, cutoff_log10, 100)),
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
      labels = seq(0, 100, by = 20),
      name = expression(-log[10]("FDR"))
    )+ 
    scale_size_continuous(
      range = c(3, 20),
      limits = c(1, 150),
      breaks = c(25, 50, 100, 150)
    ) + 
    labs(
      title = title,
      x = NULL,        
      y = NULL,
      color = expression(-log[10]("FDR")),
      size = "Gene Count"
    ) +
    guides(
      color = guide_colorbar(order = 1),  
      size = guide_legend(order = 2)      
    )+
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 16),
      legend.key.height = unit(0.7, "cm"),
      legend.key.width = unit(1.2, "cm"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_rect(fill = "white", color = NA),  
      panel.grid.major.y = element_line(color = "#D3D3D3"),        
      plot.background = element_rect(fill = "white", color = NA),   
      plot.title = element_text(face = "bold")
    )
  
}

# 3. Plot each
# Alveolar Up
plot_shared_alveolar_up <- run_ora(shared_alveolar_up, alveolar_universe, 
                                   title = "Alveolar Up - Shared (Both Sexes)", 
                                   color = "#7F69C9")

plot_male_alveolar_up <- run_ora(male_only_alveolar_up, alveolar_universe, 
                                 title = "Alveolar Up - Male Specific", 
                                 color = "deepskyblue")

plot_female_alveolar_up <- run_ora(female_only_alveolar_up, alveolar_universe, 
                                   title = "Alveolar Up - Female Specific", 
                                   color = "deeppink")

# Alveolar Down
plot_shared_alveolar_down <- run_ora(shared_alveolar_down, alveolar_universe, 
                                     title = "Alveolar Down - Shared (Both Sexes)", 
                                     color = "#7F69C9")
plot_male_alveolar_down <- run_ora(male_only_alveolar_down, alveolar_universe, 
                                   title = "Alveolar Down - Male Specific", 
                                   color = "deepskyblue")
plot_female_alveolar_down <- run_ora(female_only_alveolar_down, alveolar_universe, 
                                     title = "Alveolar Down - Female Specific", 
                                     color = "deeppink")

# Microglia Up
plot_shared_microglia_up <- run_ora(shared_microglia_up, microglia_universe, 
                                    title = "Microglia Up - Shared (Both Sexes)", 
                                    color = "#7F69C9")
plot_male_microglia_up <- run_ora(male_only_microglia_up, microglia_universe, 
                                  title = "Microglia Up - Male Specific", 
                                  color = "deepskyblue")
plot_female_microglia_up <- run_ora(female_only_microglia_up, microglia_universe, 
                                    title = "Microglia Up - Female Specific", 
                                    color = "deeppink")

# Microglia Down
plot_shared_microglia_down <- run_ora(shared_microglia_down, microglia_universe, 
                                      title = "Microglia Down - Shared (Both Sexes)", 
                                      color = "#7F69C9")
plot_male_microglia_down <- run_ora(male_only_microglia_down, microglia_universe, 
                                    title = "Microglia Down - Male Specific", 
                                    color = "deepskyblue")
plot_female_microglia_down <- run_ora(female_only_microglia_down, microglia_universe, 
                                      title = "Microglia Down - Female Specific", 
                                      color = "deeppink")




# Put all plots in a named list
plots_list <- list(
  plot_shared_alveolar_up,
  plot_male_alveolar_up,
  plot_female_alveolar_up,
  plot_shared_alveolar_down,
  plot_male_alveolar_down,
  plot_female_alveolar_down,
  plot_shared_microglia_up,
  plot_male_microglia_up,  
  plot_female_microglia_up,
  plot_shared_microglia_down,
  plot_male_microglia_down,
  plot_female_microglia_down
)

plot_names <- c(
  "shared_alveolar_up",
  "male_alveolar_up",
  "female_alveolar_up",
  "shared_alveolar_down",
  "male_alveolar_down",
  "female_alveolar_down",
  "shared_microglia_up",
  "male_microglia_up",  
  "female_microglia_up",
  "shared_microglia_down",
  "male_microglia_down",
  "female_microglia_down"
)

# Append "ORA" to each name
plot_names <- paste0(plot_names, "_ORA")

# Create output directory
output_dir <- "ORA_plots"
dir.create(output_dir, showWarnings = FALSE)

# Save function
save_plots <- function(plot_obj, filename_base) {
  if (is.null(plot_obj)) {
    message(paste("Skipping", filename_base, "- plot is NULL"))
    return()
  }
  
  pdf_file <- file.path(output_dir, paste0(filename_base, ".pdf"))
  png_file <- file.path(output_dir, paste0(filename_base, ".png"))
  
  # Save PDF
  pdf(pdf_file, width = 8, height = 14)
  print(plot_obj)
  dev.off()
  
  # Save PNG
  png(png_file, width = 8, height = 14, units = "in", res = 300)
  print(plot_obj)
  dev.off()
  
  message(paste("Saved:", pdf_file, "and", png_file))
}

# Loop through all plots and save them
for (i in seq_along(plots_list)) {
  save_plots(plots_list[[i]], plot_names[i])
}



########################################################################################################################
sink(file = paste(Sys.Date(),"Jaccard_Index_session_Info.txt", sep =""))
sessionInfo()
sink()



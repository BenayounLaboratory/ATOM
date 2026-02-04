####################################################################################
# Aging transcriptomics of mouse macrophages 
# Gene level aging gene meta-analysis - P-value combination with metaRNAseq
####################################################################################

set.seed(1234) # set seed for reproducibility


################################################################################
# 1. Read in the data
################################################################################


library(pheatmap)
library(stringr)
library(metaRNASeq)
library(dplyr)
library(clusterProfiler)
library(enrichplot)  
library(pheatmap)
library(RColorBrewer)
library(org.Mm.eg.db)
library(colorspace)
library(ggplot2)



my.out.dir <- "/Users/ellaschwab/Benayoun_Lee_Local/ATOM/PValueCombo"

#put into list
setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/2026-01-28_continuous_age_output")
deseq_files <- list.files(pattern = "all_genes_statistics\\.txt$")

my_deseq_list <- list()
# Loop through each file in the directory
for (file in deseq_files) {
  clean_name <- str_replace(file, "^\\d{4}-\\d{2}-\\d{2}_", "")
  dataset_name <- str_replace(clean_name, "_AGE_DIM_all_genes_statistics\\.txt$", "")
  dataset_name <- str_remove(dataset_name, "_$")  
  dataset <- read.csv(file, header = TRUE, row.names = NULL, sep = "\t")
  my_deseq_list[[dataset_name]] <- dataset
}

# Filter for datasets that passed QC
pass_qc <- readLines("2026-01-28_datasets_pass_qc.txt")
my_deseq_list <- my_deseq_list[names(my_deseq_list) %in% pass_qc]


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
# 2. Assess quality of data and filter lowly expressed genes 
################################################################################

setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/PValueCombo")

padj_cutoff <- 0.05
log2fc_cutoff <- 0  


deg_datasets <- my_deseq_list


# Extract p-values 
pval_list <- lapply(deg_datasets, function(df) {
  if ("pvalue" %in% names(df)) {
    return(df[["pvalue"]])
  } else {
    return(NULL)  
  }
})


# Extract logFC 
fc_list <- lapply(deg_datasets, function(df) {
  if ("log2FoldChange" %in% names(df)) {
    return(df[["log2FoldChange"]])
  } else {
    return(NULL)  
  }
})

# Extract adjust p value
adj_pval_list <- lapply(deg_datasets, function(df) {
  if ("padj" %in% names(df)) {
    return(df[["padj"]])
  } else {
    return(NULL)  
  }
})



# Set up plot grid based on number of datasets
num_datasets <- length(deg_datasets)
plot_cols <- ceiling(sqrt(num_datasets))
plot_rows <- ceiling(num_datasets / plot_cols)


par(mar = c(3, 3, 3, 1))   
par(mfrow = c(plot_rows, plot_cols))  # Set up the plot layout

# Loop through each dataset
for (name in names(deg_datasets)) {
  df <- deg_datasets[[name]]
  
  if ("pvalue" %in% names(df)) {
    hist(df$pvalue,
         breaks = 100,
         col = "grey",
         main = name,
         xlab = "p-value distribution")
  } else {
    plot.new()
    title(main = paste(name, "\n(no 'pvalue')"))
  }
}

# Reset plotting layout
par(mfrow = c(1, 1))


# Set basemean threshold to remove lowly expressed gene associated p-values  
threshold <- 50  

# Filter all datasets to remove lowly expressed genes
deg_datasets_filtered <- lapply(deg_datasets, function(df) {
  if ("baseMean" %in% names(df)) {
    df <- df %>% filter(baseMean > threshold)
  }
  return(df)  # Return the filtered or original df
})

# Set up plot grid based on number of datasets
num_datasets <- length(deg_datasets_filtered)
plot_cols <- ceiling(sqrt(num_datasets))
plot_rows <- ceiling(num_datasets / plot_cols)

par(mfrow = c(plot_rows, plot_cols))  # Set up the plot layout

# Loop through each dataset
for (name in names(deg_datasets_filtered)) {
  df <- deg_datasets_filtered[[name]]
  
  if ("pvalue" %in% names(df)) {
    hist(df$pvalue,
         breaks = 100,
         col = "grey",
         main = name,
         xlab = "filtered p-value distribution")
  } else {
    plot.new()
    title(main = paste(name, "\n(no 'pvalue')"))
  }
}




# Create a filtered p-value list
filtered_pval_list <- lapply(deg_datasets, function(df) {
  # Identify the p-value column
  pval_col <- if ("pvalue" %in% names(df)) {
    "pvalue"
  } else if ("P.Value" %in% names(df)) {
    "P.Value"
  } else {
    return(NULL)  # Skip if no p-value column
    #this means that the microarray sets won't get this filtering 
  }
  
  
  if ("baseMean" %in% names(df)) {
    df <- df %>% filter(baseMean > threshold)
  }
  
  # Return filtered p-values
  df[[pval_col]]
})



par(mar = c(3, 3, 3, 1))   # smaller margins
par(mfrow = c(plot_rows, plot_cols))  # Set up the plot layout

for (dataset_name in names(pval_list)) {
  hist(pval_list[[dataset_name]],
       breaks = 100,
       col = "grey",
       main = dataset_name,
       xlab = "Raw p-value distribution")
}


for (dataset_name in names(filtered_pval_list)) {
  hist(filtered_pval_list[[dataset_name]],
       breaks = 100,
       col = "grey",
       main = dataset_name,
       xlab = "Filtered raw p-value distribution")
}


################################################################################
# 3. P-value combination
################################################################################


# 1. Make a short dataframe for each study with the gene names and their associated p value
deg_short <- lapply(deg_datasets, function(df) {
  data.frame(
    gene = df[["row.names"]],       # extract row names as gene column
    log2FoldChange = df$log2FoldChange,  # add log fold change
    pvalue = df$pvalue,          # extract pvalue column
    adjpvalue = df$padj
  )
})



# 2. Remove rows with NA p-values or logFCs
pval_dfs_clean <- lapply(deg_short, function(df) {
  df %>% filter(!is.na(pvalue), !is.na(log2FoldChange))
})

# 3. Remove duplicated genes (keep only the first)
pval_dfs_unique <- lapply(pval_dfs_clean, function(df) {
  df %>% distinct(gene, .keep_all = TRUE)
})

# 4. Get the set of common genes across all studies
common_genes <- Reduce(intersect, lapply(pval_dfs_unique, function(df) df$gene))
length(common_genes)

# [1] 9076

# output the universe to a text file
write.table(common_genes, 
            file = "common_genes_universe.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)


# 5. Filter all datasets to the common genes and order alphabetically
pval_dfs_common <- lapply(pval_dfs_unique, function(df) {
  df %>%
    filter(gene %in% common_genes) %>%
    arrange(gene)
})

# 6. Create raw p-value list
rawpval <- lapply(pval_dfs_common, function(df) df$pvalue)

sample_counts <- read.table("sample_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Names of datasets in rawpval, in order
rawpval_names <- names(rawpval)

# Reorder sample_counts$num to match rawpval order by matching dataset names
dataset_n <- sample_counts$num[match(rawpval_names, sample_counts$dataset)]


# 7. Meta-analysis with Fisher and Inverse Normal methods
fishcomb <- fishercomb(rawpval, BHth = 0.05)
invnormcomb <- invnorm(rawpval, nrep = dataset_n, BHth = 0.05)


# 8. Create fold change matrix for directionality consistency check
fc_list <- lapply(pval_dfs_common, function(df) df$log2FoldChange)
signsFC <- do.call(cbind, lapply(fc_list, sign))

rownames(signsFC) <- pval_dfs_common[[1]]$gene

sumsigns <- rowSums(signsFC, na.rm = TRUE)
commonsgnFC <- ifelse(abs(sumsigns) == rowSums(!is.na(signsFC)), sign(sumsigns), 0) # consistent in all datasets
commonsgn_12 <- ifelse(abs(sumsigns) >= 12, sign(sumsigns), 0) # consistent in half 
commonsgn_18 <- ifelse(abs(sumsigns) >= 18, sign(sumsigns), 0) # consistent in 75%


# this will check if the rowsum is the same as the number of studies represented 
# e.g. 
# Sum  NotNACount Absequal?  Meaning
# 3.      3.        Y.         +1
# -3.     3.        Y.         -1
# 1.      3         N           0


# 9. Identify all DE genes from Fisher or Inverse Normal meta-analysis
unionDE <- unique(c(fishcomb$DEindices, invnormcomb$DEindices))

# 10. Construct a full DEresults table with consensus sign
DEresults <- data.frame(
  gene = pval_dfs_common[[1]]$gene,
  DE.fishercomb = ifelse(fishcomb$adjpval <= 0.05, 1, 0),
  DE.invnorm = ifelse(invnormcomb$adjpval <= 0.05, 1, 0),
  signFC = commonsgnFC,
  signFC_12 = commonsgn_12,
  signFC_18 = commonsgn_18
)

# 11. Add the FC matrix and logical DE (either method) only for meta-DE genes
DE_subset <- DEresults[unionDE, ]
FC_subset <- signsFC[unionDE, ]

# Use commonsgn_20 for genes consistent in ≥20/24 datasets
sign_subset <- commonsgn_18[unionDE]

# 12. Combine all info into one data.frame
FC.selecDE <- data.frame(DE_subset, FC_subset, signFC = sign_subset)

# 13. Filter for consistent DE in at least 18/24 datasets (75%)
# Keep genes where signFC is +1 or -1 
keepDE <- FC.selecDE[abs(FC.selecDE$signFC_18) == 1, ]

# Conflict = genes with mixed directions (signFC == 0)
conflictDE <- FC.selecDE[FC.selecDE$signFC_18 == 0, ]

# 14. Check dimensions and preview
dim(FC.selecDE)   # All meta-DE genes
dim(keepDE)       # Consistent DE genes (≥18/24 same sign)
dim(conflictDE)   # Inconsistent DE
head(keepDE)


fishcomb_de <- rownames(keepDE)[keepDE$DE.fishercomb == 1]
invnorm_de <- rownames(keepDE)[keepDE$DE.invnorm == 1]

# List of study columns (assuming all except the meta-analysis and signFC columns)
study_cols <- names(pval_dfs_common)

# Extract DE genes per study as named list
indstudy_de <- lapply(study_cols, function(col) {
  rownames(keepDE)[keepDE[[col]] == 1]
})
names(indstudy_de) <- study_cols

# Check the structure
str(indstudy_de)

IDD.IRR(fishcomb_de, indstudy_de)

# 11 genes detected were also detected as DE in the individual studies 

################################################################################
# 4. Heatmap showing logFC of common DEGs
################################################################################

# 1. Create a list of data frames with gene and renamed log2FoldChange
logfc_list <- lapply(names(deg_short), function(ds_name) {
  df <- deg_short[[ds_name]][, c("gene", "log2FoldChange")]
  colnames(df)[2] <- ds_name
  return(df)
})

# 2. Merge all data frames by gene using full joins
logfc_df <- Reduce(function(x, y) full_join(x, y, by = "gene"), logfc_list)

# Set gene names as rownames 
rownames(logfc_df) <- logfc_df$gene
logfc_df$gene <- NULL  

# First, make sure the genes in keepDE are in the desired order
ordered_genes <- rownames(keepDE)[order(keepDE$signFC.1)]

# Now subset and order logfc_df by those genes
subset_logFC_matrix <- logfc_df[ordered_genes, , drop = FALSE]


colnames(subset_logFC_matrix) <- ifelse(
  colnames(subset_logFC_matrix) %in% names(dataset.aliases),
  dataset.aliases[colnames(subset_logFC_matrix)],
  colnames(subset_logFC_matrix)
)

subset_logFC_matrix <- as.matrix(subset_logFC_matrix)
subset_logFC_matrix <- subset_logFC_matrix[, desired_alias_order]



# 1. Get symmetric range around 0
range_vals <- range(subset_logFC_matrix, na.rm = TRUE)
max_abs <- max(abs(range_vals))

# 2. Create separate palettes
neg_palette <- colorRampPalette(c("darkslateblue", "#A9CCE3"))(100) 
pos_palette <- colorRampPalette(c("#F5B7B1", "firebrick3"))(100) 


# 3. Combine into full palette
my_palette <- c(neg_palette, pos_palette)

# 4. Create breaks to match the color palette
breaks <- seq(-max_abs, max_abs, length.out = length(my_palette) + 1)


pdf(paste0(Sys.Date(), "_pvaluecombo_FDR5_heatmap.pdf"), width = 10, height = 5 )
# 5. Plot without scaling
pheatmap(subset_logFC_matrix,
         color = my_palette,
         breaks = breaks,
         na_col = "grey",
         cluster_cols = F,
         cluster_rows = F,
         border_color = "grey30")

dev.off()


################################################################################
# 5. Over-representation analysis with common macrophage aging genes
################################################################################

### ORA

core_macrophage_genes <- rownames(subset_logFC_matrix)
writeLines(core_macrophage_genes, "core_macrophage_genes.txt")

universe_genes <- common_genes

# Convert to Entrez IDs
universe_entrez <- bitr(universe_genes, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Convert your significant genes too
sig_genes_entrez <- bitr(core_macrophage_genes, fromType = "SYMBOL",
                         toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Warning message:
# In bitr(core_macrophage_genes, fromType = "SYMBOL", toType = "ENTREZID",  :
#           5.71% of input gene IDs are fail to map...

# Run ORA with custom universe
ego_mf <- enrichGO(gene         = sig_genes_entrez$ENTREZID,
                   universe     = universe_entrez$ENTREZID, 
                   OrgDb        = org.Mm.eg.db,
                   ont          = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)


pdf(paste0(Sys.Date(), "_pvaluecombo_FDR5_GOenrich.pdf"), width = 10, height = 5 )
dotplot(ego_mf, showCategory = 20) + ggtitle("GO ALL Enrichment")
dev.off()


pdf(paste0(Sys.Date(), "_pvaluecombo_FDR5_GOenrich.pdf"), width = 10, height = 5)

# Extract top 15 categories
df <- as.data.frame(ego_mf@result)
df <- df[df$p.adjust < 0.05, ]
df <- df[order(df$p.adjust), ][1:min(15, nrow(df)), ]

# Convert GeneRatio to numeric
if (is.character(df$GeneRatio)) {
  df$GeneRatio <- sapply(df$GeneRatio, function(x) eval(parse(text = x)))
}

# Wrap Description names
df$Description <- str_wrap(df$Description, width = 45)

# Reorder categories by GeneRatio 
df$Description <- factor(df$Description, levels = df$Description[order(df$GeneRatio)])

# Compute -log10(p.adjust)
df$neglog10 <- -log10(df$p.adjust)

# Define gradient colors
base_color <- "#7F69C9"      
light_color <- lighten(base_color, 0.7)  

color_limits <- c(1, 3)                    
color_breaks <- c(1, 1.5, 2, 2.5, 3)            

ggplot(df, aes(x = 1, y = Description)) +
  geom_point(aes(size = Count, color = neglog10)) +
  scale_color_gradient(
    low  = light_color,
    high = base_color,
    name = expression(-log[10]("FDR")),
    limits = color_limits,
    breaks = color_breaks
  ) +
  scale_size_continuous(
    range = c(0, 10),
    breaks = c(2, 4, 6),
    limits = c(0, 8),
    name = "Gene Count"
  ) +
  labs(title = "GO Enrichment", x = NULL, y = NULL) +
  guides(
    color = guide_colorbar(order = 1),
    size  = guide_legend(order = 2)
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.key.height = unit(0.7, "cm"),
    legend.key.width  = unit(1.2, "cm"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "#D3D3D3"),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 16)
  )

dev.off()

# Save the ORA results as a table

# Convert enrichGO results to a data frame
df_ego <- as.data.frame(ego_mf@result)

# filter for significant terms only
df_ego_sig <- df_ego[df_ego$p.adjust < 0.05, ]  

write.table(df_ego_sig, 
            file = file.path(my.out.dir, paste0(Sys.Date(), "_core_aging_genes_ORA.txt")), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)


#######################

sink(file = file.path(my.out.dir, paste0(Sys.Date(), "_R_session_Info_p_value_combo.txt")))
sessionInfo()
sink()



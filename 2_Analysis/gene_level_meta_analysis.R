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
library(ReactomePA)   




my.out.dir <- "/Users/ellaschwab/Benayoun_Lee_Local/ATOM/PValueCombo"

#put into list
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
pass_qc <- readLines("2026-04-19_datasets_pass_qc.txt")
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
# 2a. Assess the distribution of raw p-values
################################################################################

setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/PValueCombo")


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

# Loop through each dataset and plot the distribution of p-values
for (name in names(deg_datasets)) {
  df <- deg_datasets[[name]]
  
  if ("pvalue" %in% names(df)) {
    hist(df$pvalue,
         breaks = 100,
         col = "grey",
         main = name,
         xlab = "raw p-value distribution")
  } else {
    plot.new()
    title(main = paste(name, "\n(no 'pvalue')"))
  }
}



################################################################################
# 2b. Set raw p-values to NA where padj is NA
################################################################################

filtered_pval_list <- lapply(names(deg_datasets), function(name) {
  df <- deg_datasets[[name]]
  
  pcol <- if ("pvalue" %in% names(df)) {
    "pvalue"
  } else if ("P.Value" %in% names(df)) {
    "P.Value"
  } else {
    return(NULL)
  }
  
  p <- df[[pcol]]
  
  if ("padj" %in% names(df)) {
    p[is.na(df$padj)] <- NA
  }
  
  p
})

names(filtered_pval_list) <- names(deg_datasets)   


################################################################################
# 2c. Plot filtered raw p-value distributions
################################################################################

par(mar = c(3, 3, 3, 1))
par(mfrow = c(plot_rows, plot_cols))

for (name in names(deg_datasets)) {
  p <- filtered_pval_list[[name]]
  
  if (!is.null(p)) {
    hist(p,
         breaks = 100,
         col = "grey",
         main = name,
         xlab = "filtered raw p-value distribution")
  } else {
    plot.new()
    title(main = paste(name, "\n(no p-value column)"))
  }
}

par(mfrow = c(1, 1))

################################################################################
# 3. P-value combination
################################################################################

# 1. Make a short dataframe for each study with gene, log2FC, and filtered p-values
deg_short <- lapply(names(deg_datasets), function(ds_name) {
  df <- deg_datasets[[ds_name]]
  
  pcol <- if ("pvalue" %in% names(df)) {
    "pvalue"
  } else if ("P.Value" %in% names(df)) {
    "P.Value"
  } else {
    stop("No p-value column found")
  }
  
  p <- df[[pcol]]
  if ("padj" %in% names(df)) {
    p[is.na(df$padj)] <- NA
  }
  
  data.frame(
    gene = if ("row.names" %in% names(df)) df[["row.names"]] else rownames(df),
    log2FoldChange = if ("log2FoldChange" %in% names(df)) df$log2FoldChange else NA,
    pvalue = p,
    adjpvalue = df$padj
  )
})

names(deg_short) <- names(deg_datasets)


# Remove rows with NA p-values or logFCs
pval_dfs_clean <- lapply(deg_short, function(df) {
  df %>% filter(!is.na(pvalue), !is.na(log2FoldChange))
})

# Remove duplicated genes (keep only the first)
pval_dfs_unique <- lapply(pval_dfs_clean, function(df) {
  df %>% distinct(gene, .keep_all = TRUE)
})


# Get the set of common genes across all studies
common_genes <- sort(Reduce(intersect, lapply(pval_dfs_unique, function(df) df$gene)))
length(common_genes)

# [1] 9225


# output the universe to a text file
write.table(common_genes, 
            file = "common_genes_universe.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)


# Filter all datasets to the common genes and order alphabetically
pval_dfs_common <- lapply(pval_dfs_unique, function(df) {
  df %>%
    filter(gene %in% common_genes) %>%
    arrange(match(gene, common_genes))
})

# Create raw p-value list
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
fc_list <- lapply(pval_dfs_common, function(df) {
  df$log2FoldChange
})


signsFC <- do.call(cbind, lapply(fc_list, sign))

# Assign correct row names
rownames(signsFC) <- common_genes


# Count +1 per gene
pos_count <- rowSums(signsFC == 1)

# Count -1 per gene
neg_count <- rowSums(signsFC == -1)

# Combine into a table
sign_counts <- data.frame(
  gene = rownames(signsFC),
  pos_count = pos_count,
  neg_count = neg_count
)

# Filter for genes that are the same sign in at >= 18/24 datasets
genes_18_of_24 <- sign_counts %>%
  filter(pos_count >= 18 | neg_count >= 18)


# Identify all DE genes from Fisher or Inverse Normal meta-analysis
unionDE <- unique(c(fishcomb$DEindices, invnormcomb$DEindices))

#  Construct a full DEresults table with consensus sign
DEresults <- data.frame(
  gene = common_genes,
  DE.fishercomb = ifelse(fishcomb$adjpval <= 0.05, 1, 0),
  DE.invnorm = ifelse(invnormcomb$adjpval <= 0.05, 1, 0)
)

DEresults_filtered <- DEresults %>%
  filter(DE.fishercomb == 1 | DE.invnorm == 1) %>%   # keep sig in either
  inner_join(genes_18_of_24, by = "gene")           # then apply sign filter

# Plot sign-consistency among genes significant in both meta-analysis methods

# Get genes significant in either method
DE_genes <- DEresults$gene[
  DEresults$DE.fishercomb == 1 | DEresults$DE.invnorm == 1
]

# Filter sign_counts to only those genes
sign_counts_DE <- sign_counts %>%
  filter(gene %in% DE_genes)

# Recompute max agreement
sign_counts_DE <- sign_counts_DE %>%
  mutate(max_same_sign = pmax(pos_count, neg_count))

# Threshold counts
thresholds <- 12:24

gene_counts <- sapply(thresholds, function(t) {
  sum(sign_counts_DE$max_same_sign >= t)
})

plot_df <- data.frame(
  threshold = thresholds,
  gene_count = gene_counts
)

write.csv(plot_df, paste0(Sys.Date(), "_significant_genes_sensitivity_table.csv"), row.names = FALSE)

pdf(paste0(Sys.Date(), "_significant_genes_sensitivity_barplot.pdf"), width = 8, height = 8 )

library(ggplot2)

ggplot(plot_df, aes(x = factor(threshold), y = gene_count)) +
  geom_col(fill = "steelblue") + 
  geom_text(aes(label = gene_count), vjust = -0.3, size = 4) +
  scale_y_continuous(
    limits = c(0, 9000),
    breaks = seq(0, 9000, by = 1000),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Minimum number of datasets with same sign",
    y = "Number of DE genes",
    title = "Sign consistency among DE genes"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),   # remove grid
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black")
  )

dev.off()

FC.selecDE <- DEresults %>%
  filter(DE.fishercomb == 1 | DE.invnorm == 1) %>%   # meta-DE genes in either method
  left_join(sign_counts, by = "gene") %>%          
  mutate(
    signFC = case_when(
      pos_count >= 18 ~ 1,
      neg_count >= 18 ~ -1,
      TRUE ~ 0
    )
  )

keepDE <- FC.selecDE %>% filter(signFC != 0)
conflictDE <- FC.selecDE %>% filter(signFC == 0)


# Check dimensions and preview
dim(FC.selecDE)   # All meta-DE genes
dim(keepDE)       # Consistent DE genes (≥18/24 same sign)
dim(conflictDE)   # Inconsistent DE
head(keepDE)



################################################################################
# 4. Heatmap showing logFC of common DEGs
################################################################################

# Build a merged log2FC table
logfc_list <- lapply(names(deg_short), function(ds_name) {
  df <- deg_short[[ds_name]][, c("gene", "log2FoldChange")]
  colnames(df)[2] <- ds_name
  df
})

logfc_df <- Reduce(function(x, y) full_join(x, y, by = "gene"), logfc_list)

rownames(logfc_df) <- logfc_df$gene
logfc_df$gene <- NULL




################################################################################
# 593 Gene Heatmap 
################################################################################

keepDE_ordered_by_consistency <- keepDE %>%
  mutate(max_same_sign = pmax(pos_count, neg_count)) %>%
  arrange(signFC, desc(max_same_sign), gene)

ordered_genes <- keepDE_ordered_by_consistency$gene

# Subset and order
subset_logFC_matrix <- logfc_df[ordered_genes, , drop = FALSE]

# Rename columns
new_colnames <- colnames(subset_logFC_matrix)
idx <- match(new_colnames, names(dataset.aliases))
new_colnames[!is.na(idx)] <- unname(unlist(dataset.aliases)[idx[!is.na(idx)]])
colnames(subset_logFC_matrix) <- new_colnames

# Reorder columns
subset_logFC_matrix <- as.matrix(subset_logFC_matrix)
subset_logFC_matrix <- subset_logFC_matrix[, desired_alias_order, drop = FALSE]

# Symmetric color scale
range_vals <- range(subset_logFC_matrix, na.rm = TRUE)
max_abs <- max(abs(range_vals))
neg_palette <- colorRampPalette(c("darkslateblue", "#A9CCE3"))(100)
pos_palette <- colorRampPalette(c("#F5B7B1", "firebrick3"))(100)
my_palette <- c(neg_palette, pos_palette)
breaks <- seq(-0.5, 0.5, length.out = length(my_palette) + 1)



pdf(paste0(Sys.Date(), "_pvaluecombo_ALL_FDR5_heatmap.pdf"), width = 30, height = 10)
pheatmap(
  subset_logFC_matrix,
  color         = my_palette,
  breaks        = breaks,
  legend_breaks  = seq(-0.4, 0.4, by = 0.2),
  legend_labels  = seq(-0.4, 0.4, by = 0.2),
  na_col        = "grey",
  cluster_cols  = FALSE,
  cluster_rows  = FALSE,
  border_color  = NA,            
  show_rownames = FALSE          # hides gene names
)
dev.off()

################################################################################
# Supplementary: 35 consistent genes (>=21/24) heatmap, split UP vs DOWN
################################################################################

# Filter to highly consistent genes, split by direction
consistent_up <- keepDE_ordered_by_consistency %>%
  filter(max_same_sign >= 21, signFC == 1) %>%
  arrange(desc(max_same_sign), gene)

consistent_down <- keepDE_ordered_by_consistency %>%
  filter(max_same_sign >= 21, signFC == -1) %>%
  arrange(desc(max_same_sign), gene)

# Combine: UP first, then DOWN
ordered_genes <- c(consistent_up$gene, consistent_down$gene)

# Subset
subset_logFC_matrix <- logfc_df[ordered_genes, , drop = FALSE]

# Rename columns via alias map
new_colnames <- colnames(subset_logFC_matrix)
idx <- match(new_colnames, names(dataset.aliases))
new_colnames[!is.na(idx)] <- unname(unlist(dataset.aliases)[idx[!is.na(idx)]])
colnames(subset_logFC_matrix) <- new_colnames

# Reorder columns
subset_logFC_matrix <- as.matrix(subset_logFC_matrix)
subset_logFC_matrix <- subset_logFC_matrix[, desired_alias_order, drop = FALSE]


# Dynamic height based on number of genes
n_genes <- length(ordered_genes)
pdf_height <- max(8, n_genes * 0.18)

pdf(paste0(Sys.Date(), "_pvaluecombo_CONSISTENT_21of24_UPDOWN_heatmap.pdf"),
    width = 10, height = pdf_height)
pheatmap(
  subset_logFC_matrix,
  color             = my_palette,
  breaks            = breaks,
  legend_breaks     = seq(-0.4, 0.4, by = 0.2),
  legend_labels     = seq(-0.4, 0.4, by = 0.2),
  na_col            = "grey",
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  border_color      = "grey30",
  show_rownames     = TRUE,
  fontsize_row      = 7,
  gaps_row          = nrow(consistent_up)
)
dev.off()


################################################################################
# 5. Over-representation analysis (GO, Reactome, KEGG)
################################################################################

universe_genes  <- common_genes
universe_entrez <- bitr(universe_genes, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Mm.eg.db)


run_ora_plot <- function(gene_symbols, direction, high_color, out_dir = my.out.dir) {
  
  writeLines(gene_symbols,
             file.path(out_dir, paste0("core_macrophage_genes_", direction, ".txt")))
  
  sig_entrez <- bitr(gene_symbols, fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  # ----- Run enrichment against each database ---------------------------------
  enrich_results <- list(
    GO = enrichGO(gene          = sig_entrez$ENTREZID,
                  universe      = universe_entrez$ENTREZID,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1,         # collect all results so we can show what didn't pass significance
                  qvalueCutoff  = 1,        
                  readable      = TRUE),
    
    Reactome = enrichPathway(gene          = sig_entrez$ENTREZID,
                             universe      = universe_entrez$ENTREZID,
                             organism      = "mouse",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 1,    
                             qvalueCutoff  = 1,    
                             readable      = TRUE),
    
    KEGG = enrichKEGG(gene          = sig_entrez$ENTREZID,
                      universe      = universe_entrez$ENTREZID,
                      organism      = "mmu",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,           
                      qvalueCutoff  = 1)          
  )
  
  # KEGG doesn't have a readable option, so convert ENTREZ IDs to symbols manually
  if (!is.null(enrich_results$KEGG) && nrow(as.data.frame(enrich_results$KEGG)) > 0) {
    enrich_results$KEGG <- setReadable(enrich_results$KEGG,
                                       OrgDb   = org.Mm.eg.db,
                                       keyType = "ENTREZID")
  }
  
  # ----- Loop over databases: save tables and make plots ----------------------
  for (db in names(enrich_results)) {
    ego <- enrich_results[[db]]
    
    if (is.null(ego) || nrow(as.data.frame(ego@result)) == 0) {
      message("No results returned for ", db, " (", direction, ")")
      next
    }
    
    # Save full results table (all terms, including non-significant)
    df_full <- as.data.frame(ego@result)
    write.table(df_full,
                file = file.path(out_dir,
                                 paste0(Sys.Date(), "_core_aging_genes_ORA_",
                                        db, "_", direction, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Filter for plotting
    df <- df_full[df_full$p.adjust < 0.10, ] # fix
    
    if (nrow(df) == 0) {
      message("No significant ", db, " terms for ", direction)
      next
    }
    
    df <- df[order(df$p.adjust), ][1:min(10, nrow(df)), ]
    
    if (is.character(df$GeneRatio)) {
      df$GeneRatio <- sapply(df$GeneRatio, function(x) eval(parse(text = x)))
    }
    df$Description <- str_wrap(df$Description, width = 45)
    df$Description <- factor(df$Description, levels = df$Description[order(df$GeneRatio)])
    df$neglog10    <- -log10(df$p.adjust)
    
    light_color   <- lighten(high_color, 0.7)
    min_count     <- min(df$Count,    na.rm = TRUE)
    max_count     <- max(df$Count,    na.rm = TRUE)
    size_breaks   <- pretty(c(min_count, max_count), n = 4)
    color_limits <- c(1.0, 5.5) 
    color_breaks <- 1:6
    
    pdf(file.path(out_dir,
                  paste0(Sys.Date(), "_pvaluecombo_FDR5_",
                         db, "enrich_", direction, ".pdf")),
        width = 10, height = 12)
    
    p <- ggplot(df, aes(x = 1, y = Description)) +
      geom_point(aes(size = Count, color = neglog10)) +
      scale_color_gradient(
        low    = light_color,
        high   = high_color,
        name   = expression(-log[10]("FDR")),
        limits = color_limits,
        breaks = color_breaks
      ) +
      scale_size_continuous(
        range  = c(2, 8),
        breaks = size_breaks,
        limits = c(min_count, max_count),
        name   = "Gene Count"
      ) +
      labs(title = paste0(db, " Enrichment — ", direction, "-regulated"),
           x = NULL, y = NULL) +
      guides(
        color = guide_colorbar(order = 1),
        size  = guide_legend(order = 2)
      ) +
      theme(
        axis.title.x       = element_blank(),
        axis.text.x        = element_blank(),
        axis.ticks.x       = element_blank(),
        axis.text.y        = element_text(size = 14),
        legend.key.height  = unit(0.7, "cm"),
        legend.key.width   = unit(1.2, "cm"),
        panel.border       = element_rect(color = "black", fill = NA, size = 1),
        panel.background   = element_rect(fill = "white", color = NA),
        panel.grid.major.y = element_line(color = "#D3D3D3"),
        plot.background    = element_rect(fill = "white", color = NA),
        plot.title         = element_text(face = "bold", size = 16)
      )
    print(p)
    dev.off()
  }
  
  invisible(enrich_results)
}

# Split genes by direction using the consistency counts
up_genes   <- keepDE_ordered_by_consistency$gene[keepDE_ordered_by_consistency$pos_count >= 18]
down_genes <- keepDE_ordered_by_consistency$gene[keepDE_ordered_by_consistency$neg_count >= 18]

cat("Up genes  (pos_count >= 18):", length(up_genes),   "\n")
cat("Down genes (neg_count >= 18):", length(down_genes), "\n")

ora_up   <- run_ora_plot(up_genes,   direction = "up",   high_color = "firebrick3")
ora_down <- run_ora_plot(down_genes, direction = "down", high_color = "darkslateblue")



#######################

sink(file = file.path(my.out.dir, paste0(Sys.Date(), "_R_session_Info_p_value_combo.txt")))
sessionInfo()
sink()






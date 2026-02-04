####################################################################################
# Aging transcriptomics of mouse macrophages 
# Strip plot 
####################################################################################

library(stringr)
set.seed(1234) # set seed for reproducibility


################################################################################
# 1. Read in the data
################################################################################

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


alias.to.dataset <- setNames(names(dataset.aliases), unlist(dataset.aliases))
ordered_dataset_names <- alias.to.dataset[desired_alias_order]



################################################################################
# 3. Create the jitterplot
################################################################################

# code adapted from https://github.com/brunetlab/Leeman_et_al_2017/blob/master/kallisto_deseq2/Fig4A_stripplot_cell_type_colors.R

setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/StripPlot")

# Get corresponding sex colors
ct.unique <- data.frame("label" = names(my_deseq_list)) 
ct.unique$cols <- c(rep("deepskyblue", 17), rep("deeppink", 7))

# Sort genes in each dataset by p-value
sex.results <- lapply(my_deseq_list, function(x) x[order(x$padj),])
names(sex.results) <- names(my_deseq_list)

# Reorder datasets by niche
sex.results <- sex.results[ordered_dataset_names]

# Assign colors based on significance and direction
cols <- list()
xlab <- character(length(sex.results))

for (i in seq_along(sex.results)) {
  df <- sex.results[[i]]
  sig.idx <- which(df$padj < 0.05)
  
  # Start with light gray for all
  col.vec <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), nrow(df))
  
  # Red for significant up
  col.vec[sig.idx[df$log2FoldChange[sig.idx] > 0]] <- "firebrick3"
  
  # Blue for significant down
  col.vec[sig.idx[df$log2FoldChange[sig.idx] < 0]] <- "darkslateblue"
  
  cols[[i]] <- col.vec
  
  # Label with alias and number of significant genes
  alias <- dataset.aliases[[names(sex.results)[i]]]
  xlab[i] <- paste0(alias, "\n(", length(sig.idx), ")")
}

names(cols) <- names(sex.results)

# Begin plotting
pdf(paste0(Sys.Date(), "_stripplot_DESeq2_full_length.pdf"), width = 16, height = 6)

# Margins
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))

# Empty canvas
plot(x = 1, y = 1, type = "n",
     xlim = c(0.5, length(sex.results) + 0.5),
     ylim = c(-1.5, 1.5),
     axes = FALSE,
     xlab = "",
     ylab = "Log2FC")

# Get coordinates of the plot area
plot_coords <- par("usr")

# Highlight up/down zones
polygon(x = c(plot_coords[1], plot_coords[2], plot_coords[2], plot_coords[1]),
        y = c(0, 0, plot_coords[4], plot_coords[4]),
        col = "#CD262633", border = NA)

polygon(x = c(plot_coords[1], plot_coords[2], plot_coords[2], plot_coords[1]),
        y = c(0, 0, plot_coords[3], plot_coords[3]),
        col = "#483D8B33", border = NA)

# Add horizontal lines
abline(h = 0)
abline(h = seq(-2, 2, by = 0.5)[-3], lty = "dotted", col = "grey89")

# Plot points per dataset
for (i in seq_along(sex.results)) {
  set.seed(1234)
  df <- sex.results[[i]]
  points(x = jitter(rep(i, nrow(df)), amount = 0.2),
         y = rev(df$log2FoldChange),
         pch = 16,
         col = rev(cols[[i]]),
         bg  = rev(cols[[i]]),
         cex = 0.4)
}

# Add axes
axis(1, at = 1:length(sex.results), tick = FALSE, las = 2,
     lwd = 0, labels = xlab, cex.axis = 0.7)

axis(2, las = 1, at = seq(-1.5, 1.5, 0.5))

box()
dev.off()



###############################################################################################

# Separate male and female dataset names based on your new_order vector
male.datasets <- c(
  "GSE93202_Spleen", "GSE93202_VAT", "GSE98401_Microglia", "GSE134397_Alveolar_CTL",
  "GSE131869_Microglia_M", "GSE128830_Peritoneal", "GSE137028_Microglia",
  "GSE154832_eWAT", "PRJNA682234_Callus", "GSE267529_Microglia_M",
  "GSE145295_Alveolar", "GSE142580_SkM", "GSE190689_Alveolar", 
  "PRJNA800823_SkM", "PRJNA1173774_BMDM", "PRJNA524906_Microglia", 
  "PRJNA816431_Callus"
)

female.datasets <- c(
  "GSE199763_SkinWound", "GSE199879_Spleen_Red_Pulp", "GSE156762_Microglia",
  "GSE267529_Microglia_F", "GSE131869_Microglia_F",
  "BMDM_NIA_F", "GSE134397_Alveolar"
)

# Subset the sex.results and cols lists for males and females
male_results <- sex.results[male.datasets]
female_results <- sex.results[female.datasets]

male_cols <- cols[male.datasets]
female_cols <- cols[female.datasets]

male_n <- sapply(male_results, nrow)
female_n <- sapply(female_results, nrow)

male_xlab <- character(length = length(male_results))
for(i in seq_along(male_results)){
  ind.sig.i <- male_results[[i]]$padj < 0.05
  male_xlab[i] <- paste(names(male_results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}

female_xlab <- character(length = length(female_results))
for(i in seq_along(female_results)){
  ind.sig.i <- female_results[[i]]$padj < 0.05
  female_xlab[i] <- paste(names(female_results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}

pdf(paste0(Sys.Date(),"_stripplot_DESeq2_sex_split.pdf"), width = 16, height = 6)

par(mfrow = c(1, 2))  # side by side
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))

# -------- Male plot --------
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, length(male_results) + 0.5),
     ylim = c(-1.5, 1.5),
     axes = FALSE,
     xlab = "",
     ylab = "Log2FC"
)

x_range <- c(0.5, length(male_results) + 0.5)
y_range <- c(-1.5, 1.5)

polygon(
  c(x_range[1], x_range[2], x_range[2], x_range[1]),
  c(0, 0, y_range[2], y_range[2]),
  col = "#CD26264D", border = NA
)
polygon(
  c(x_range[1], x_range[2], x_range[2], x_range[1]),
  c(0, 0, y_range[1], y_range[1]),
  col = "#483D8B4D", border = NA 
)

abline(h = 0)
abline(h = seq(-1.5, 1.5, by = 0.5)[-which(seq(-1.5, 1.5, 0.5) == 0)], lty = "dotted", col = "grey")

for(i in seq_along(male_results)){
  set.seed(1234)
  points(x = jitter(rep(i, male_n[i]), amount = 0.2),
         y = rev(male_results[[i]]$log2FoldChange),
         pch = 16,
         col = rev(male_cols[[i]]),
         bg = rev(male_cols[[i]]),
         cex = 0.4)
}
axis(1, at = seq_along(male_results), tick = FALSE, las = 2, lwd = 0, labels = male_xlab, cex.axis = 0.7)
axis(2, las = 1, at = seq(-1.5, 1.5, 0.5))
box()
title(main = "Male Datasets")


# -------- Female plot --------
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, length(female_results) + 0.5),
     ylim = c(-1.5, 1.5), 
     axes = FALSE,
     xlab = "",
     ylab = ""
)

x_range <- c(0.5, length(female_results) + 0.5)
y_range <- c(-1.5, 1.5)  

polygon(
  c(x_range[1], x_range[2], x_range[2], x_range[1]),
  c(0, 0, y_range[2], y_range[2]),
  col = "#CD26264D", border = NA
)
polygon(
  c(x_range[1], x_range[2], x_range[2], x_range[1]),
  c(0, 0, y_range[1], y_range[1]),
  col = "#483D8B4D", border = NA 
)

abline(h = 0)
abline(h = seq(-1.5, 1.5, by = 0.5)[-which(seq(-1.5, 1.5, 0.5) == 0)], lty = "dotted", col = "grey")

for(i in seq_along(female_results)){
  set.seed(1234)
  points(x = jitter(rep(i, female_n[i]), amount = 0.2),
         y = rev(female_results[[i]]$log2FoldChange),
         pch = 16,
         col = rev(female_cols[[i]]),
         bg = rev(female_cols[[i]]),
         cex = 0.4)
}
axis(1, at = seq_along(female_results), tick = FALSE, las = 2, lwd = 0, labels = female_xlab, cex.axis = 0.7)
axis(2, las = 1, at = seq(-1.5, 1.5, 0.5))  
box()
title(main = "Female Datasets")


dev.off()


###############################################################################################

sink(file = paste(Sys.Date(),"ATOM_stripplot_RNAseq_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()


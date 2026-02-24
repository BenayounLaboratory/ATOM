####################################################################################
# Aging transcriptomics of mouse macrophages 
# Pre-processing - check Xist/Ddx3y, normalize, read depth check, replicates, DEGs
####################################################################################

set.seed(1234) # set seed for reproducibility

library('sva')
library('limma')
library(DESeq2)
library(ggplot2)
library(dplyr)


setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/final_filtered_counts")

################################################################################
# 1. Read in the data
################################################################################


#Males
GSE128830_Peritoneal <- read.csv("2025-03-07_GSE128830_Peritoneal_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE205395_Skeletal <- read.csv("2025-03-07_GSE205395_Skeletal_Muscle_Macrophages_aging_Male_counts.txt", sep = "\t", header = T, skip = 1) #this is actually a F set 
GSE93202_Spleen <- read.csv("2021-06-04_GSE93202_Spleen_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE93202_VAT <- read.csv("2021-06-04_GSE93202_VAT_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE98249_BAM <- read.csv("2021-06-04_GSE98249_BAM_macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE98249_BMM <- read.csv("2021-06-04_GSE98249_BMM_macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE98401_Brain <- read.csv("2021-06-04_GSE98401_Brain_microglia_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE134397_Alveolar_CTL <- read.csv("2021-06-04_GSE134397_Alveolar_Macrophages_CTL_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE137028_Microglia <- read.csv("2021-06-04_GSE137028_Microglia_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE154832_eWAT <- read.csv("2021-06-04_GSE154832_eWAT_Phagocytic_SVF_aging_counts.txt", sep = "\t", header = T, skip = 1)
PRJNA682234_Callus <- read.csv("2021-06-04_PRJNA682234_Callus_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)

PRJNA800823_SkM <- read.csv("2025-08-07_PRJNA800823_Skeletal_Muscle_macrophages_aging_counts.txt",  sep = "\t", header = T)
PRJNA1173774_BMDM <- read.csv("2025-12-11_PRJNA1173774_BMDM_aging_counts_mm10.txt", sep = "\t", header = T)

#Females
GSE199763_SkinWound <- read.csv("2025-03-07_GSE199763_SkinWound_Macrophages_aging_Female_counts.txt", sep = "\t", header = T, skip = 1)
GSE199879_Spleen_Red_Pulp <- read.csv("2025-03-07_GSE199879_Spleen_Red_Pulp_Macrophages_aging_Female_counts.txt", sep = "\t", header = T, skip = 1)
GSE132882_Nerve <- read.csv("2021-06-04_GSE132882_Nerve_Macrophage_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE156762_Microglia <- read.csv("2021-08-06_GSE156762_Microglia_JCI_counts.txt", sep = "\t", header = T, skip = 1)

#MixedSex
GSE124829_ImmGen_Peritoneal <- read.csv("2021-06-04_GSE124829_ImmGen_Peritoneal_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE131869_Microglia <- read.csv("2021-06-04_GSE131869_Microglia_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE267529_Microglia <- read.csv("2021-07-27_GSE267529_Microglia_Goodridge_counts.txt", sep = "\t", header = T, skip = 1)
BMDM_NIA <- read.csv("2024-03-28_BMDM_NIA_aging_CLEAN_counts.txt", sep = "\t", header = T, skip = 1)

#Unknown Sex
GSE190689_Alveolar <- read.csv("2025-03-07_GSE190689_Alveolar_Macrophages_aging_UnkSex_counts.txt", sep = "\t", header = T, skip = 1)
GSE124872_Alveolar <- read.csv("2021-06-04_GSE124872_Alveolar_Macrophages_aging_counts.txt" , sep = "\t" , header = T , skip = 1)
GSE142580_SkM <- read.csv("2021-06-04_GSE142580_SkM_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE145295_Alveolar <- read.csv("2021-06-04_GSE145295_Alveolar_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE134397_Alveolar <- read.csv("2025-03-07_GSE134397_Alveolar_Macrophages_aging_Mixed_Sex_counts.txt", sep = "\t", header = T, skip = 1)

PRJNA847895_BMDM <- read.csv("2025-12-11_PRJNA847895_BMDM_aging_counts_mm10.txt", sep = "\t", header = T)
PRJNA524906_microglia <- read.csv("2025-12-11_PRJNA524906_microglia_aging_counts_mm10.txt", sep = "\t", header = T)
PRJNA816431_callus <- read.csv("2025-12-11_PRJNA816431_callus_aging_counts_mm10.txt", sep = "\t", header = T)



################################################################################
# 2. Perform visual checks on unknown datasets to sort into sex
################################################################################


print(GSE190689_Alveolar[GSE190689_Alveolar$Geneid == "Xist",-c(2:6)])
print(GSE190689_Alveolar[GSE190689_Alveolar$Geneid == "Ddx3y",-c(2:6)])
# one individual has unclear sex due to low Ddx3y, remove
GSE190689_Alveolar <- GSE190689_Alveolar[,-9]
# the rest are male, file with males

print(GSE124872_Alveolar[GSE124872_Alveolar$Geneid == "Xist",-c(2:6)])
print(GSE124872_Alveolar[GSE124872_Alveolar$Geneid == "Ddx3y",-c(2:6)])
# one individual is a female, the rest are males, remove female and file with males
GSE124872_Alveolar <- GSE124872_Alveolar[,-11]

print(GSE142580_SkM[GSE142580_SkM$Geneid == "Xist",-c(2:6)])
print(GSE142580_SkM[GSE142580_SkM$Geneid == "Ddx3y",-c(2:6)])
# all are male, file with males

print(GSE145295_Alveolar[GSE145295_Alveolar$Geneid == "Xist",-c(2:6)])
print(GSE145295_Alveolar[GSE145295_Alveolar$Geneid == "Ddx3y",-c(2:6)])
# all are male, file with males

print(GSE134397_Alveolar[GSE134397_Alveolar$Geneid == "Xist",-c(2:6)])
print(GSE134397_Alveolar[GSE134397_Alveolar$Geneid == "Ddx3y",-c(2:6)])
# all are female, file with females

print(PRJNA847895_BMDM[PRJNA847895_BMDM$Geneid == "Xist",-c(2:6)])
print(PRJNA847895_BMDM[PRJNA847895_BMDM$Geneid == "Ddx3y",-c(2:6)])
# all are male, file with males

print(PRJNA524906_microglia[PRJNA524906_microglia$Geneid == "Xist",-c(2:6)])
print(PRJNA524906_microglia[PRJNA524906_microglia$Geneid == "Ddx3y",-c(2:6)])
# all are male, file with males

print(PRJNA816431_callus[PRJNA816431_callus$Geneid == "Xist",-c(2:6)])
print(PRJNA816431_callus[PRJNA816431_callus$Geneid == "Ddx3y",-c(2:6)])
# all are male, file with males

print(GSE205395_Skeletal[GSE205395_Skeletal$Geneid == "Xist",-c(2:6)])
print(GSE205395_Skeletal[GSE205395_Skeletal$Geneid == "Ddx3y",-c(2:6)])
# all are actually female, file with females, sex is misreported



#######################################################################
# Handle any potential mislabelling

#labelled male
print(GSE205395_Skeletal[GSE205395_Skeletal$Geneid == "Xist",-c(2:6)])
print(GSE205395_Skeletal[GSE205395_Skeletal$Geneid == "Ddx3y",-c(2:6)])
#no Ddx3y expression, high Xist expression in all samples, so these are likely female
#file with females

#labelled male
print(GSE98249_BMM[GSE98249_BMM$Geneid == "Xist",-c(2:6)]) # young males, old females 
print(GSE98249_BMM[GSE98249_BMM$Geneid == "Ddx3y",-c(2:6)])
# remove from analysis

#labelled male
print(GSE98249_BAM[GSE98249_BAM$Geneid == "Xist",-c(2:6)]) # young males, old females 
print(GSE98249_BAM[GSE98249_BAM$Geneid == "Ddx3y",-c(2:6)])
# remove from analysis

print(GSE131869_Microglia[GSE131869_Microglia$Geneid == "Xist",-c(2:6)])
print(GSE131869_Microglia[GSE131869_Microglia$Geneid == "Ddx3y",-c(2:6)])
# one individual has unclear sex due to high Ddx3y, remove
GSE131869_Microglia <- GSE131869_Microglia[
  , colnames(GSE131869_Microglia) !=
    "GSM3823691_Microglia_2m_female_STAR_Aligned.sortedByCoord.out.bam"
]

########################################################################
# Separate males and females from the mixed sex datasets

my.both.data.list <- list("GSE124829_ImmGen_Peritoneal"  = GSE124829_ImmGen_Peritoneal,
                          "GSE131869_Microglia"          = GSE131869_Microglia,
                          "GSE267529_Microglia"          = GSE267529_Microglia,
                          "BMDM_NIA"                     = BMDM_NIA)

for(i in 1:length(my.both.data.list)){
  dataset_name <- names(my.both.data.list[i])
  my.data <- my.both.data.list[[i]]
  
  #look for female columns (only works with this naming convention)
  femalecols <- grep("_female" , colnames(my.data), value = FALSE, ignore.case = TRUE)
  if(length(femalecols) == 0){ #for BMDM_NIA column name convention
    femalecols <- grep("_YF" , colnames(my.data), value = FALSE)
    oldfemalecols <- grep("_OF" , colnames(my.data), value = FALSE)
    femalecols <- append(femalecols,oldfemalecols)
  }
  #subset data
  female_data <- my.data[,c(1,2,3,4,5,6, femalecols)]
  male_data <- my.data[,-femalecols]
  
  assign(paste0(dataset_name, "_F"), female_data)
  assign(paste0(dataset_name, "_M"), male_data)
  
}

########################################################################
my.male.data.list <- list("GSE93202_Spleen"         = GSE93202_Spleen,
                          "GSE93202_VAT"            = GSE93202_VAT,
                          "GSE98401_Microglia"          = GSE98401_Brain,
                          "GSE134397_Alveolar_CTL"      = GSE134397_Alveolar_CTL,
                          "GSE124829_ImmGen_Peritoneal_M"  = GSE124829_ImmGen_Peritoneal_M,
                          "GSE131869_Microglia_M"   = GSE131869_Microglia_M,
                          "GSE128830_Peritoneal" = GSE128830_Peritoneal,
                          "GSE137028_Microglia" = GSE137028_Microglia,
                          "GSE154832_eWAT" = GSE154832_eWAT,
                          "PRJNA682234_Callus" = PRJNA682234_Callus,
                          "GSE267529_Microglia_M" = GSE267529_Microglia_M,
                          "BMDM_NIA_M" = BMDM_NIA_M, 
                          "GSE145295_Alveolar" = GSE145295_Alveolar,
                          "GSE142580_SkM" = GSE142580_SkM,
                          "GSE124872_Alveolar" = GSE124872_Alveolar,
                          "GSE190689_Alveolar" = GSE190689_Alveolar,
                          "PRJNA524906_Microglia" = PRJNA524906_microglia,
                          "PRJNA816431_Callus" = PRJNA816431_callus,
                          "PRJNA800823_SkM" = PRJNA800823_SkM,
                          "PRJNA1173774_BMDM" = PRJNA1173774_BMDM,
                          "PRJNA847895_BMDM" = PRJNA847895_BMDM
)

#"GSE98249_BAM"            = GSE98249_BAM, # removed due to old F, young M
#"GSE98249_BMM"            = GSE98249_BMM, # removed due to old F, young M

my.female.data.list <- list("GSE199763_SkinWound" = GSE199763_SkinWound,
                            "GSE199879_Spleen_Red_Pulp" = GSE199879_Spleen_Red_Pulp,
                            "GSE156762_Microglia" = GSE156762_Microglia,
                            "GSE267529_Microglia_F" = GSE267529_Microglia_F,
                            "GSE132882_Nerve"  = GSE132882_Nerve,
                            "GSE124829_ImmGen_Peritoneal_F" = GSE124829_ImmGen_Peritoneal_F,
                            "GSE131869_Microglia_F"   = GSE131869_Microglia_F,
                            "BMDM_NIA_F" = BMDM_NIA_F,
                            "GSE134397_Alveolar" = GSE134397_Alveolar,
                            "GSE205395_Skeletal"   = GSE205395_Skeletal)


# Construct the path with today's date
output_dir <- file.path("/Users/ellaschwab/Benayoun_Lee_Local/ATOM", paste0(Sys.Date(), "_continuous_age_output"))

# Create the output directory 
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set it as the working directory
setwd(output_dir)

my.percent.difs <- list()
my.degs <- list()



################################################################################
# 3. Function to pre-process the count files and perform DGE analysis 
################################################################################


# Function to process RNA-seq datasets and do differential gene expression analysis
process_rnaseq_dataset <- function(data.list, sex, outprefix.base = Sys.Date()) {
  
  # Initialize storage for output
  dataset_age_df <- data.frame(dataset = character(), age_group = numeric(), stringsAsFactors = FALSE)
  percent_difs <- list()
  deg_counts <- list()
  
  # Extract months from column names function
  extract_months <- function(colname) {
    if (grepl("_(\\d+\\.?\\d*)m_", colname)) {
      # Match _Xm_ or _X.Xm_
      match <- regmatches(colname, regexpr("_(\\d+\\.?\\d*)m_", colname))
      as.numeric(gsub("[^0-9\\.]", "", match))
    } 
    else if (grepl("_Old_", colname)) {
      25  # Set age manually for Goodridge microglia old
    }
    else if (grepl("_Young_", colname)) {
      4   # Set age manually for Goodridge microglia young
    }
    else if (grepl("_OM", colname)) {
      20  # Set age manually for BMDM NIA Old Male
    }
    else if (grepl("_YM", colname)) {
      4   # Set age manually for BMDM NIA Young Male
    }
    else if (grepl("_OF", colname)) {
      20  # Set age manually for BMDM NIA Old Female
    }
    else if (grepl("_YF", colname)) {
      4   # Set age manually for BMDM NIA Young Female
    }
    else if (grepl("_(\\d+)_months_", colname)) {
      # match pattern for _X_months_
      match <- regmatches(colname, regexpr("_(\\d+)_months_", colname))
      as.numeric(gsub("[^0-9\\.]", "", match))
    } 
    else if (grepl("_(\\d+)months_", colname)) {
      # match pattern for _Xmonths_
      match <- regmatches(colname, regexpr("_(\\d+)months_", colname))
      as.numeric(gsub("[^0-9\\.]", "", match))
    } 
    else if (grepl("_(\\d+)m(\\d+)_", colname)) {
      # match pattern for _XmX_
      match <- regmatches(colname, regexec("_(\\d+)m(\\d+)_", colname))[[1]]
      as.numeric(match[2])
    } 
    else {
      NA
    }
  }
  
  # Main processing loop
  for(i in 1:length(data.list)) {
    my.outprefix <- paste(outprefix.base, names(data.list[i]), sep = "_")
    print(my.outprefix)
    curr_name <- names(data.list[i])
    my.initial <- data.list[[i]]
    
    # Remove extra columns (accommodates the format of pseudobulked single cell sets)
    if (colnames(my.initial)[2] == "Chr") {
      my.data <- my.initial[, -c(2:6)]
    } else {
      my.data <- my.initial
    }
    
    # Remove genes with no reads 
    my.matrix <- my.data[rowSums(my.data[, -1]) > 0, ]
    my.matrix <- my.matrix[, -1]
    my.goodlist <- which(rowSums(my.data[, -1]) > 0)
    rownames(my.matrix) <- my.data$Geneid[my.goodlist]
    
    # Extract ages
    months <- sapply(colnames(my.matrix), extract_months)
    
    # Store age information
    sample_cols <- colnames(my.data)[-1]
    ages <- sapply(sample_cols, extract_months)
    temp_df <- data.frame(
      dataset = rep(curr_name, length(ages)),
      age_group = ages,
      stringsAsFactors = FALSE
    )
    dataset_age_df <- rbind(dataset_age_df, temp_df)
    
    # Check number of replicates per age group
    age_counts <- table(ages)
    
    # Print the replicate counts for visibility
    print(age_counts)
    
    # Skip dataset if any age group has < 2 replicates
    if (any(age_counts < 2)) {
      message("Skipping dataset ", curr_name, ": fewer than 2 replicates in at least one age group")
      next
    }
    
    # Get min/max age
    youngest_age <- min(months)
    oldest_age <- max(months)
    young_columns <- which(months == youngest_age)
    old_columns <- which(months == oldest_age)
    
    # Read depth QC
    young.rd <- colSums(my.matrix[young_columns])
    old.rd <- colSums(my.matrix[old_columns])
    avg.young.rd <- sum(young.rd) / length(young_columns)
    avg.old.rd <- sum(old.rd) / length(old_columns)
    percent.dif <- (avg.old.rd - avg.young.rd) / avg.old.rd
    percent_difs[[curr_name]] <- percent.dif
    
    # Skip if read depth discrepancy
    if (percent.dif < -0.5 || percent.dif > 0.5) {
      print(paste(curr_name, "contains discrepancy in raw read depth between old and young"))
      # dont skip, just a warning 
    }
    
    # Remove low read depth samples
    my.best <- colSums(my.matrix) > 1000000 
    my.filtered.matrix <- my.matrix[, my.best]
    
    # Update age array
    months <- sapply(colnames(my.filtered.matrix), extract_months)
    age_array <- as.numeric(months)
    
    # Build design matrix and check replicates
    dataDesign <- data.frame(
      row.names = colnames(my.filtered.matrix),
      age = age_array
    )
    
    # SVA processing
    mod1 <- model.matrix(~ age, data = dataDesign)
    n.sv.be <- num.sv(my.filtered.matrix, mod1, method = "be")
    
    if (n.sv.be > 0) {
      my.svseq <- svaseq(as.matrix(my.filtered.matrix), mod1, n.sv = n.sv.be, constant = 0.1)
      my.clean <- removeBatchEffect(
        log2(my.filtered.matrix + 0.1),
        batch = NULL,
        covariates = cbind(my.svseq$sv),
        design = mod1
      )
      my.filtered.sva <- round(2^my.clean - 0.1)
    } else {
      my.filtered.sva <- my.filtered.matrix
    }
    
    # DESeq2 analysis
    dds <- DESeqDataSetFromMatrix(
      countData = my.filtered.sva,
      colData = dataDesign,
      design = ~ age
    )
    
    
    dds.deseq <- DESeq(dds)
    
    # Plot dispersion
    pdf(paste(my.outprefix, "dispersion_plot.pdf", sep = "_"))
    plotDispEsts(dds.deseq)
    dev.off()
    
    # Get normalized counts
    tissue.cts <- getVarianceStabilizedData(dds.deseq)
    write.csv(as.data.frame(tissue.cts), file = paste0(my.outprefix, "_norm_counts.csv"))
    
    
    # Check normalized read depth for BMDM datasets
    if (grepl("BMDM_NIA", curr_name)) {
      post_months <- sapply(colnames(tissue.cts), extract_months)
      young_columns_post <- which(post_months == min(post_months))
      old_columns_post <- which(post_months == max(post_months))
      
      young.rd.post <- colSums(tissue.cts[, young_columns_post])
      old.rd.post <- colSums(tissue.cts[, old_columns_post, drop = FALSE])
      avg.young.rd.post <- sum(young.rd.post) / length(young_columns_post)
      avg.old.rd.post <- sum(old.rd.post) / length(old_columns_post)
      post.percent.dif <- (avg.old.rd.post - avg.young.rd.post) / avg.old.rd.post
      print(paste("Post-SVA read depth difference:", post.percent.dif))
    }
    
    # Create color scheme for plots
    months_plot <- sapply(colnames(tissue.cts), extract_months)
    unique_ages <- sort(unique(months_plot))
    scaled_ages <- (unique_ages - min(unique_ages)) / (max(unique_ages) - min(unique_ages))
    
    young_color <- rgb(186, 85, 211, maxColorValue = 255)
    old_color <- rgb(85, 107, 47, maxColorValue = 255)
    color_gradient <- colorRampPalette(c(young_color, old_color))
    gradient_colors <- color_gradient(length(unique_ages))
    
    age_to_color <- setNames(gradient_colors, unique_ages)
    my.colors <- age_to_color[as.character(months_plot)]
    legend_ages <- sort(unique(months_plot))
    legend_colors <- age_to_color[as.character(legend_ages)]
    
    # MDS plot
    mds.result <- cmdscale(1 - cor(tissue.cts, method = "spearman"), k = 2)
    pdf(paste(my.outprefix, "MDS_plot.pdf", sep = "_"))
    plot(mds.result[, 1], mds.result[, 2],
         xlab = "MDS dimension 1", ylab = "MDS dimension 2",
         main = "Multi-dimensional Scaling",
         cex = 3, pch = 16, col = my.colors,
         xlim = c(-0.06, 0.06), ylim = c(-0.06, 0.06),
         cex.lab = 1.5, cex.axis = 1.5)
    legend("topright", legend = legend_ages, fill = legend_colors,
           title = "Age (months)", border = NA, cex = 1.2)
    dev.off()
    
    # Normalized counts boxplot
    pdf(paste(my.outprefix, "_Normalized_counts_boxplot.pdf", sep = ""))
    par(mar = c(5, 4, 4, 8))
    boxplot(tissue.cts, col = my.colors, cex = 0.5,
            ylab = "Log2 DESeq2 Normalized counts", las = 2)
    par(xpd = TRUE)
    legend("topright", inset = c(-0.3, 0), legend = legend_ages,
           fill = legend_colors, title = "Age (months)", cex = 0.8, bty = "n")
    dev.off()
    
    # Xist/Ddx3y barplot
    xist_present <- "Xist" %in% rownames(tissue.cts)
    ddx3y_present <- "Ddx3y" %in% rownames(tissue.cts)
    
    pdf(paste(my.outprefix, "_Xist_Ddx3y_Barplot.pdf", sep = ""))
    par(mar = c(5, 4, 4, 8))
    
    if (xist_present && ddx3y_present) {
      xist_counts <- tissue.cts["Xist", ]
      ddx3y_counts <- tissue.cts["Ddx3y", ]
      bar_heights <- rbind(xist_counts, ddx3y_counts)
      bar_colors <- rbind(my.colors, my.colors)
      bar_densities <- rbind(rep(20, length(xist_counts)), rep(NA, length(ddx3y_counts)))
      
      barplot(bar_heights, beside = TRUE, col = bar_colors, density = bar_densities,
              main = paste(curr_name, "Xist vs Ddx3y Expression"),
              ylab = "Normalized log2(counts)", las = 2)
    } else if (xist_present) {
      xist_counts <- tissue.cts["Xist", ]
      barplot(xist_counts, col = my.colors, density = 20,
              main = paste(curr_name, "Xist Expression"),
              ylab = "Normalized log2(counts)", las = 2)
    } else if (ddx3y_present) {
      ddx3y_counts <- tissue.cts["Ddx3y", ]
      barplot(ddx3y_counts, col = my.colors,
              main = paste(curr_name, "Ddx3y Expression"),
              ylab = "Normalized log2(counts)", las = 2)
    }
    
    par(xpd = TRUE)
    legend("topright", inset = c(-0.3, 0), legend = legend_ages,
           fill = legend_colors, title = "Age (months)", cex = 0.8, bty = "n")
    dev.off()
    
    # DEG analysis
    res.age <- results(dds.deseq, name = "age")
    res.age <- res.age[!is.na(res.age$padj), ]
    genes.age <- rownames(res.age)[res.age$padj < 0.05]
    deg_counts[[curr_name]] <- length(genes.age)
    
    # Save results
    save(res.age, file = paste(my.outprefix, "_AGE.RData", sep = "_"))
    write.table(tissue.cts, file = paste(my.outprefix, "_log2_counts_matrix_DEseq2_SVA.txt", sep = "_"),
                sep = "\t", row.names = TRUE, quote = FALSE)
    write.table(res.age, file = paste(my.outprefix, "_AGE_DIM_all_genes_statistics.txt", sep = "_"),
                sep = "\t", row.names = TRUE, quote = FALSE)
       write.csv(res.age[order(res.age$padj), ], 
              file = paste0(my.outprefix, "_AGE_DIM_all_genes_statistics_sorted.csv"), 
              row.names = TRUE)
    write.table(res.age[genes.age, ], file = paste(my.outprefix, "_AGE_DIM_FDR5_genes_statistics.txt", sep = "_"),
                sep = "\t", row.names = TRUE, quote = FALSE)
    
    
  }
  
  # Return results
  return(list(
    age_info = dataset_age_df,
    percent_differences = percent_difs,
    deg_counts = deg_counts
  ))
}


################################################################################
# 4. Run the function for male and female datasets separately
################################################################################

# Initialize empty lists to store results for each dataset
male_results_list <- list()
female_results_list <- list()


# Process each male dataset one by one
for(dataset_name in names(my.male.data.list)) {
  cat("Processing male dataset:", dataset_name, "\n")
  
  # Create a single-item list for this dataset
  single_dataset_list <- list(my.male.data.list[[dataset_name]])
  names(single_dataset_list) <- dataset_name
  
  # Process one dataset at a time
  result <- process_rnaseq_dataset(single_dataset_list, sex = "male")
  
  # Store the result
  male_results_list[[dataset_name]] <- result
}


# Process each female dataset one by one  
for(dataset_name in names(my.female.data.list)) {
  cat("Processing female dataset:", dataset_name, "\n")
  
  # Create a single-item list for this dataset
  single_dataset_list <- list(my.female.data.list[[dataset_name]])
  names(single_dataset_list) <- dataset_name
  
  # Process just this one dataset
  result <- process_rnaseq_dataset(single_dataset_list, sex = "female")
  
  # Store the result
  female_results_list[[dataset_name]] <- result
}


################################################################################
# 5. Create QC summary tables
################################################################################

generate_summary_tables <- function(results_list, sex_label) {
  library(dplyr)
  
  # Combine all age_info into one data frame
  age_info_table <- do.call(rbind, lapply(results_list, function(x) x$age_info))
  
  # Collapse age info
  collapsed_age_info <- age_info_table %>%
    group_by(dataset) %>%
    summarise(age_groups = paste(age_group, collapse = ",")) %>%
    ungroup()
  
  # DEG counts table
  deg_counts_table <- data.frame(
    dataset = names(results_list),
    deg_count = sapply(results_list, function(x) {
      deg_entry <- x$deg_counts
      if (length(deg_entry) == 0) return(NA)
      return(as.integer(deg_entry[[1]]))
    })
  ) %>%
    arrange(deg_count)  # Sort here
  
  # Percent differences table
  percent_diff_table <- data.frame(
    dataset = names(results_list),
    percent_difference = sapply(results_list, function(x) {
      perc_entry <- x$percent_differences
      if (length(perc_entry) == 0) return(NA)
      return(as.numeric(perc_entry[[1]]))
    })
  )
  
  # Save outputs to CSV
  write.csv(collapsed_age_info, paste0(Sys.Date(), "_collapsed_age_info_", sex_label, ".csv"), row.names = FALSE)
  write.csv(deg_counts_table, paste0(Sys.Date(), "_deg_counts_table_", sex_label, ".csv"), row.names = FALSE)
  write.csv(percent_diff_table, paste0(Sys.Date(), "_percent_diff_table_", sex_label, ".csv"), row.names = FALSE)
  
  # Return as list 
  return(list(
    age_info = collapsed_age_info,
    deg_counts = deg_counts_table,
    percent_differences = percent_diff_table
  ))
}

male_tables <- generate_summary_tables(male_results_list, "male")
female_tables <- generate_summary_tables(female_results_list, "female")

# Filter male datasets: non-NA and deg_count >= 100
male_deg_pass <- male_tables$deg_counts %>%
  filter(!is.na(deg_count) & deg_count >= 100) %>%
  mutate(sex = "Male")

# Filter female datasets: non-NA and deg_count >= 100
female_deg_pass <- female_tables$deg_counts %>%
  filter(!is.na(deg_count) & deg_count >= 100) %>%
  mutate(sex = "Female")

# Combine the two
combined_passed_datasets <- bind_rows(male_deg_pass, female_deg_pass) %>%
  dplyr::select(dataset, sex)


# Save to text file
writeLines(rownames(combined_passed_datasets), paste0(Sys.Date(), "_datasets_pass_qc.txt"))


#################################################################################################
# Save package versions
sink(file = paste(Sys.Date(),"Preprocessing_Session_Info.txt", sep =""))
sessionInfo()
sink()

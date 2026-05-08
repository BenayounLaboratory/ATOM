####################################################################################
# Aging transcriptomics of mouse macrophages 
# Pre-processing - check replicates, Xist/Ddx3y, normalize, read depth check, DEGs
####################################################################################

set.seed(1234) # set seed for reproducibility

library('sva')
library('limma')
library(DESeq2)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/final_filtered_counts")

################################################################################
# 1. Read in the data
################################################################################

#Males
GSE128830_Peritoneal    <- read.csv("2025-03-07_GSE128830_Peritoneal_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE93202_Spleen         <- read.csv("2021-06-04_GSE93202_Spleen_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE93202_VAT            <- read.csv("2021-06-04_GSE93202_VAT_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE98249_BAM            <- read.csv("2021-06-04_GSE98249_BAM_macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE98249_BMM            <- read.csv("2021-06-04_GSE98249_BMM_macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE98401_Brain          <- read.csv("2021-06-04_GSE98401_Brain_microglia_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE134397_Alveolar_CTL  <- read.csv("2021-06-04_GSE134397_Alveolar_Macrophages_CTL_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE137028_Microglia     <- read.csv("2021-06-04_GSE137028_Microglia_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE154832_eWAT          <- read.csv("2021-06-04_GSE154832_eWAT_Phagocytic_SVF_aging_counts.txt", sep = "\t", header = T, skip = 1)
PRJNA682234_Callus      <- read.csv("2021-06-04_PRJNA682234_Callus_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
PRJNA800823_SkM         <- read.csv("2025-08-07_PRJNA800823_Skeletal_Muscle_macrophages_aging_counts.txt", sep = "\t", header = T)
PRJNA1173774_BMDM       <- read.csv("2025-12-11_PRJNA1173774_BMDM_aging_counts_mm10.txt", sep = "\t", header = T)

#Females
GSE199763_SkinWound         <- read.csv("2025-03-07_GSE199763_SkinWound_Macrophages_aging_Female_counts.txt", sep = "\t", header = T, skip = 1)
GSE199879_Spleen_Red_Pulp   <- read.csv("2025-03-07_GSE199879_Spleen_Red_Pulp_Macrophages_aging_Female_counts.txt", sep = "\t", header = T, skip = 1)
GSE132882_Nerve             <- read.csv("2021-06-04_GSE132882_Nerve_Macrophage_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE156762_Microglia         <- read.csv("2021-08-06_GSE156762_Microglia_JCI_counts.txt", sep = "\t", header = T, skip = 1)

#MixedSex
GSE124829_ImmGen_Peritoneal <- read.csv("2021-06-04_GSE124829_ImmGen_Peritoneal_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE131869_Microglia         <- read.csv("2021-06-04_GSE131869_Microglia_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE267529_Microglia         <- read.csv("2021-07-27_GSE267529_Microglia_Goodridge_counts.txt", sep = "\t", header = T, skip = 1)
BMDM_NIA                    <- read.csv("2024-03-28_BMDM_NIA_aging_CLEAN_counts.txt", sep = "\t", header = T, skip = 1)

#Unknown Sex
GSE190689_Alveolar  <- read.csv("2025-03-07_GSE190689_Alveolar_Macrophages_aging_UnkSex_counts.txt", sep = "\t", header = T, skip = 1)
GSE124872_Alveolar  <- read.csv("2021-06-04_GSE124872_Alveolar_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE142580_SkM       <- read.csv("2021-06-04_GSE142580_SkM_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE145295_Alveolar  <- read.csv("2021-06-04_GSE145295_Alveolar_Macrophages_aging_counts.txt", sep = "\t", header = T, skip = 1)
GSE134397_Alveolar  <- read.csv("2025-03-07_GSE134397_Alveolar_Macrophages_aging_Mixed_Sex_counts.txt", sep = "\t", header = T, skip = 1)
PRJNA847895_BMDM    <- read.csv("2025-12-11_PRJNA847895_BMDM_aging_counts_mm10.txt", sep = "\t", header = T)
PRJNA524906_microglia <- read.csv("2025-12-11_PRJNA524906_microglia_aging_counts_mm10.txt", sep = "\t", header = T)
PRJNA816431_callus  <- read.csv("2025-12-11_PRJNA816431_callus_aging_counts_mm10.txt", sep = "\t", header = T)
GSE205395_SkM_NEW <- read.csv("2026-04-16_GSE205395_SkM_macrophages_injured_counts.txt", sep = "\t", header = T)

sex_key <- data.frame(
  dataset = c(
    "GSE128830_Peritoneal",
    "GSE93202_Spleen",
    "GSE93202_VAT",
    "GSE98249_BAM",
    "GSE98249_BMM",
    "GSE98401_Microglia",
    "GSE134397_Alveolar_CTL",
    "GSE137028_Microglia",
    "GSE154832_eWAT",
    "PRJNA682234_Callus",
    "PRJNA800823_SkM",
    "PRJNA1173774_BMDM",
    "GSE199763_SkinWound",
    "GSE199879_Spleen_Red_Pulp",
    "GSE132882_Nerve",
    "GSE156762_Microglia",
    
    "GSE124829_ImmGen_Peritoneal_M",
    "GSE124829_ImmGen_Peritoneal_F",
    "GSE131869_Microglia_M",
    "GSE131869_Microglia_F",
    "GSE267529_Microglia_M",
    "GSE267529_Microglia_F",
    "BMDM_NIA_M",
    "BMDM_NIA_F",
    "GSE190689_Alveolar",
    "GSE124872_Alveolar",
    "GSE142580_SkM",
    "GSE145295_Alveolar",
    "GSE134397_Alveolar",
    "PRJNA847895_BMDM",
    "PRJNA524906_Microglia",
    "PRJNA816431_Callus",
    "GSE205395_SkM_NEW"
  ),
  reported_sex = c(
    "Male",
    "Male",
    "Male",
    "Male",
    "Male",
    "Male",
    "Male",
    "Male",
    "Male",
    "Male",
    "Male",
    "Male",
    "Female",
    "Female",
    "Female",
    "Female",
    # Split datasets
    "Mixed",
    "Mixed",
    "Mixed",
    "Mixed",
    "Mixed",
    "Mixed",
    "Mixed",
    "Mixed",
    "Unknown",
    "Unknown",
    "Unknown",
    "Unknown",
    "Unknown",
    "Unknown",
    "Unknown",
    "Unknown",
    "Female"
  ),
  stringsAsFactors = FALSE
)


# Construct the path with today's date
outdir <- file.path("/Users/ellaschwab/Benayoun_Lee_Local/ATOM", paste0(Sys.Date(), "_continuous_age_output_REFACTORED"))

# Create the output directory 
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

setwd(outdir)

################################################################################
# 2. Separate males and females from the mixed sex datasets
################################################################################

my.both.data.list <- list("GSE124829_ImmGen_Peritoneal" = GSE124829_ImmGen_Peritoneal,
                          "GSE131869_Microglia"         = GSE131869_Microglia,
                          "GSE267529_Microglia"         = GSE267529_Microglia,
                          "BMDM_NIA"                    = BMDM_NIA)

for (i in 1:length(my.both.data.list)) {
  dataset_name <- names(my.both.data.list[i])
  my.data      <- my.both.data.list[[i]]
  
  femalecols <- grep("_female", colnames(my.data), value = FALSE, ignore.case = TRUE)
  if (length(femalecols) == 0) {
    femalecols  <- grep("_YF", colnames(my.data), value = FALSE)
    femalecols  <- append(femalecols, grep("_OF", colnames(my.data), value = FALSE))
  }
  
  female_data <- my.data[, c(1, 2, 3, 4, 5, 6, femalecols)]
  male_data   <- my.data[, -femalecols]
  
  assign(paste0(dataset_name, "_F"), female_data)
  assign(paste0(dataset_name, "_M"), male_data)
}

################################################################################
# 3. Build sex-stratified data lists
################################################################################

my.male.data.list <- list(
  "GSE93202_Spleen"                  = GSE93202_Spleen,
  "GSE93202_VAT"                     = GSE93202_VAT,
  "GSE98401_Microglia"               = GSE98401_Brain,
  "GSE134397_Alveolar_CTL"           = GSE134397_Alveolar_CTL,
  "GSE124829_ImmGen_Peritoneal_M"    = GSE124829_ImmGen_Peritoneal_M,
  "GSE131869_Microglia_M"            = GSE131869_Microglia_M,
  "GSE128830_Peritoneal"             = GSE128830_Peritoneal,
  "GSE137028_Microglia"              = GSE137028_Microglia,
  "GSE154832_eWAT"                   = GSE154832_eWAT,
  "PRJNA682234_Callus"               = PRJNA682234_Callus,
  "GSE267529_Microglia_M"            = GSE267529_Microglia_M,
  "BMDM_NIA_M"                       = BMDM_NIA_M,
  "GSE145295_Alveolar"               = GSE145295_Alveolar,
  "GSE142580_SkM"                    = GSE142580_SkM,
  "GSE124872_Alveolar"               = GSE124872_Alveolar,
  "GSE190689_Alveolar"               = GSE190689_Alveolar,
  "PRJNA524906_Microglia"            = PRJNA524906_microglia,
  "PRJNA816431_Callus"               = PRJNA816431_callus,
  "PRJNA800823_SkM"                  = PRJNA800823_SkM,
  "PRJNA1173774_BMDM"                = PRJNA1173774_BMDM,
  "GSE98249_BAM"                     = GSE98249_BAM,   # problematic
  "GSE98249_BMM"                     = GSE98249_BMM    # problematic
)

my.female.data.list <- list(
  "GSE199763_SkinWound"              = GSE199763_SkinWound,
  "GSE199879_Spleen_Red_Pulp"        = GSE199879_Spleen_Red_Pulp,
  "GSE156762_Microglia"              = GSE156762_Microglia,
  "GSE267529_Microglia_F"            = GSE267529_Microglia_F,
  "GSE132882_Nerve"                  = GSE132882_Nerve,
  "GSE124829_ImmGen_Peritoneal_F"    = GSE124829_ImmGen_Peritoneal_F,
  "GSE131869_Microglia_F"            = GSE131869_Microglia_F,
  "BMDM_NIA_F"                       = BMDM_NIA_F,
  "GSE134397_Alveolar"               = GSE134397_Alveolar,
  "PRJNA847895_BMDM"                 = PRJNA847895_BMDM,
  "GSE205395_SkM_NEW"                = GSE205395_SkM_NEW
)

################################################################################
# 4. Helper functions
################################################################################

extract_months <- function(colname) {
  if (grepl("_(\\d+\\.?\\d*)m_", colname)) {
    match <- regmatches(colname, regexpr("_(\\d+\\.?\\d*)m_", colname))
    as.numeric(gsub("[^0-9\\.]", "", match))
  } else if (grepl("_Old_", colname)) {
    25
  } else if (grepl("_Young_", colname)) {
    4
  } else if (grepl("_OM", colname)) {
    20
  } else if (grepl("_YM", colname)) {
    4
  } else if (grepl("_OF", colname)) {
    20
  } else if (grepl("_YF", colname)) {
    4
  } else if (grepl("_(\\d+)_months_", colname)) {
    match <- regmatches(colname, regexpr("_(\\d+)_months_", colname))
    as.numeric(gsub("[^0-9\\.]", "", match))
  } else if (grepl("_(\\d+)months_", colname)) {
    match <- regmatches(colname, regexpr("_(\\d+)months_", colname))
    as.numeric(gsub("[^0-9\\.]", "", match))
  } else if (grepl("_(\\d+)m(\\d+)_", colname)) {
    match <- regmatches(colname, regexec("_(\\d+)m(\\d+)_", colname))[[1]]
    as.numeric(match[2])
  } else {
    NA
  }
}

get_gene <- function(mat, gene) {
  if (gene %in% rownames(mat)) {
    return(mat[gene, ])
  } else if (toupper(gene) %in% toupper(rownames(mat))) {
    return(mat[toupper(rownames(mat)) == toupper(gene), ][1, ])
  } else {
    return(rep(NA, ncol(mat)))
  }
}

################################################################################
# 5. Main processing function — runs up to and including VST transformation.
#    Returns a named list; each element is a list with:
#      $vst_matrix  - VST-transformed count matrix (genes x samples)
#      $age         - numeric age vector aligned to vst_matrix columns
#      $sex         - sex label passed in
#      $dds         - fitted DESeqDataSet 
################################################################################

process_rnaseq_dataset <- function(data.list, sex, outdir, outprefix.base = Sys.Date()) {
  
  results_list <- list()
  
  for (i in seq_along(data.list)) {
    
    curr_name    <- names(data.list)[i]
    my.outprefix <- paste(outprefix.base, curr_name, sep = "_")
    message("Processing: ", my.outprefix)
    
    my.initial <- data.list[[i]]
    
    # Remove featureCounts annotation columns if present
    my.data <- if (colnames(my.initial)[2] == "Chr") my.initial[, -c(2:6)] else my.initial
    
    # Remove genes with no reads; set rownames to Geneid
    keep_rows  <- rowSums(my.data[, -1]) > 0
    my.matrix  <- my.data[keep_rows, -1]
    rownames(my.matrix) <- my.data$Geneid[keep_rows]
    
    # Extract ages from column names
    months <- sapply(colnames(my.matrix), extract_months)
    
    # Check replicates (need >= 2 per age group)
    age_counts <- table(months)
    print(age_counts)
    if (any(age_counts < 2)) {
      message("Skipping ", curr_name, ": insufficient replicates")
      next
    }
    
    # Filter samples with low library depth
    removed_samples <- colnames(my.matrix)[colSums(my.matrix) <= 1e6]
    if (length(removed_samples) > 0) {
      message("  Removing low-depth samples from ", curr_name, ": ",
              paste(removed_samples, collapse = ", "))
    }
    my.filtered.matrix <- my.matrix[, colSums(my.matrix) > 1e6]
    age_array <- as.numeric(sapply(colnames(my.filtered.matrix), extract_months))
    
    dataDesign <- data.frame(row.names = colnames(my.filtered.matrix), age = age_array)
    
    # ------------------------------------------------------------------
    # SVA: estimate and remove unwanted variation
    # ------------------------------------------------------------------
    mod1     <- model.matrix(~ age, data = dataDesign)
    n.sv.be  <- num.sv(my.filtered.matrix, mod1, method = "be")
    
    if (n.sv.be > 0) {
      my.svseq       <- svaseq(as.matrix(my.filtered.matrix), mod1, n.sv = n.sv.be, constant = 0.1)
      my.clean       <- removeBatchEffect(log2(my.filtered.matrix + 0.1),
                                          batch      = NULL,
                                          covariates = cbind(my.svseq$sv),
                                          design     = mod1)
      my.filtered.sva <- round(2^my.clean - 0.1)
    } else {
      my.filtered.sva <- my.filtered.matrix
    }
    
    # ------------------------------------------------------------------
    # DESeq2
    # ------------------------------------------------------------------
    dds <- DESeqDataSetFromMatrix(countData = my.filtered.sva,
                                  colData   = dataDesign,
                                  design    = ~ age)
    dds.deseq <- DESeq(dds)
    
    # Dispersion plot
    pdf(file.path(outdir, paste0(my.outprefix, "_dispersion_plot.pdf")))
    plotDispEsts(dds.deseq)
    dev.off()
    
    # ------------------------------------------------------------------
    # VST transformation
    # ------------------------------------------------------------------
    vst_mat <- getVarianceStabilizedData(dds.deseq)
    
    # Save VST matrix as a plain CSV for easy re-loading
    vst_outfile <- file.path(outdir, paste0(my.outprefix, "_VST_matrix.csv"))
    write.csv(vst_mat, vst_outfile)
    message("  VST matrix saved to: ", vst_outfile)
    
    # ------------------------------------------------------------------
    # Collect results 
    # ------------------------------------------------------------------
    results_list[[curr_name]] <- list(
      vst_matrix = vst_mat,
      age        = age_array,
      sex        = sex,
      dds        = dds.deseq
    )
  }
  
  return(results_list)
}

###################################################################################
# 6. DEG calculation and read depth check function — second pass over the results list produced
#    by process_rnaseq_dataset.  Adds $res_df to each entry in-place and
#    writes per-dataset CSV files.
#    Returns the same results_list with $res_df appended to every entry.
###################################################################################

calculate_degs_and_depth <- function(results_list, outdir, outprefix.base = Sys.Date()) {
  
  rd_summary_list <- list()
  
  for (curr_name in names(results_list)) {
    
    entry        <- results_list[[curr_name]]
    my.outprefix     <- paste(outprefix.base, curr_name, sep = "_")   # used for read depth
    my.deg.outprefix <- paste(Sys.Date(), curr_name, sep = "_")        # used for DEG files
    message("Calculating DEGs: ", my.outprefix)
    
    # ------------------------------------------------------------------
    # Re-derive age from VST matrix column names to stay in sync with
    # any samples removed during sex-mislabeling filtering 
    # ------------------------------------------------------------------
    vst_mat  <- entry$vst_matrix
    age_arr  <- as.numeric(sapply(colnames(vst_mat), extract_months))
    
    # ------------------------------------------------------------------
    # Normalized read depth QC (on VST matrix)
    # ------------------------------------------------------------------
    young_age      <- min(age_arr)
    old_age        <- max(age_arr)
    young_cols_vst <- which(age_arr == young_age)
    old_cols_vst   <- which(age_arr == old_age)
    
    avg.young.rd.vst <- mean(colSums(vst_mat[, young_cols_vst, drop = FALSE]))
    avg.old.rd.vst   <- mean(colSums(vst_mat[, old_cols_vst,   drop = FALSE]))
    percent.dif.vst  <- (avg.old.rd.vst - avg.young.rd.vst) / avg.old.rd.vst
    
    message("  Normalized read depth (VST colSums) — young (", young_age, "m): ",
            round(avg.young.rd.vst, 2), " | old (", old_age, "m): ",
            round(avg.old.rd.vst, 2))
    message("  Young vs old % difference: ", round(percent.dif.vst * 100, 2), "%")
    
    if (abs(percent.dif.vst) > 0.5) {
      warning("  WARNING: >50% normalized read depth discrepancy between young and old in ", curr_name)
    }
    
    rd_summary_list[[curr_name]] <- data.frame(
      dataset               = curr_name,
      sex                   = entry$sex,
      young_age_months      = young_age,
      old_age_months        = old_age,
      avg_young_rd          = round(avg.young.rd.vst, 2),
      avg_old_rd            = round(avg.old.rd.vst, 2),
      pct_diff_young_vs_old = round(percent.dif.vst * 100, 2),
      stringsAsFactors      = FALSE
    )
    
    # ------------------------------------------------------------------
    # DEG calculation
    # ------------------------------------------------------------------
    res_df <- as.data.frame(results(entry$dds))
    res_df <- res_df[!is.na(res_df$padj), ]
    sig    <- res_df[res_df$padj < 0.05, ]
    
    write.table(res_df,
                file.path(outdir, paste0(my.deg.outprefix, "_AGE_DIM_all_genes_statistics.txt")),
                sep = "\t", row.names = TRUE, quote = FALSE)
    
    write.table(sig,
                file.path(outdir, paste0(my.deg.outprefix, "_AGE_DIM_FDR5_genes_statistics.txt")),
                sep = "\t", row.names = TRUE, quote = FALSE)
    
    write.csv(res_df[order(res_df$padj), ],
              file.path(outdir, paste0(my.deg.outprefix, "_AGE_DIM_all_genes_statistics_sorted.csv")),
              row.names = TRUE)
    
    results_list[[curr_name]]$res_df <- res_df
  }
  
  rd_summary_df <- do.call(rbind, rd_summary_list)
  write.csv(rd_summary_df,
            file.path(outdir, paste0(outprefix.base, "_read_depth_summary.csv")),
            row.names = FALSE)
  message("Read depth summary saved.")
  
  return(results_list)
}

################################################################################
# 7. Xist / Ddx3y heatmap function  
#    Accepts the list returned by process_rnaseq_dataset, or a single element.
################################################################################

plot_sex_marker_heatmaps <- function(results_list, outdir, sex_key, outprefix.base = Sys.Date(),
                                     cellwidth = 3, cellheight = 8) {
  
  library(ComplexHeatmap)
  library(circlize)
  
  base_names <- unique(sub("_(M|F)$", "", names(results_list)))
  plot_list  <- list()
  
  for (bn in base_names) {
    member_names <- names(results_list)[sub("_(M|F)$", "", names(results_list)) == bn]
    member_names <- member_names[order(match(
      sapply(member_names, function(n) results_list[[n]]$sex),
      c("female", "male", "unknown")
    ))]
    
    mat_list <- list()
    author_sex_list <- list()
    age_list <- list()
    sex_list <- list()
    
    for (curr_name in member_names) {
      entry   <- results_list[[curr_name]]
      vst_mat <- entry$vst_matrix
      age_arr <- as.numeric(sapply(colnames(vst_mat), extract_months))
      sex     <- entry$sex
      
      xist  <- get_gene(vst_mat, "Xist")
      ddx3y <- get_gene(vst_mat, "Ddx3y")
      
      mat     <- rbind(Xist = xist, Ddx3y = ddx3y)
      ord     <- order(age_arr)
      mat     <- mat[, ord, drop = FALSE]
      age_arr <- age_arr[ord]
      
      mat_list[[curr_name]] <- mat
      age_list[[curr_name]] <- age_arr
      sex_list[[curr_name]] <- rep(sex, length(age_arr))
      # inside the for (curr_name in member_names) loop, after sex_list append:
      author_sex_list[[curr_name]] <- rep(
        sex_key$reported_sex[sex_key$dataset == curr_name],
        length(age_arr)
      )
    }
    
    heat_mat    <- do.call(cbind, mat_list)
    age_ordered <- unlist(age_list)
    sex_labels  <- unlist(sex_list)
    author_sex_labels <- unlist(author_sex_list)
    
    xist_mat  <- heat_mat["Xist",  , drop = FALSE]
    ddx3y_mat <- heat_mat["Ddx3y", , drop = FALSE]
    
    unique_ages  <- sort(unique(age_ordered[!is.na(age_ordered)]))
    age_to_color <- setNames(
      colorRampPalette(c(rgb(186, 85,  211, maxColorValue = 255),
                         rgb( 85, 107,  47, maxColorValue = 255)))(length(unique_ages)),
      as.character(unique_ages)
    )
    
    
    ha <- HeatmapAnnotation(
      Age = factor(as.character(age_ordered), levels = as.character(unique_ages)),
      Sex = factor(sex_labels, levels = c("female", "male", "unknown")),
      AuthorSex = factor(author_sex_labels, levels = c("Female", "Male", "Mixed", "Unknown")),
      col = list(
        Age = age_to_color,
        Sex = c("male" = "deepskyblue", "female" = "deeppink","Unknown" = "grey80"),
        AuthorSex = c("Male" = "deepskyblue", "Female" = "deeppink", 
                      "Mixed" = "purple", "Unknown" = "grey80")
      )
    )
    
    # Separate color scales per gene anchored to their own VST range
    xist_min  <- if (all(is.na(xist_mat)))  0 else min(xist_mat,  na.rm = TRUE)
    xist_max  <- if (all(is.na(xist_mat)))  0 else max(xist_mat,  na.rm = TRUE)
    ddx3y_min <- if (all(is.na(ddx3y_mat))) 0 else min(ddx3y_mat, na.rm = TRUE)
    ddx3y_max <- if (all(is.na(ddx3y_mat))) 0 else max(ddx3y_mat, na.rm = TRUE)
    
    if (xist_min  == xist_max)  { xist_min  <- xist_min  - 0.01; xist_max  <- xist_max  + 0.01 }
    if (ddx3y_min == ddx3y_max) { ddx3y_min <- ddx3y_min - 0.01; ddx3y_max <- ddx3y_max + 0.01 }
    
    xist_breaks  <- seq(floor(xist_min),  ceiling(xist_max),  by = 1)
    xist_cols    <- colorRampPalette(c("#FFF0F0", "#FFBBBB"))(length(xist_breaks))
    xist_col_fun <- colorRamp2(xist_breaks, xist_cols)
    
    ddx3y_breaks  <- seq(floor(ddx3y_min), ceiling(ddx3y_max), by = 1)
    ddx3y_cols    <- colorRampPalette(c("#F0F5FF", "#BBCCEE"))(length(ddx3y_breaks))
    ddx3y_col_fun <- colorRamp2(ddx3y_breaks, ddx3y_cols)
    
    ht_xist <- Heatmap(
      xist_mat,
      col               = xist_col_fun,
      name              = "Xist VST",
      top_annotation    = ha,
      cluster_rows      = FALSE,
      cluster_columns   = FALSE,
      show_column_names = FALSE,
      row_names_gp      = gpar(fontsize = 10),
      border            = TRUE,
      width             = unit(cellwidth  * ncol(xist_mat),  "mm"),
      height            = unit(cellheight * nrow(xist_mat),  "mm"),
      column_title      = bn,
      na_col            = "grey80"
    )
    
    ht_ddx3y <- Heatmap(
      ddx3y_mat,
      col               = ddx3y_col_fun,
      name              = "Ddx3y VST",
      cluster_rows      = FALSE,
      cluster_columns   = FALSE,
      show_column_names = TRUE,
      column_names_gp   = gpar(fontsize = 5),
      column_names_rot  = 45,
      row_names_gp      = gpar(fontsize = 10),
      border            = TRUE,
      width             = unit(cellwidth  * ncol(ddx3y_mat), "mm"),
      height            = unit(cellheight * nrow(ddx3y_mat), "mm"),
      column_title      = "",
      na_col            = "grey80"
    )
    
    plot_list[[bn]] <- grid.grabExpr(draw(ht_xist %v% ht_ddx3y))
    
    message("  Plot collected for: ", bn)
  }
  
  pdf_path <- file.path(outdir, paste0(outprefix.base, "_Xist_Ddx3y_VST_heatmaps_all.pdf"))
  pdf(pdf_path, width = 14, height = 8)
  
  for (k in seq(1, length(plot_list), by = 2)) {
    grid.newpage()
    if (k + 1 <= length(plot_list)) {
      grid.draw(arrangeGrob(plot_list[[k]], plot_list[[k + 1]], ncol = 2))
    } else {
      grid.draw(arrangeGrob(plot_list[[k]], ncol = 1))
    }
  }
  
  dev.off()
  message("All heatmaps saved to: ", pdf_path)
}


################################################################################
# 8. Run pipeline
################################################################################

# --- Pass 1: normalize, SVA, DESeq2, VST and check sex ---
male_results   <- process_rnaseq_dataset(my.male.data.list,   "male",   outdir)
female_results <- process_rnaseq_dataset(my.female.data.list, "female", outdir)

#######################################################################################
# 9. Handle sex mislabeling — remove bad samples from raw data, then refit the model 
#######################################################################################

#plot_sex_marker_heatmaps(male_results, outdir, sex_key = sex_key)

all_results <- c(male_results, female_results)
plot_sex_marker_heatmaps(all_results, outdir, sex_key = sex_key)


# GSE98249_BAM and GSE98249_BMM have young males, old females — exclude entirely
my.male.data.list.clean <- my.male.data.list[!names(my.male.data.list) %in% c("GSE98249_BAM", "GSE98249_BMM")]

# GSE131869_Microglia_F — one sample with unclear sex
my.female.data.list.clean <- my.female.data.list
my.female.data.list.clean$GSE131869_Microglia_F <- 
  my.female.data.list$GSE131869_Microglia_F[,
                                            colnames(my.female.data.list$GSE131869_Microglia_F) != 
                                              "GSM3823691_Microglia_2m_female_STAR_Aligned.sortedByCoord.out.bam"]

# GSE124872_Alveolar — one sample with unclear sex
my.male.data.list.clean$GSE124872_Alveolar <- 
  my.male.data.list$GSE124872_Alveolar[,
                                       colnames(my.male.data.list$GSE124872_Alveolar) != 
                                         "GSM3557713_Alveolar_Macrophages_3m_NA_STAR_Aligned.sortedByCoord.out.bam"]

# --- Pass 2: Filter mislabeled samples, calculate new VST transformed counts and use for downstream processing
set.seed(1234)
male_results_filtered   <- process_rnaseq_dataset(my.male.data.list.clean, "male",   outdir,
                                                  outprefix.base = paste0(Sys.Date(), "_clean"))
set.seed(1234)
female_results_filtered <- process_rnaseq_dataset(my.female.data.list.clean, "female", outdir,
                                                  outprefix.base = paste0(Sys.Date(), "_clean"))

################################################################################
# 10. DEG calculation and read depth check
################################################################################

male_results_filtered   <- calculate_degs_and_depth(male_results_filtered,   outdir, outprefix.base = paste0(Sys.Date(), "_male"))
female_results_filtered <- calculate_degs_and_depth(female_results_filtered, outdir, outprefix.base = paste0(Sys.Date(), "_female"))



################################################################################
# 11. Compile DEG counts across all datasets
################################################################################

compile_deg_summary <- function(results_list) {
  do.call(rbind, lapply(names(results_list), function(curr_name) {
    entry  <- results_list[[curr_name]]
    res_df <- entry$res_df
    n_degs <- sum(!is.na(res_df$padj) & res_df$padj < 0.05)
    data.frame(dataset = curr_name,
               sex     = entry$sex,
               n_DEGs  = n_degs,
               stringsAsFactors = FALSE)
  }))
}

deg_summary <- rbind(compile_deg_summary(male_results_filtered),
                     compile_deg_summary(female_results_filtered))

print(deg_summary)
write.csv(deg_summary, file.path(outdir, paste0(Sys.Date(), "_DEG_summary.csv")), row.names = FALSE)

# Split by sex
male_deg_summary   <- deg_summary[deg_summary$sex == "male",   ]
female_deg_summary <- deg_summary[deg_summary$sex == "female", ]

# Sort each (smallest to largest)
male_deg_summary   <- male_deg_summary[order(male_deg_summary$n_DEGs),     ]
female_deg_summary <- female_deg_summary[order(female_deg_summary$n_DEGs), ]

cat("\n--- MALE DATASETS (sorted by DEGs) ---\n")
print(male_deg_summary)

cat("\n--- FEMALE DATASETS (sorted by DEGs) ---\n")
print(female_deg_summary)

datasets_pass_qc <- deg_summary$dataset[deg_summary$n_DEGs >= 100]

writeLines(
  datasets_pass_qc,
  file.path(outdir, paste0(Sys.Date(), "_datasets_pass_qc.txt"))
)

#################################################################################################
# Save package versions
sink(file = paste(Sys.Date(),"Preprocessing_Session_Info.txt", sep =""))
sessionInfo()
sink()




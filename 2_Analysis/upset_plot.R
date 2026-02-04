####################################################################################
# Aging transcriptomics of mouse macrophages 
# Upset plots showing common genes in each niche group by sex
####################################################################################

library(stringr)
library(ComplexHeatmap)
library(grid)
library(ComplexUpset) 



################################################################################
# 1. Read in the data
################################################################################

#put into list
setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/2026-01-26_continuous_age_output")
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
pass_qc <- readLines("2026-01-26_datasets_pass_qc.txt")
my_deseq_list <- my_deseq_list[names(my_deseq_list) %in% pass_qc]


################################################################################
# 2. Make Upset Plots
################################################################################


setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/UpsetPlot")
my.out.dir <- "/Users/ellaschwab/Benayoun_Lee_Local/ATOM/UpsetPlot"


padj_cutoff <- 0.05


sex_map_alias <- c(
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
  "GSE195507_SkM" = "Male",          
  "GSE279654_BMDM" = "Male",         
  "GSE127542_Microglia" = "Male",    
  "GSE198666_Callus" = "Male",       
  
  # Females
  "GSE199763_SkinWound" = "Female",
  "GSE199879_Spleen_Red_Pulp" = "Female",
  "GSE156762_Microglia" = "Female",
  "GSE267529_Microglia_F" = "Female",
  "GSE131869_Microglia_F" = "Female",
  "PRJNA1029936" = "Female",         
  "GSE134397_Alveolar" = "Female"
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
  "PRJNA682234_Callus"        = "PRJNA682234_Callus", #no GEO accession
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

niche_map <- list(
  
  Microglia = c(
    "GSE98401_Microglia",
    "GSE131869_Microglia_M",
    "GSE131869_Microglia_F",
    "GSE137028_Microglia",
    "GSE156762_Microglia",
    "GSE267529_Microglia_M",
    "GSE267529_Microglia_F",
    "GSE127542_Microglia"     
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
  
  Peritoneal = c(
    "GSE128830_Peritoneal"
  ),
  
  Callus = c(
    "PRJNA682234_Callus",     
    "GSE198666_Callus"        
  ),
  
  Skin = c(
    "GSE199763_SkinWound"
  ),
  
  Skeletal = c(
    "GSE142580_SkM",
    "GSE195507_SkM"           
  ),
  
  BMDM = c(
    "PRJNA1029936",           
    "GSE279654_BMDM"          
  )
)


# Filter each dataset for padj < 0.05
filtered_deseq_list <- lapply(my_deseq_list, function(df) {
  df[df$padj < padj_cutoff & !is.na(df$padj), ]
})

# Rename filtered_deseq_list to match aliases
names(filtered_deseq_list) <- sapply(names(filtered_deseq_list), function(x) {
  if (x %in% names(dataset.aliases)) {
    dataset.aliases[[x]]   
  } else {
    x                      
  }
})



# Loop over each niche
for (niche_name in names(niche_map)) {
  
  datasets_in_niche <- niche_map[[niche_name]]
  
  # Only continue if there are at least 2 datasets
  if (length(datasets_in_niche) < 2) next
  
  # Make a named list of gene vectors
  gene_sets <- lapply(datasets_in_niche, function(ds) {
    if (ds %in% names(filtered_deseq_list)) {
      filtered_deseq_list[[ds]]$row.name  # replace with your gene ID column
    } else {
      character(0)
    }
  })
  names(gene_sets) <- datasets_in_niche
  
  # Skip niche if all datasets have no genes
  if (all(sapply(gene_sets, length) == 0)) next
  
  # Make combination matrix
  m <- make_comb_mat(gene_sets, mode = "intersect")
  
  # Map sex to datasets in the matrix
  datasets_in_matrix <- set_name(m)
  dataset_sex <- sex_map_alias[datasets_in_matrix]
  sex_colors <- c("Male" = "deepskyblue", "Female" = "deeppink")
  
  # Right annotation with only sex
  row_ann <- rowAnnotation(
    Sex = dataset_sex,
    col = list(Sex = sex_colors),
    annotation_name_side = "top"
  )
  
  # Plot UpSet exactly like your example
  p <- UpSet(
    m,
    right_annotation = row_ann
  )
  
  # Save as PDF
  pdf(file.path(my.out.dir, paste0(Sys.Date(), "_UpSet_", niche_name, ".pdf")),
      width = 8, height = 5)
  draw(p)
  dev.off()
}

#######################

sink(file = file.path(my.out.dir, paste0(Sys.Date(), "_R_session_Info_upset_plots.txt")))
sessionInfo()
sink()



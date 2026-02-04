####################################################################################
# Aging transcriptomics of mouse macrophages 
# Transcription Factor Analysis with DecoupleR
####################################################################################

set.seed(1234) # set seed for reproducibility

################################################################################
# 1. Read in the data
################################################################################

library(decoupleR)
library(OmnipathR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(stringr)
library(ComplexHeatmap)
library(fgsea)
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


# rename the datasets to the correct name
names(my_deseq_list) <- dataset.aliases[names(my_deseq_list)]

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




################################################################################
# 2. Set up TF network for mouse
################################################################################

padj_cutoff <- 0.05
n_tfs <- 20
my.out.dir <- "/Users/ellaschwab/Benayoun_Lee_Local/ATOM/TF_output"



### Get TF network data from mouse model
mm.net <- decoupleR::get_collectri(organism='mouse', split_complexes=FALSE) 
mm.net

# # A tibble: 39,961 × 3
# source target   mor
# <chr>  <chr>  <dbl>
#   1 Myc    Tert       1
# 2 Spi1   Bglap2     1
# 3 Spi1   Bglap      1
# 4 Spi1   Bglap3     1
# 5 Smad3  Jun        1
# 6 Smad4  Jun        1
# 7 Stat5a Il2        1
# 8 Stat5b Il2        1
# 9 Rela   Fas        1
# 10 F6Y2T6 Nr0b1      1
# # ℹ 39,951 more rows
# # ℹ Use `print(n = ...)` to see more rows



#################################################
# 3. Infer TFs for each dataset separately
#################################################
# infer pathway activities from the t-values of the DEGs with aging
# for each dataset

Top_TFs_res.mph <- data.frame(matrix(0,0,7))
colnames(Top_TFs_res.mph) <- c("statistic", "source","condition","score","p_value","Reg_Rank","Dataset")


# Set rownames using the 'row.names' column for each element in the list
my_deseq_list <- lapply(my_deseq_list, function(df) {
  rownames(df) <- df$row.names
  return(df)
})

# Loop over DEseq2 results
for (i in 1:length(my_deseq_list)) {
  
  # Run fgsea scoring
  contrast_acts <- run_fgsea(mat     = my_deseq_list[[i]][, 'stat', drop=FALSE],
                             net     = mm.net                                  , 
                             minsize = 5                                             )
  
  # We select the norm_fgsea activities and then we show changes in activity with aging
  # Filter norm_fgsea
  tf.scores         <- contrast_acts[contrast_acts$statistic == 'norm_fgsea',]
  tf.sig            <- tf.scores[tf.scores$p_value < 0.10,]
  tf.sig$Reg_Rank   <- NA
  tf.sig$Dataset  <- names(my_deseq_list)[i]
  
  # Only keep TF regulons who TF is expressed in dataset
  tf.sig <- tf.sig[tf.sig$source %in% row.names(my_deseq_list[[i]]),]
  
  # Add rank information
  # Ranking by decoupleR score for both the upreg and downreg TFs
  msk <- tf.sig$score > 0
  tf.sig$Reg_Rank[msk ]  <- round(rank(-tf.sig[msk, 'score']))
  tf.sig$Reg_Rank[!msk]  <- round(rank(-abs(tf.sig[!msk, 'score'])))
  
  Top_TFs_res.mph <- rbind(Top_TFs_res.mph,tf.sig)
  
  # Filter top significant TFs in both signs
  tfs <- tf.sig %>%
    arrange(Reg_Rank) %>%
    head(n_tfs) %>%
    pull(source)
  
  f_contrast_acts <- tf.sig %>% filter(source %in% tfs)
  
  # Plot top each direction
  tf.reg.plot <- ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + geom_bar(aes(fill = score), stat = "identity") 
  tf.reg.plot <- tf.reg.plot + scale_fill_gradient2(low = "darkslateblue", mid = "whitesmoke", high = "firebrick3", midpoint = 0, limits = c(-5,5)) 
  tf.reg.plot <- tf.reg.plot  + theme(axis.title = element_text(face = "bold", size = 12),
                                      axis.text.x = element_text(angle = 45, hjust = 1, size =10),
                                      axis.text.y = element_text(size =10, face= "bold"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank() ) 
  tf.reg.plot <- tf.reg.plot + xlab("Enriched TF regulons") + ylab("decoupleR score") + theme(plot.title = element_text(hjust = 0.5))
  tf.reg.plot <- tf.reg.plot + ylim(c(-5,5)) + ggtitle(names(my_deseq_list[i])) + coord_flip()
  # tf.reg.plot
  
  pdf(file = file.path(my.out.dir, paste0(Sys.Date(), "_BarPLot_decoupleR_fgsea_Top_", n_tfs, "_TF_regulons_", names(my_deseq_list[i]), "Aging_FDR05.pdf")), 
      width = 4, height = 5)
  print(tf.reg.plot)
  dev.off()
  
}


write.table(Top_TFs_res.mph, 
            file = file.path(my.out.dir, paste0(Sys.Date(), "_decoupleR_fgsea_TFregulons_by_Dataset_MphAging_FDR10.txt")), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)


###############################################################################
# 4. Find recurrent TFs in more than half the datasets at FDR < 5%
###############################################################################

# Utility function to concatenate vector elements into a comma-separated string
concat <- function(x) paste0(x, collapse = ",")

# 1. Filter significant results
Top_TFs_res.mph.5 <- Top_TFs_res.mph %>%
  filter(p_value < 0.05)

# 2. Deduplicate — keep top scoring TF per dataset
Top_TFs_res.mph.5 <- Top_TFs_res.mph.5 %>%
  group_by(source, Dataset) %>%
  slice_max(abs(score), n = 1, with_ties = FALSE) %>%
  ungroup()

# 3. Add sign column (direction of enrichment)
Top_TFs_res.mph.5 <- Top_TFs_res.mph.5 %>%
  mutate(sign = sign(score))

# 4. Summarize TFs across datasets
my.sum.mph <- Top_TFs_res.mph.5 %>%
  group_by(TF_regulon = source) %>%
  summarise(
    N_Data = n_distinct(Dataset),
    Signs = concat(sign),
    Datasets = concat(unique(Dataset)),
    .groups = "drop"
  )

# 5. Filter recurrent TFs (in more than 12 datasets)
my.rec.mph <- my.sum.mph %>%
  filter(N_Data > 12)

# Step 6: Determine direction consistency
my.rec.mph <- my.rec.mph %>%
  rowwise() %>%
  mutate(
    Direction = {
      tf_signs <- Top_TFs_res.mph.5 %>%
        filter(source == TF_regulon) %>%
        pull(sign)
      if (length(unique(tf_signs)) == 1) "consistent" else "inconsistent"
    }
  ) %>%
  ungroup()


### all marked inconsistent

###############################################################################
# 5. Plot TF enrichment
###############################################################################

# make a table for plotting
my.datasets <- sort(unique(Top_TFs_res.mph$Dataset))

# my.datasets
# [1] "GSE127542_Microglia"       "GSE128830_Peritoneal"      "GSE131869_Microglia_F"     "GSE131869_Microglia_M"    
# [5] "GSE134397_Alveolar"        "GSE134397_Alveolar_CTL"    "GSE137028_Microglia"       "GSE142580_SkM"            
# [9] "GSE145295_Alveolar"        "GSE154832_eWAT"            "GSE156762_Microglia"       "GSE190689_Alveolar"       
# [13] "GSE195507_SkM"             "GSE198666_Callus"          "GSE199763_SkinWound"       "GSE199879_Spleen_Red_Pulp"
# [17] "GSE267529_Microglia_F"     "GSE267529_Microglia_M"     "GSE279654_BMDM"            "GSE93202_Spleen"          
# [21] "GSE93202_VAT"              "GSE98401_Microglia"        "PRJNA1029936"              "PRJNA682234_Callus"   


# grab FDR10 for plotting
mph.top.tf.data <- Top_TFs_res.mph[Top_TFs_res.mph$source %in% my.rec.mph$TF_regulon,]

my.plot.table.mph <- data.frame(matrix(0,nrow(my.rec.mph),length(my.datasets)))
colnames(my.plot.table.mph) <- my.datasets

rownames(my.plot.table.mph) <- my.rec.mph$TF_regulon


for (i in 1:nrow(my.rec.mph)) {
  
  # regulon
  my.reg  <- my.rec.mph$TF_regulon[i]
  
  # regulon data in mph 
  mph.reg <- mph.top.tf.data[mph.top.tf.data$source %in% my.reg, ]
  
  # Deduplicate before making the plotting matrix
  mph.reg <- distinct(mph.reg)
  
  
  for (j in 1:length(my.datasets)) {
    
    # subset for cell type
    mph.reg.ct <- mph.reg[mph.reg$Dataset %in% my.datasets[j],]
    
    # go over mph data
    if (nrow(mph.reg.ct) > 0) {
      if ((mph.reg.ct$score > 0) && (mph.reg.ct$p_value < 0.05))        {
        my.plot.table.mph[i,j]  <- 1
      } else if ((mph.reg.ct$score > 0) && (mph.reg.ct$p_value < 0.1))  {
        my.plot.table.mph[i,j]   <- 0.5
      } else if ((mph.reg.ct$score < 0) && (mph.reg.ct$p_value < 0.05)) {
        my.plot.table.mph[i,j]  <- -1
      } else if ((mph.reg.ct$score < 0) && (mph.reg.ct$p_value < 0.1))  {
        my.plot.table.mph[i,j]   <- -0.5
      }
    }
    
  }
  
}

# Reorder the columns
my.plot.table.mph <- my.plot.table.mph[, desired_alias_order]

# Count number of times TF is up or down
up_count <- rowSums(my.plot.table.mph > 0.5)    # strictly up
down_count <- rowSums(my.plot.table.mph < -0.5) # strictly down

# Combine into a data frame
tf_consistency <- data.frame(up_count, down_count)
tf_consistency

# sort rows by consistency 
consistency_score <- abs(up_count - down_count)

row_order <- order(consistency_score, decreasing = TRUE)

# reorder table and consistency info
my.plot.table.mph <- my.plot.table.mph[row_order, ]
tf_consistency    <- tf_consistency[row_order, ]

write.csv(
  my.plot.table.mph,
  file = file.path(
    my.out.dir,
    paste0(Sys.Date(), "_recurrent_TFs_sorted_by_consistency.csv")
  ),
  row.names = TRUE
)

top_n <- 20
my.plot.table.mph <- my.plot.table.mph[1:top_n, ]
tf_consistency    <- tf_consistency[1:top_n, ]

mph.heat <- Heatmap(my.plot.table.mph, border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
                    cluster_rows = F, cluster_columns = F,
                    column_title = "recurrent decoupleR (macrophage aging)")



# Define a named color mapping
heatmap_colors <- c(
  "-1"   = "darkslateblue",       # FDR 0.05 down
  "-0.5" = "#A9CCE3",  # FDR 0.1 down
  "0"   = "white",      # not significant
  "0.5" = "#F5B7B1",  # FDR 0.1 up
  "1"   = "firebrick3"         # FDR 0.05 up
)

# Convert names to numeric for the scale
col_fun <- colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  heatmap_colors
)

# Build the heatmap with discrete legend
mph.heat <- Heatmap(
  my.plot.table.mph,
  border = TRUE,
  rect_gp = gpar(col = "grey", lwd = 0.5),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  heatmap_legend_param = list(
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("FDR<0.05 down", "FDR<0.1 down", "NS", "FDR<0.1 up", "FDR<0.05 up"),
    title = "Significance",
    nrow = 5,            # one legend item per row
    grid_height = unit(15, "mm"),
    grid_width = unit(6, "mm"),
    row_names_side = "left",
    border = TRUE
  )
)

mph.heat


pdf(file.path(my.out.dir, paste0(Sys.Date(), "_Recurrent_decoupleR_Regulons_byDataset_FDR5_FDR10_overlap.pdf")), width = 8, height = 15)
mph.heat
dev.off()

###############################################################################
# 6. Find recurrent TFs in only Microglia and Alveolar Macrophage sets 
###############################################################################

# Utility function to concatenate vector elements into a comma-separated string
concat <- function(x) paste0(x, collapse = ",")

# 1. Filter significant results
Top_TFs_res.mph.5 <- Top_TFs_res.mph %>%
  filter(p_value < 0.05)

# 2. Deduplicate — keep top scoring TF per dataset
Top_TFs_res.mph.5 <- Top_TFs_res.mph.5 %>%
  group_by(source, Dataset) %>%
  slice_max(abs(score), n = 1, with_ties = FALSE) %>%
  ungroup()

# 3. Add sign column (direction of enrichment)
Top_TFs_res.mph.5 <- Top_TFs_res.mph.5 %>%
  mutate(sign = sign(score))

# 4. extract the microglia datasets and alveolar macrophage datasets separately
microglia_TFs <- Top_TFs_res.mph.5 %>%
  dplyr::filter(Dataset %in% niche_map$Microglia)

alveolar_TFs <- Top_TFs_res.mph.5 %>%
  dplyr::filter(Dataset %in% niche_map$Alveolar)


# 4. Summarize TFs across datasets
my.sum.mph.microglia <- microglia_TFs  %>%
  group_by(TF_regulon = source) %>%
  summarise(
    N_Data = n_distinct(Dataset),
    Signs = concat(sign),
    Datasets = concat(unique(Dataset)),
    .groups = "drop"
  )

my.sum.mph.alveolar <- alveolar_TFs  %>%
  group_by(TF_regulon = source) %>%
  summarise(
    N_Data = n_distinct(Dataset),
    Signs = concat(sign),
    Datasets = concat(unique(Dataset)),
    .groups = "drop"
  )


# 5. Filter recurrent TFs (in more than 12 datasets)
my.rec.mph.microglia <- my.sum.mph.microglia %>%
  filter(N_Data == 8)

my.rec.mph.alveolar <- my.sum.mph.alveolar %>%
  filter(N_Data == 4)

# Step 6: Determine direction consistency
my.rec.mph.microglia  <- my.rec.mph.microglia  %>%
  rowwise() %>%
  mutate(
    Direction = {
      tf_signs <- Top_TFs_res.mph.5 %>%
        filter(source == TF_regulon) %>%
        pull(sign)
      if (length(unique(tf_signs)) == 1) "consistent" else "inconsistent"
    }
  ) %>%
  ungroup()

my.rec.mph.alveolar  <- my.rec.mph.alveolar  %>%
  rowwise() %>%
  mutate(
    Direction = {
      tf_signs <- Top_TFs_res.mph.5 %>%
        filter(source == TF_regulon) %>%
        pull(sign)
      if (length(unique(tf_signs)) == 1) "consistent" else "inconsistent"
    }
  ) %>%
  ungroup()

###############################################################################
# 7. Plot top TF enrichment for microglia and alveolar datasets
###############################################################################

build_tf_heatmap_for_niche <- function(top_tf_res,
                                       recurrent_tf_df,
                                       niche_datasets,
                                       desired_order = NULL,
                                       title = "TF Heatmap") {
  
  # Dataset list for this niche
  my.datasets <- sort(unique(niche_datasets))
  
  # Filter for relevant TFs
  niche.tf.data <- top_tf_res[top_tf_res$source %in% recurrent_tf_df$TF_regulon, ]
  
  # Initialize matrix
  my.plot.table <- data.frame(matrix(0,
                                     nrow = nrow(recurrent_tf_df),
                                     ncol = length(my.datasets)))
  colnames(my.plot.table) <- my.datasets
  rownames(my.plot.table) <- recurrent_tf_df$TF_regulon
  
  # Fill matrix
  for (i in 1:nrow(recurrent_tf_df)) {
    
    tf <- recurrent_tf_df$TF_regulon[i]
    
    tf.sub <- niche.tf.data %>%
      filter(source == tf) %>%
      distinct()
    
    for (j in seq_along(my.datasets)) {
      ds <- my.datasets[j]
      
      tf.ds <- tf.sub %>% filter(Dataset == ds)
      
      if (nrow(tf.ds) > 0) {
        sc <- tf.ds$score
        p  <- tf.ds$p_value
        
        if (sc > 0 && p < 0.05)        my.plot.table[i,j] <- 1
        else if (sc > 0 && p < 0.1)   my.plot.table[i,j] <- 0.5
        else if (sc < 0 && p < 0.05)  my.plot.table[i,j] <- -1
        else if (sc < 0 && p < 0.1)   my.plot.table[i,j] <- -0.5
      }
    }
  }
  
  # Reorder columns if needed
  if (!is.null(desired_order)) {
    my.plot.table <- my.plot.table[, desired_order, drop = FALSE]
  }
  
  # Consistency counts
  up_count   <- rowSums(my.plot.table > 0.5)
  down_count <- rowSums(my.plot.table < -0.5)
  tf_consistency <- data.frame(up_count, down_count)
  
  # Colors
  heatmap_colors <- c(
    "-1"   = "darkslateblue",
    "-0.5" = "#A9CCE3",
    "0"    = "white",
    "0.5"  = "#F5B7B1",
    "1"    = "firebrick3"
  )
  
  col_fun <- colorRamp2(
    c(-1, -0.5, 0, 0.5, 1),
    heatmap_colors
  )
  
  # Build heatmap
  hm <- Heatmap(
    my.plot.table,
    border = TRUE,
    rect_gp = gpar(col = "grey", lwd = 0.5),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = col_fun,
    column_title = title,
    heatmap_legend_param = list(
      at = c(-1, -0.5, 0, 0.5, 1),
      labels = c("FDR<0.05 down", "FDR<0.1 down", "NS",
                 "FDR<0.1 up", "FDR<0.05 up"),
      title = "Significance",
      nrow = 5,
      grid_height = unit(15, "mm"),
      grid_width = unit(6, "mm"),
      border = TRUE
    )
  )
  
  list(
    matrix = my.plot.table,
    consistency = tf_consistency,
    heatmap = hm
  )
}

microglia_datasets <- niche_map$Microglia

res_microglia <- build_tf_heatmap_for_niche(
  top_tf_res = Top_TFs_res.mph.5,
  recurrent_tf_df = my.rec.mph.microglia,
  niche_datasets = microglia_datasets,
  desired_order = microglia_datasets, 
  title = "Microglia TF enrichment (FDR<0.05 recurrent)"
)

pdf(file.path(my.out.dir, paste0(Sys.Date(), "_Recurrent_decoupleR_Regulons_Microglia_FDR5_overlap.pdf")), width = 8, height = 15)
res_microglia$heatmap
dev.off()


alveolar_datasets <- niche_map$Alveolar

res_alveolar <- build_tf_heatmap_for_niche(
  top_tf_res = Top_TFs_res.mph.5,
  recurrent_tf_df = my.rec.mph.alveolar,
  niche_datasets = alveolar_datasets,
  desired_order = alveolar_datasets, 
  title = "Alveolar TF enrichment (FDR<0.05 recurrent)"
)

pdf(file.path(my.out.dir, paste0(Sys.Date(), "_Recurrent_decoupleR_Regulons_Alveolar_FDR5_overlap.pdf")), width = 8, height = 15)
res_alveolar$heatmap
dev.off()

#######################
sink(file = file.path(my.out.dir, paste0(Sys.Date(), "_R_session_Info_decoupleR.txt")))
sessionInfo()
sink()




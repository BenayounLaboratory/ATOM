####################################################################################
# Aging transcriptomics of mouse macrophages
# TF-gene-pathway network
####################################################################################

set.seed(1234)

################################################################################
# 0. Packages
################################################################################

library(dplyr)
library(ggplot2)
library(igraph)
library(ggraph)
library(scales)
library(msigdbr)
library(stringr)
library(tidyr)

################################################################################
# 1. Read in DESeq2 results and build DE gene list
################################################################################

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


padj_cutoff <- 0.15 
# use a more permissive cutoff as the genes contributing to 
# Mitch enrichment are subjected to a specific P value threshold 

de_genes <- unique(na.omit(unlist(lapply(my_deseq_list, function(df) {
  df %>%
    filter(!is.na(padj), padj < padj_cutoff) %>%
    pull(`row.names`)   
}))))


################################################################################
# 2. TF network and pathway enrichment outputs
################################################################################

mm.net <- readRDS("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/NetworkOutput/collectri_mouse_working_copy.rds")

gobp <- read.csv(
  "/Users/ellaschwab/Benayoun_Lee_Local/ATOM/MitchOutput/GOBP_mitch_results_effect.csv",
  stringsAsFactors = FALSE
)

kegg <- read.csv(
  "/Users/ellaschwab/Benayoun_Lee_Local/ATOM/MitchOutput/kegg_mitch_results_effect.csv",
  stringsAsFactors = FALSE
)

reactome <- read.csv(
  "/Users/ellaschwab/Benayoun_Lee_Local/ATOM/MitchOutput/reactome_mitch_results_effect.csv",
  stringsAsFactors = FALSE
)

################################################################################
# 3. Mouse gene set collections
################################################################################

mouse_reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
mouse_kegg     <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG_LEGACY")
mouse_GOBP     <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")

genesets_mouse_reactome <- split(mouse_reactome$gene_symbol, mouse_reactome$gs_name)
genesets_mouse_kegg      <- split(mouse_kegg$gene_symbol, mouse_kegg$gs_description)
genesets_mouse_GOBP      <- split(mouse_GOBP$gene_symbol, mouse_GOBP$gs_name)

################################################################################
# 4. TFs of interest
################################################################################

tf_list <- c(
  "Jun", "Cebpb", "Cebpg", "Sp3", "Spi1", "Cebpd", "Egr1", "Irf1", "Hif1a",
  "Nfkb2", "Rela", "Atf3", "Fos", "Irf2", "Relb", "Usf1", "Bcl6", "Creb1",
  "Jund", "Nfkb1"
)

alveolar_tf_list <- c(
  "Cebpb",
  "Cebpg",
  "E2f1",
  "E2f2",
  "E2f4",
  "Hif1a",
  "Hsf1",
  "Irf1",
  "Jun",
  "Lyl1",
  "Myc",
  "Nfe2l2",
  "Nfkb1",
  "Rela",
  "Sp3",
  "Spi1",
  "Usf1"
)

microglia_tf_list <- c(
  "Atf3",
  "Cebpg",
  "Irf1",
  "Irf2",
  "Jun",
  "Nfkb1",
  "Nr3c1",
  "Rela",
  "Sp3",
  "Stat5a"
)



################################################################################
# 1. Helper: consistent pathways
################################################################################

get_consistent_up_down_heatmaps <- function(
    enrichment_result,
    score_range = c(-1, 1),
    top_n = 10,
    consistency_count = 18
) {
  score_cols <- grep("^s\\.", colnames(enrichment_result), value = TRUE)
  score_cols <- setdiff(score_cols, "s.dist")
  if (length(score_cols) == 0)
    stop("No columns starting with 's.' found in enrichment_result.")
  
  enrichment_result[, score_cols] <- lapply(
    enrichment_result[, score_cols, drop = FALSE],
    function(x) as.numeric(as.character(x))
  )
  
  consistency_sign <- apply(enrichment_result[, score_cols, drop = FALSE], 1, function(x) {
    pos <- sum(x > 0, na.rm = TRUE)
    neg <- sum(x < 0, na.rm = TRUE)
    if (pos >= consistency_count) return(1)
    if (neg >= consistency_count) return(-1)
    return(NA)
  })
  enrichment_result$consistency_sign <- consistency_sign
  
  consistent_up   <- enrichment_result[enrichment_result$consistency_sign ==  1, ]
  consistent_down <- enrichment_result[enrichment_result$consistency_sign == -1, ]
  
  consistent_up <- consistent_up[!is.na(consistent_up$set), ]
  consistent_down <- consistent_down[!is.na(consistent_down$set), ]
  
  top_up   <- head(consistent_up,   n = min(top_n, nrow(consistent_up)))
  top_down <- head(consistent_down, n = min(top_n, nrow(consistent_down)))
  
  message(sprintf(
    "Found %d consistent up and %d consistent down pathways (showing top %d each).",
    nrow(consistent_up), nrow(consistent_down), top_n
  ))
  
  list(
    up = top_up,
    down = top_down,
    up_df = consistent_up,
    down_df = consistent_down
  )
}

################################################################################
# 2. Helper: pathway to gene table
################################################################################

make_pathway_df <- function(sets, direction, genesets, de_genes) {
  genesets[sets] %>%
    stack() %>%
    as.data.frame() %>%
    setNames(c("gene", "pathway")) %>%
    filter(gene %in% de_genes) %>%
    distinct() %>%
    mutate(direction = direction)
}

################################################################################
# 3. Helper: TF-pathway overlaps
################################################################################

make_heat_df <- function(tf_targets, pathway_df) {
  tf_targets %>%
    inner_join(pathway_df, by = "gene", relationship = "many-to-many") %>%
    count(tf, pathway, direction, name = "n_genes")
}

################################################################################
# 4. Helper: make one heatmap per database
################################################################################

make_tf_pathway_heatmap <- function(
    enrichment_result,
    genesets,
    de_genes,
    tf_targets,
    tf_list,  
    out_file,
    top_n = 10,
    consistency_count = 18,
    lim
) {
  hm <- get_consistent_up_down_heatmaps(
    enrichment_result = enrichment_result,
    top_n = top_n,
    consistency_count = consistency_count
  )
  
  top_up_sets   <- hm$up_df$set[seq_len(min(top_n, nrow(hm$up_df)))]
  top_down_sets <- hm$down_df$set[seq_len(min(top_n, nrow(hm$down_df)))]
  
  up_pathway_gene_df <- make_pathway_df(top_up_sets, "UP", genesets, de_genes)
  down_pathway_gene_df <- make_pathway_df(top_down_sets, "DOWN", genesets, de_genes)
  
  heat_df <- bind_rows(
    make_heat_df(tf_targets, up_pathway_gene_df),
    make_heat_df(tf_targets, down_pathway_gene_df)
  ) %>%
    mutate(signed_n_genes = ifelse(direction == "UP", n_genes, -n_genes))
  
  pathway_order <- c(
    up_pathway_gene_df$pathway %>% unique(),
    down_pathway_gene_df$pathway %>% unique()
  )
  
  # TFs on the x-axis in the exact order defined by tf_list regardless
  # of how the rows appear in heat_df
  heat_df <- heat_df %>%
    mutate(
      pathway = factor(pathway, levels = rev(unique(pathway_order))),
      tf      = factor(tf, levels = tf_list)
    )
  
  # Force every TF in tf_list to appear on the x-axis, even TFs that have no
  # overlapping target genes in any of the selected pathways 
  heat_df <- heat_df %>%
    tidyr::complete(
      tf      = factor(tf_list, levels = tf_list),
      pathway,
      fill    = list(signed_n_genes = 0)
    )
  
  p_heat <- ggplot(heat_df, aes(x = tf, y = pathway, fill = signed_n_genes)) +
    geom_tile(color = "black", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-lim, lim),
      breaks = seq(-lim, lim, by = 5),
      name = "number of genes"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed() +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 11),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    labs(x = "TF", y = "Pathway")
  
  ggsave(out_file, plot = p_heat, width = 14, height = 8)
  
  return(list(
    heatmap = p_heat,
    heat_df = heat_df,
    pathway_up = up_pathway_gene_df,
    pathway_down = down_pathway_gene_df,
    consistent = hm
  ))
}

################################################################################
# 5. Prepare TF targets once
################################################################################

mm_df <- as.data.frame(mm.net)

tf_targets <- mm_df %>%
  filter(source %in% tf_list, target %in% de_genes) %>%
  distinct(source, target) %>%
  rename(tf = source, gene = target)

################################################################################
# 6. Loop over GOBP / KEGG / Reactome
################################################################################

setwd("/Users/ellaschwab/Benayoun_Lee_Local/ATOM/NetworkOutput")


pathway_configs <- list(
  GOBP = list(
    enrichment = gobp,
    genesets = genesets_mouse_GOBP,
    outfile = paste0(Sys.Date(), "_TF_pathway_heatmap_GOBP.pdf"),
    lim = 15
  ),
  KEGG = list(
    enrichment = kegg,
    genesets = genesets_mouse_kegg,
    outfile = paste0(Sys.Date(), "_TF_pathway_heatmap_KEGG.pdf"),
    lim = 25
  ),
  Reactome = list(
    enrichment = reactome,
    genesets = genesets_mouse_reactome,
    outfile = paste0(Sys.Date(), "_TF_pathway_heatmap_Reactome.pdf"),
    lim = 15
  )
)


results <- list()

for (nm in names(pathway_configs)) {
  cfg <- pathway_configs[[nm]]
  message("Processing ", nm, " ...")
  
  results[[nm]] <- make_tf_pathway_heatmap(
    enrichment_result = cfg$enrichment,
    genesets = cfg$genesets,
    de_genes = de_genes,
    tf_targets = tf_targets,
    tf_list           = tf_list,  
    out_file = cfg$outfile,
    top_n = 10,
    consistency_count = 18,
    lim = cfg$lim
  )
}

################################################################################
# 7. Additional analyses for tissue-specific TF lists 
################################################################################

alveolar_down   <- read.csv("Alveolar_Down_Shared_Both_Sexes_ORA_results.csv")
alveolar_up     <- read.csv("Alveolar_Up_Shared_Both_Sexes_ORA_results.csv")
microglia_down  <- read.csv("Microglia_Down_Shared_Both_Sexes_ORA_results.csv")
microglia_up    <- read.csv("Microglia_Up_Shared_Both_Sexes_ORA_results.csv")


################################################################################
# Helper: turn one ORA result into a long pathway-gene table
################################################################################

ora_to_pathway_df <- function(ora_df, direction, top_n = 10) {
  ora_df <- ora_df[order(ora_df$p.adjust), ]
  ora_df <- head(ora_df, top_n)
  
  do.call(rbind, lapply(seq_len(nrow(ora_df)), function(i) {
    genes <- strsplit(ora_df$geneID[i], "/", fixed = TRUE)[[1]]
    data.frame(
      pathway   = ora_df$Description[i],
      gene      = genes,
      direction = direction,
      stringsAsFactors = FALSE
    )
  }))
}

################################################################################
# Helper: build TF x pathway heatmap from up + down ORA CSVs
################################################################################

make_tf_pathway_heatmap_from_ora <- function(
    ora_up, ora_down, tf_list, mm_df, out_file, top_n = 10, lim = 15
) {
  up_df   <- ora_to_pathway_df(ora_up,   "UP",   top_n)
  down_df <- ora_to_pathway_df(ora_down, "DOWN", top_n)
  
  all_genes <- unique(c(up_df$gene, down_df$gene))
  
  tf_targets <- mm_df %>%
    filter(source %in% tf_list, target %in% all_genes) %>%
    distinct(source, target) %>%
    rename(tf = source, gene = target)
  
  pathway_df <- bind_rows(up_df, down_df)
  
  heat_df <- tf_targets %>%
    inner_join(pathway_df, by = "gene", relationship = "many-to-many") %>%
    count(tf, pathway, direction, name = "n_genes") %>%
    mutate(signed_n_genes = ifelse(direction == "UP", n_genes, -n_genes))
  
  pathway_order <- c(unique(up_df$pathway), unique(down_df$pathway))
  
  heat_df <- heat_df %>%
    mutate(
      pathway = factor(pathway, levels = rev(pathway_order)),
      tf      = factor(tf, levels = tf_list[tf_list %in% unique(tf)])
    ) %>%
    tidyr::complete(tf, pathway, fill = list(signed_n_genes = 0))
  
  p_heat <- ggplot(heat_df, aes(x = tf, y = pathway, fill = signed_n_genes)) +
    geom_tile(color = "black", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = c(-lim, lim),
      breaks = seq(-lim, lim, by = 5),
      name = "number of genes",
      na.value = "white"   # <-- add this
    )+ 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed() +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid  = element_blank(),
      axis.ticks  = element_blank(),
      axis.title  = element_text(size = 11),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    labs(x = "TF", y = "Pathway")
  
  ggsave(out_file, plot = p_heat, width = 14, height = 8)
  
  invisible(list(heatmap = p_heat, heat_df = heat_df,
                 pathway_up = up_df, pathway_down = down_df,
                 tf_targets = tf_targets))
}

################################################################################
# 8. Run for alveolar and microglia
################################################################################

results_alveolar <- make_tf_pathway_heatmap_from_ora(
  ora_up   = alveolar_up,
  ora_down = alveolar_down,
  tf_list  = alveolar_tf_list,
  mm_df    = mm_df,
  out_file = paste0(Sys.Date(), "_TF_shared_pathway_heatmap_alveolar.pdf"),
  top_n    = 10,
  lim      = 15
)

results_microglia <- make_tf_pathway_heatmap_from_ora(
  ora_up   = microglia_up,
  ora_down = microglia_down,
  tf_list  = microglia_tf_list,
  mm_df    = mm_df,
  out_file = paste0(Sys.Date(), "_TF_shared_pathway_heatmap_microglia.pdf"),
  top_n    = 10,
  lim      = 15
)



#####################################################################################################
# Supplementary: Make TF-gene network for alveolar and microglia mapping back to common aging pathways
#####################################################################################################

# -- Alveolar --
tf_targets_alveolar <- mm_df %>%
  filter(source %in% alveolar_tf_list, target %in% de_genes) %>%
  distinct(source, target) %>%
  rename(tf = source, gene = target)

pathway_configs_alveolar <- list(
  GOBP = list(
    enrichment = gobp,
    genesets   = genesets_mouse_GOBP,
    outfile    = paste0(Sys.Date(), "_TF_pathway_all_heatmap_GOBP_alveolar.pdf"),
    lim        = 15
  ),
  KEGG = list(
    enrichment = kegg,
    genesets   = genesets_mouse_kegg,
    outfile    = paste0(Sys.Date(), "_TF_pathway_all_heatmap_KEGG_alveolar.pdf"),
    lim        = 25
  ),
  Reactome = list(
    enrichment = reactome,
    genesets   = genesets_mouse_reactome,
    outfile    = paste0(Sys.Date(), "_TF_pathway_all_heatmap_Reactome_alveolar.pdf"),
    lim        = 15
  )
)

results_alveolar <- list()

for (nm in names(pathway_configs_alveolar)) {
  cfg <- pathway_configs_alveolar[[nm]]
  message("Processing ", nm, " (alveolar) ...")
  
  results_alveolar[[nm]] <- make_tf_pathway_heatmap(
    enrichment_result = cfg$enrichment,
    genesets          = cfg$genesets,
    de_genes          = de_genes,
    tf_targets        = tf_targets_alveolar,
    tf_list           = alveolar_tf_list,
    out_file          = cfg$outfile,
    top_n             = 10,
    consistency_count = 18,
    lim               = cfg$lim
  )
}

# -- Microglia --
tf_targets_microglia <- mm_df %>%
  filter(source %in% microglia_tf_list, target %in% de_genes) %>%
  distinct(source, target) %>%
  rename(tf = source, gene = target)

pathway_configs_microglia <- list(
  GOBP = list(
    enrichment = gobp,
    genesets   = genesets_mouse_GOBP,
    outfile    = paste0(Sys.Date(), "_TF_pathway_all_heatmap_GOBP_microglia.pdf"),
    lim        = 15
  ),
  KEGG = list(
    enrichment = kegg,
    genesets   = genesets_mouse_kegg,
    outfile    = paste0(Sys.Date(), "_TF_pathway_all_heatmap_KEGG_microglia.pdf"),
    lim        = 25
  ),
  Reactome = list(
    enrichment = reactome,
    genesets   = genesets_mouse_reactome,
    outfile    = paste0(Sys.Date(), "_TF_pathway_all_heatmap_Reactome_microglia.pdf"),
    lim        = 15
  )
)

results_microglia <- list()

for (nm in names(pathway_configs_microglia)) {
  cfg <- pathway_configs_microglia[[nm]]
  message("Processing ", nm, " (microglia) ...")
  
  results_microglia[[nm]] <- make_tf_pathway_heatmap(
    enrichment_result = cfg$enrichment,
    genesets          = cfg$genesets,
    de_genes          = de_genes,
    tf_targets        = tf_targets_microglia,
    tf_list           = microglia_tf_list,
    out_file          = cfg$outfile,
    top_n             = 10,
    consistency_count = 18,
    lim               = cfg$lim
  )
}

################################################################################
# Save package versions
sink(file = paste(Sys.Date(),"TF_Network_Session_Info.txt", sep =""))
sessionInfo()
sink()



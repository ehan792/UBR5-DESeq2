############################################################
# UBR5 KD RNA-seq pipeline
# Runs:
#   1. Human KD vs Control
#   2. Mouse KD vs Control
#
# Input files remain in:
#   Data/HumanKDCounts.xlsx
#   Data/MouseKDCounts.xlsx
#
# Outputs are written to:
#   KD/results/
#   KD/figures/
############################################################

set.seed(1)

#############################################
# Package download
#############################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "DESeq2",
  "apeglm",
  "org.Mm.eg.db",
  "org.Hs.eg.db",
  "clusterProfiler",
  "ReactomePA",
  "msigdbr",
  "fgsea",
  "AnnotationDbi",
  "EnhancedVolcano",
  "pheatmap",
  "enrichplot"
)

cran_pkgs <- c(
  "tidyverse",
  "here",
  "janitor",
  "readxl"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

############################################################
# Package loading
############################################################

library(tidyverse)
library(here)
library(janitor)
library(readxl)

library(DESeq2)
library(apeglm)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(msigdbr)
library(fgsea)
library(AnnotationDbi)
library(EnhancedVolcano)
library(pheatmap)
library(enrichplot)

############################################################
# Project-specific output root
############################################################
# All KD outputs will be written inside:
#   KD/results/
#   KD/figures/
#
# The Data folder remains in the main project folder:
#   Data/HumanKDCounts.xlsx
#   Data/MouseKDCounts.xlsx
############################################################

output_root <- here("KD")

############################################################
# Analysis configuration
############################################################
# Each list entry defines one independent KD analysis.
#
# Important fields:
#   count_file:
#     Input Excel count file.
#
#   contrast_name:
#     Used for folder names and output filenames.
#
#   control_regex / kd_regex:
#     Patterns used to assign samples to Control or KD based on column names.
#
#   orgdb, kegg_organism, reactome_organism:
#     Species-specific databases for GO, KEGG, and Reactome.
#
#   msig_db_species, msig_species, hallmark_collection:
#     Species-specific MSigDB Hallmark settings.
############################################################

analysis_configs <- list(
  human_KD_vs_Control = list(
    analysis_label = "Human UBR5 KD",
    contrast_name = "human_KD_vs_Control",
    count_file = here("Data", "HumanKDCounts.xlsx"),
    control_group = "Control",
    kd_group = "KD",
    control_regex = "jh_2_002",
    kd_regex = "shrubr5",
    orgdb = org.Hs.eg.db,
    kegg_organism = "hsa",
    reactome_organism = "human",
    msig_db_species = "HS",
    msig_species = "Homo sapiens",
    hallmark_collection = "H"
  ),
  
  mouse_KD_vs_Control = list(
    analysis_label = "Mouse UBR5 KD",
    contrast_name = "mouse_KD_vs_Control",
    count_file = here("Data", "MouseKDCounts.xlsx"),
    control_group = "Control",
    kd_group = "KD",
    control_regex = "jw23\\.3",
    kd_regex = "shrubr5",
    orgdb = org.Mm.eg.db,
    kegg_organism = "mmu",
    reactome_organism = "mouse",
    msig_db_species = "MM",
    msig_species = "Mus musculus",
    hallmark_collection = "MH"
  )
)

contrast_names <- names(analysis_configs)

############################################################
# Folder organization
############################################################

output_dirs <- c(
  file.path(output_root, "results"),
  file.path(output_root, "results", "deseq2"),
  file.path(output_root, "results", "gsea"),
  file.path(output_root, "figures"),
  file.path(output_root, "figures", "qc"),
  file.path(output_root, "figures", "volcano"),
  file.path(output_root, "figures", "gsea")
)

for (contrast_name in contrast_names) {
  output_dirs <- c(
    output_dirs,
    file.path(output_root, "results", "deseq2", contrast_name),
    file.path(output_root, "results", "gsea", contrast_name),
    file.path(output_root, "figures", "volcano", contrast_name),
    file.path(output_root, "figures", "gsea", contrast_name)
  )
}

walk(output_dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

############################################################
# Helper functions for output paths
############################################################

get_deseq_results_dir <- function(output_prefix) {
  file.path(output_root, "results", "deseq2", output_prefix)
}

get_gsea_results_dir <- function(output_prefix) {
  file.path(output_root, "results", "gsea", output_prefix)
}

get_volcano_fig_dir <- function(output_prefix) {
  file.path(output_root, "figures", "volcano", output_prefix)
}

get_gsea_fig_dir <- function(output_prefix) {
  file.path(output_root, "figures", "gsea", output_prefix)
}

get_qc_fig_dir <- function() {
  file.path(output_root, "figures", "qc")
}

############################################################
# Data import and preprocessing function
############################################################
# This function reads one KD count file and returns:
#   count_matrix:
#     integer count matrix with genes as rows and samples as columns.
#
#   gene_annot:
#     gene annotation table used later for joining gene names.
#
#   sample_info:
#     DESeq2 metadata table assigning samples to Control or KD.
#
# Important formatting notes for the KD files:
#   - Files are Excel files.
#   - Sample count columns begin with "sample.".
#   - Gene IDs are in "ensembl_gene_id".
#   - Entrez IDs are named "entrezgene", so this script renames
#     that column to "entrezgene_id" for consistency with your
#     CRISPR script and GSEA function.
############################################################

read_kd_count_file <- function(config) {
  
  counts_raw <- readxl::read_excel(config$count_file) %>%
    as.data.frame()
  
  # Standardize Entrez column name to match downstream code.
  if ("entrezgene" %in% colnames(counts_raw) && !"entrezgene_id" %in% colnames(counts_raw)) {
    counts_raw <- counts_raw %>%
      rename(entrezgene_id = entrezgene)
  }
  
  # Annotation columns differ slightly from the CRISPR CSV file.
  # Use only columns actually present in the Excel files.
  possible_annotation_cols <- c(
    "ensembl_gene_id",
    "entrezgene_id",
    "external_gene_name",
    "chromosome_name",
    "start_position",
    "end_position",
    "gene_biotype",
    "external_gene_source",
    "transcript_count",
    "description"
  )
  
  annotation_cols <- intersect(possible_annotation_cols, colnames(counts_raw))
  
  count_cols <- grep("^sample\\.", colnames(counts_raw), value = TRUE)
  
  if (length(count_cols) == 0) {
    stop("No sample count columns found. Expected columns beginning with 'sample.'.")
  }
  
  gene_annot <- counts_raw %>%
    select(all_of(annotation_cols)) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)
  
  count_matrix <- counts_raw %>%
    select(ensembl_gene_id, all_of(count_cols)) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE) %>%
    column_to_rownames("ensembl_gene_id") %>%
    as.matrix()
  
  # DESeq2 requires integer counts.
  # If Excel imports counts as numeric, round safely before integer conversion.
  count_matrix <- round(count_matrix)
  storage.mode(count_matrix) <- "integer"
  
  if (any(is.na(count_matrix))) {
    stop("Count matrix contains NA values.")
  }
  
  if (any(count_matrix < 0)) {
    stop("Count matrix contains negative counts.")
  }
  
  sample_info <- tibble(
    sample = colnames(count_matrix),
    condition = case_when(
      str_detect(sample, config$kd_regex) ~ config$kd_group,
      str_detect(sample, config$control_regex) ~ config$control_group,
      TRUE ~ NA_character_
    )
  ) %>%
    mutate(
      condition = factor(
        condition,
        levels = c(config$control_group, config$kd_group)
      )
    ) %>%
    column_to_rownames("sample")
  
  stopifnot(all(rownames(sample_info) == colnames(count_matrix)))
  stopifnot(!any(is.na(sample_info$condition)))
  
  return(list(
    counts_raw = counts_raw,
    gene_annot = gene_annot,
    count_matrix = count_matrix,
    sample_info = sample_info
  ))
}

############################################################
# QC function
############################################################
# This creates:
#   - normalized count table for all 6 samples
#   - PCA plot
#   - grouped sample distance heatmap
#   - clustered sample distance heatmap
#
# Filtering for QC:
#   keep genes with >=10 counts in at least 3 of the 6 samples.
#
# This matches the same basic low-count threshold used for DESeq2.
############################################################

run_qc <- function(dds_all, sample_info, output_prefix, analysis_label) {
  
  keep_qc <- rowSums(counts(dds_all) >= 10) >= 3
  dds_qc <- dds_all[keep_qc, ]
  
  dds_qc <- DESeq(dds_qc)
  
  norm_counts_qc <- counts(dds_qc, normalized = TRUE)
  
  write.csv(
    as.data.frame(norm_counts_qc) %>%
      rownames_to_column("ensembl_gene_id"),
    file.path(
      get_deseq_results_dir(output_prefix),
      paste0(output_prefix, "_normalized_counts_QC_all_samples.csv")
    ),
    row.names = FALSE
  )
  
  clean_sample_names <- setNames(
    paste0(
      sample_info$condition,
      " ",
      ave(
        seq_along(sample_info$condition),
        sample_info$condition,
        FUN = seq_along
      )
    ),
    rownames(sample_info)
  )
  
  vsd <- vst(dds_qc, blind = FALSE)
  
  pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  pca_data$clean_name <- clean_sample_names[pca_data$name]
  
  percent_var <- round(100 * attr(pca_data, "percentVar"), 2)
  
  p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = clean_name)) +
    geom_point(size = 4) +
    geom_text(vjust = -1, size = 3, show.legend = FALSE) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggtitle(paste0(analysis_label, ": PCA using VST counts")) +
    theme_bw()
  
  ggsave(
    file.path(get_qc_fig_dir(), paste0(output_prefix, "_PCA_vst_condition.png")),
    p_pca,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  sample_dists <- dist(t(assay(vsd)))
  sample_dist_mat <- as.matrix(sample_dists)
  
  rownames(sample_dist_mat) <- colnames(vsd)
  colnames(sample_dist_mat) <- colnames(vsd)
  
  sample_order_df <- sample_info %>%
    rownames_to_column("sample") %>%
    arrange(condition, sample)
  
  sample_order <- sample_order_df$sample
  
  sample_dist_mat_ordered <- sample_dist_mat[sample_order, sample_order]
  
  rownames(sample_dist_mat_ordered) <- clean_sample_names[rownames(sample_dist_mat_ordered)]
  colnames(sample_dist_mat_ordered) <- clean_sample_names[colnames(sample_dist_mat_ordered)]
  
  group_sizes <- table(sample_order_df$condition)
  group_gaps <- cumsum(as.integer(group_sizes))
  group_gaps <- group_gaps[-length(group_gaps)]
  
  png(
    file.path(get_qc_fig_dir(), paste0(output_prefix, "_sample_distance_heatmap_grouped.png")),
    width = 1800,
    height = 1600,
    res = 250
  )
  
  pheatmap(
    sample_dist_mat_ordered,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_row = group_gaps,
    gaps_col = group_gaps,
    main = paste0(analysis_label, ": sample distances, grouped by condition"),
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 45
  )
  
  dev.off()
  
  png(
    file.path(get_qc_fig_dir(), paste0(output_prefix, "_sample_distance_heatmap_clustered.png")),
    width = 1800,
    height = 1600,
    res = 250
  )
  
  pheatmap(
    sample_dist_mat_ordered,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = paste0(analysis_label, ": sample distances, clustered"),
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 45
  )
  
  dev.off()
  
  return(list(
    dds_qc = dds_qc,
    vsd = vsd,
    pca_data = pca_data
  ))
}

############################################################
# Contrast-specific DESeq2 analysis function
############################################################
# This is the KD equivalent of your CRISPR contrast function.
#
# Here there is one comparison per dataset:
#   KD vs Control
#
# Important filtering parameters:
#   min_count = 10:
#     A sample counts as expressing a gene if raw count >= 10.
#
#   min_samples = 3:
#     A gene is retained if at least 3 of the 6 samples pass min_count.
#
# Because each KD dataset has 3 KD and 3 Control samples, this means
# a gene can pass if it is consistently expressed in either phenotype.
############################################################

run_deseq_for_kd <- function(
    dds_all,
    output_prefix,
    gene_annot,
    treatment_group = "KD",
    control_group = "Control",
    min_count = 10,
    min_samples = 3
) {
  
  dds_contrast <- dds_all
  
  dds_contrast$condition <- droplevels(dds_contrast$condition)
  dds_contrast$condition <- relevel(dds_contrast$condition, ref = control_group)
  
  keep <- rowSums(counts(dds_contrast) >= min_count) >= min_samples
  dds_contrast <- dds_contrast[keep, ]
  
  message(output_prefix, ": kept ", sum(keep), " genes after filtering.")
  
  filter_summary <- tibble(
    contrast = output_prefix,
    treatment_group = treatment_group,
    control_group = control_group,
    samples_used = paste(colnames(dds_contrast), collapse = ";"),
    number_of_samples_used = ncol(dds_contrast),
    min_count = min_count,
    min_samples = min_samples,
    genes_before_filtering = length(keep),
    genes_after_filtering = sum(keep),
    genes_removed_by_filtering = length(keep) - sum(keep)
  )
  
  write.csv(
    filter_summary,
    file.path(
      get_deseq_results_dir(output_prefix),
      paste0(output_prefix, "_filter_summary.csv")
    ),
    row.names = FALSE
  )
  
  dds_contrast <- DESeq(dds_contrast)
  
  norm_counts <- counts(dds_contrast, normalized = TRUE)
  
  write.csv(
    as.data.frame(norm_counts) %>%
      rownames_to_column("ensembl_gene_id"),
    file.path(
      get_deseq_results_dir(output_prefix),
      paste0(output_prefix, "_normalized_counts.csv")
    ),
    row.names = FALSE
  )
  
  res <- results(
    dds_contrast,
    contrast = c("condition", treatment_group, control_group),
    alpha = 0.05
  )
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(gene_annot, by = "ensembl_gene_id") %>%
    arrange(padj)
  
  write.csv(
    res_df,
    file.path(
      get_deseq_results_dir(output_prefix),
      paste0(output_prefix, "_DESeq2_results.csv")
    ),
    row.names = FALSE
  )
  
  message("Available coefficients for ", output_prefix, ":")
  print(resultsNames(dds_contrast))
  
  coef_name <- paste0("condition_", treatment_group, "_vs_", control_group)
  
  if (!coef_name %in% resultsNames(dds_contrast)) {
    stop(
      "Coefficient not found: ", coef_name,
      "\nAvailable coefficients are: ",
      paste(resultsNames(dds_contrast), collapse = ", ")
    )
  }
  
  res_shrunk <- lfcShrink(
    dds_contrast,
    coef = coef_name,
    type = "apeglm"
  )
  
  res_shrunk_df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(gene_annot, by = "ensembl_gene_id") %>%
    arrange(padj)
  
  write.csv(
    res_shrunk_df,
    file.path(
      get_deseq_results_dir(output_prefix),
      paste0(output_prefix, "_DESeq2_shrunkLFC.csv")
    ),
    row.names = FALSE
  )
  
  summary_df <- tibble(
    contrast = output_prefix,
    comparison = paste(treatment_group, "vs", control_group),
    min_count_filter = min_count,
    min_samples_filter = min_samples,
    genes_tested_after_filtering = nrow(res_df),
    significant_padj = sum(res_df$padj < 0.05, na.rm = TRUE),
    up_padj = sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE),
    down_padj = sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE),
    up_padj_lfc1 = sum(res_df$padj < 0.05 & res_df$log2FoldChange >= 1, na.rm = TRUE),
    down_padj_lfc1 = sum(res_df$padj < 0.05 & res_df$log2FoldChange <= -1, na.rm = TRUE)
  )
  
  write.csv(
    summary_df,
    file.path(
      get_deseq_results_dir(output_prefix),
      paste0(output_prefix, "_DE_summary.csv")
    ),
    row.names = FALSE
  )
  
  return(list(
    dds = dds_contrast,
    res = res_df,
    shrunk = res_shrunk_df,
    summary = summary_df,
    filter_summary = filter_summary
  ))
}

############################################################
# Volcano plot function
############################################################

plot_volcano <- function(
    res_shrunk_df,
    output_prefix,
    title,
    pCutoff = 0.05,
    FCcutoff = 1
) {
  
  png(
    file.path(
      get_volcano_fig_dir(output_prefix),
      paste0("volcano_", output_prefix, ".png")
    ),
    width = 2000,
    height = 1800,
    res = 250
  )
  
  print(
    EnhancedVolcano(
      res_shrunk_df,
      lab = res_shrunk_df$external_gene_name,
      x = "log2FoldChange",
      y = "padj",
      title = title,
      subtitle = "DESeq2 with apeglm-shrunk log2FC",
      pCutoff = pCutoff,
      FCcutoff = FCcutoff
    )
  )
  
  dev.off()
}

############################################################
# Prepare ranked gene list for GSEA
############################################################
# Uses the unshrunk DESeq2 Wald statistic.
# Positive statistic means higher in KD than Control.
############################################################

make_ranked_list <- function(res_df) {
  
  ranked <- res_df %>%
    filter(
      !is.na(stat),
      is.finite(stat),
      !is.na(entrezgene_id),
      entrezgene_id != "",
      entrezgene_id != "NA"
    ) %>%
    mutate(entrezgene_id = as.character(entrezgene_id)) %>%
    group_by(entrezgene_id) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(stat))
  
  gene_list <- ranked$stat
  names(gene_list) <- ranked$entrezgene_id
  
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gene_list <- gene_list[
    !is.na(gene_list) &
      is.finite(gene_list) &
      !is.na(names(gene_list)) &
      names(gene_list) != "" &
      names(gene_list) != "NA"
  ]
  
  gene_list <- gene_list[!duplicated(names(gene_list))]
  
  return(gene_list)
}

############################################################
# GSEA functions
############################################################

get_msig_entrez_column <- function(msig_df) {
  if ("ncbi_gene" %in% colnames(msig_df)) {
    return("ncbi_gene")
  } else if ("entrez_gene" %in% colnames(msig_df)) {
    return("entrez_gene")
  } else {
    stop("Could not find Entrez gene column in msigdbr output.")
  }
}

run_hallmark_gsea <- function(gene_list, output_prefix, config) {
  
  hallmark_sets <- msigdbr(
    db_species = config$msig_db_species,
    species = config$msig_species,
    collection = config$hallmark_collection
  )
  
  entrez_col <- get_msig_entrez_column(hallmark_sets)
  
  hallmark_sets <- hallmark_sets %>%
    dplyr::select(gs_name, all_of(entrez_col)) %>%
    dplyr::filter(!is.na(.data[[entrez_col]]))
  
  pathways <- split(
    hallmark_sets[[entrez_col]],
    hallmark_sets$gs_name
  )
  
  fgsea_res <- fgsea(
    pathways = pathways,
    stats = gene_list,
    minSize = 15,
    maxSize = 500
  ) %>%
    dplyr::arrange(padj)
  
  fgsea_export <- fgsea_res %>%
    dplyr::mutate(
      leadingEdge = sapply(leadingEdge, paste, collapse = ";")
    )
  
  write.csv(
    fgsea_export,
    file.path(
      get_gsea_results_dir(output_prefix),
      paste0(output_prefix, "_GSEA_Hallmark.csv")
    ),
    row.names = FALSE
  )
  
  return(list(
    results = fgsea_res,
    export = fgsea_export,
    pathways = pathways
  ))
}

run_kegg_gsea <- function(gene_list, output_prefix, config) {
  
  kegg_res <- gseKEGG(
    geneList = gene_list,
    organism = config$kegg_organism,
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = FALSE
  )
  
  kegg_df <- as.data.frame(kegg_res) %>%
    arrange(p.adjust)
  
  write.csv(
    kegg_df,
    file.path(
      get_gsea_results_dir(output_prefix),
      paste0(output_prefix, "_GSEA_KEGG.csv")
    ),
    row.names = FALSE
  )
  
  return(kegg_res)
}

run_reactome_gsea <- function(gene_list, output_prefix, config) {
  
  reactome_res <- gsePathway(
    geneList = gene_list,
    organism = config$reactome_organism,
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = FALSE
  )
  
  reactome_df <- as.data.frame(reactome_res) %>%
    arrange(p.adjust)
  
  write.csv(
    reactome_df,
    file.path(
      get_gsea_results_dir(output_prefix),
      paste0(output_prefix, "_GSEA_Reactome.csv")
    ),
    row.names = FALSE
  )
  
  return(reactome_res)
}

run_go_bp_gsea <- function(gene_list, output_prefix, config) {
  
  go_res <- gseGO(
    geneList = gene_list,
    OrgDb = config$orgdb,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = FALSE,
    nPermSimple = 10000
  )
  
  go_df <- as.data.frame(go_res) %>%
    filter(
      !is.na(ID),
      !is.na(pvalue),
      !is.na(p.adjust),
      !is.na(NES),
      is.finite(pvalue),
      is.finite(p.adjust),
      is.finite(NES)
    ) %>%
    arrange(p.adjust)
  
  write.csv(
    go_df,
    file.path(
      get_gsea_results_dir(output_prefix),
      paste0(output_prefix, "_GSEA_GO_BP.csv")
    ),
    row.names = FALSE
  )
  
  return(go_res)
}

############################################################
# GSEA plotting functions
############################################################

plot_fgsea_dotplot <- function(fgsea_df, output_prefix, title, top_n = 20) {
  
  plot_df <- fgsea_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    slice_head(n = top_n) %>%
    mutate(
      pathway = factor(pathway, levels = rev(pathway)),
      direction = ifelse(NES > 0, "Up", "Down")
    )
  
  if (nrow(plot_df) == 0) {
    message("No Hallmark GSEA results to plot for: ", output_prefix)
    return(NULL)
  }
  
  p <- ggplot(plot_df, aes(x = NES, y = pathway, size = size, color = padj)) +
    geom_point() +
    theme_bw() +
    labs(
      title = title,
      x = "Normalized enrichment score, NES",
      y = NULL,
      size = "Gene set size",
      color = "Adjusted p-value"
    )
  
  ggsave(
    file.path(
      get_gsea_fig_dir(output_prefix),
      paste0(output_prefix, "_Hallmark_dotplot.png")
    ),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

plot_top_fgsea_enrichment <- function(gsea_obj, gene_list, output_prefix) {
  
  fgsea_df <- gsea_obj$export
  pathways <- gsea_obj$pathways
  
  top_up <- fgsea_df %>%
    filter(!is.na(padj), NES > 0) %>%
    arrange(padj) %>%
    slice_head(n = 1) %>%
    pull(pathway)
  
  top_down <- fgsea_df %>%
    filter(!is.na(padj), NES < 0) %>%
    arrange(padj) %>%
    slice_head(n = 1) %>%
    pull(pathway)
  
  if (length(top_up) == 1) {
    p_up <- plotEnrichment(
      pathways[[top_up]],
      gene_list
    ) +
      labs(title = paste0(output_prefix, ": ", top_up))
    
    ggsave(
      file.path(
        get_gsea_fig_dir(output_prefix),
        paste0(output_prefix, "_Hallmark_top_up_enrichment.png")
      ),
      p_up,
      width = 7,
      height = 5,
      dpi = 300
    )
  }
  
  if (length(top_down) == 1) {
    p_down <- plotEnrichment(
      pathways[[top_down]],
      gene_list
    ) +
      labs(title = paste0(output_prefix, ": ", top_down))
    
    ggsave(
      file.path(
        get_gsea_fig_dir(output_prefix),
        paste0(output_prefix, "_Hallmark_top_down_enrichment.png")
      ),
      p_down,
      width = 7,
      height = 5,
      dpi = 300
    )
  }
}

save_gsea_dotplot <- function(
    gsea_result,
    output_prefix,
    database_name,
    title,
    show_n = 10,
    label_width = 40,
    fig_width = 11,
    fig_height = 7,
    split_direction = TRUE,
    text_size = 8
) {
  
  gsea_df <- as.data.frame(gsea_result)
  
  if (nrow(gsea_df) == 0) {
    message("No GSEA results to plot for: ", output_prefix, " ", database_name)
    return(NULL)
  }
  
  if (split_direction) {
    p <- enrichplot::dotplot(
      gsea_result,
      showCategory = show_n,
      split = ".sign",
      label_format = label_width
    ) +
      facet_grid(. ~ .sign)
  } else {
    p <- enrichplot::dotplot(
      gsea_result,
      showCategory = show_n,
      label_format = label_width
    )
  }
  
  p <- p +
    ggtitle(title) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = text_size),
      axis.text.x = element_text(size = 9),
      plot.title = element_text(size = 14),
      strip.text = element_text(size = 11)
    )
  
  ggsave(
    file.path(
      get_gsea_fig_dir(output_prefix),
      paste0(output_prefix, "_", database_name, "_dotplot.png")
    ),
    p,
    width = fig_width,
    height = fig_height,
    dpi = 300
  )
  
  return(p)
}

save_top_gseaplot <- function(gsea_result, output_prefix, database_name, title_prefix) {
  
  gsea_df <- as.data.frame(gsea_result) %>%
    filter(
      !is.na(ID),
      !is.na(p.adjust),
      is.finite(p.adjust)
    )
  
  if (nrow(gsea_df) == 0) {
    message("No valid GSEA results to plot for: ", output_prefix, " ", database_name)
    return(NULL)
  }
  
  top_id <- gsea_df %>%
    arrange(p.adjust) %>%
    slice_head(n = 1) %>%
    pull(ID)
  
  p <- enrichplot::gseaplot2(
    gsea_result,
    geneSetID = top_id,
    title = paste0(title_prefix, ": ", top_id)
  )
  
  ggsave(
    file.path(
      get_gsea_fig_dir(output_prefix),
      paste0(output_prefix, "_", database_name, "_top_gseaplot.png")
    ),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

############################################################
# Run full KD pipeline for one dataset
############################################################

run_kd_pipeline <- function(config) {
  
  output_prefix <- config$contrast_name
  
  message("Starting analysis: ", config$analysis_label)
  
  imported <- read_kd_count_file(config)
  
  dds_all <- DESeqDataSetFromMatrix(
    countData = imported$count_matrix,
    colData = imported$sample_info,
    design = ~ condition
  )
  
  qc_results <- run_qc(
    dds_all = dds_all,
    sample_info = imported$sample_info,
    output_prefix = output_prefix,
    analysis_label = config$analysis_label
  )
  
  deseq_results <- run_deseq_for_kd(
    dds_all = dds_all,
    output_prefix = output_prefix,
    gene_annot = imported$gene_annot,
    treatment_group = config$kd_group,
    control_group = config$control_group,
    min_count = 10,
    min_samples = 3
  )
  
  plot_volcano(
    res_shrunk_df = deseq_results$shrunk,
    output_prefix = output_prefix,
    title = paste0(config$analysis_label, ": KD vs Control"),
    pCutoff = 0.05,
    FCcutoff = 1
  )
  
  gene_list <- make_ranked_list(deseq_results$res)
  
  gsea_hallmark <- run_hallmark_gsea(gene_list, output_prefix, config)
  gsea_kegg <- run_kegg_gsea(gene_list, output_prefix, config)
  gsea_reactome <- run_reactome_gsea(gene_list, output_prefix, config)
  gsea_go <- run_go_bp_gsea(gene_list, output_prefix, config)
  
  plot_fgsea_dotplot(
    gsea_hallmark$export,
    output_prefix,
    paste0("Hallmark GSEA: ", config$analysis_label)
  )
  
  plot_top_fgsea_enrichment(
    gsea_hallmark,
    gene_list,
    output_prefix
  )
  
  gsea_plot_settings <- list(
    KEGG = list(
      result = gsea_kegg,
      title = paste0("KEGG GSEA: ", config$analysis_label),
      show_n = 10,
      label_width = 40,
      fig_width = 11,
      fig_height = 7,
      text_size = 8
    ),
    Reactome = list(
      result = gsea_reactome,
      title = paste0("Reactome GSEA: ", config$analysis_label),
      show_n = 10,
      label_width = 45,
      fig_width = 13,
      fig_height = 7,
      text_size = 8
    ),
    GO_BP = list(
      result = gsea_go,
      title = paste0("GO Biological Process GSEA: ", config$analysis_label),
      show_n = 8,
      label_width = 35,
      fig_width = 14,
      fig_height = 8,
      text_size = 7
    )
  )
  
  for (database_name in names(gsea_plot_settings)) {
    
    settings <- gsea_plot_settings[[database_name]]
    
    save_gsea_dotplot(
      gsea_result = settings$result,
      output_prefix = output_prefix,
      database_name = database_name,
      title = settings$title,
      show_n = settings$show_n,
      label_width = settings$label_width,
      fig_width = settings$fig_width,
      fig_height = settings$fig_height,
      text_size = settings$text_size
    )
    
    save_top_gseaplot(
      gsea_result = settings$result,
      output_prefix = output_prefix,
      database_name = database_name,
      title_prefix = paste0("Top ", database_name, " GSEA: ", config$analysis_label)
    )
  }
  
  return(list(
    imported = imported,
    qc = qc_results,
    deseq = deseq_results,
    gene_list = gene_list,
    gsea = list(
      hallmark = gsea_hallmark,
      kegg = gsea_kegg,
      reactome = gsea_reactome,
      go_bp = gsea_go
    )
  ))
}

############################################################
# Run both KD analyses
############################################################

human_KD_results <- run_kd_pipeline(analysis_configs$human_KD_vs_Control)
mouse_KD_results <- run_kd_pipeline(analysis_configs$mouse_KD_vs_Control)

############################################################
# Save session info for reproducibility
############################################################

writeLines(
  capture.output(sessionInfo()),
  file.path(output_root, "results", "sessionInfo.txt")
)
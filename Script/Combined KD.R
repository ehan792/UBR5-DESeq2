############################################################
# RNA-seq DESeq2 + GSEA pipeline for Mouse UBR5 KD and Human UBR5 KD
#
# Input files expected:
#   Data/MouseKDCounts.xlsx
#   Data/HumanKDCounts.xlsx
#
# Output folders created:
#   Mouse_UBR5_KD/Results/
#   Mouse_UBR5_KD/Figures/
#   Human_UBR5_KD/Results/
#   Human_UBR5_KD/Figures/
#
# Main contrast:
#   KD vs Control
############################################################

set.seed(1)

############################################################
# 0. Package install/load
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_pkgs <- c(
  "tidyverse",
  "here",
  "janitor",
  "readxl"
)

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
# 1. Experiment configurations
############################################################

analysis_configs <- list(
  
  mouse_UBR5_KD = list(
    analysis_name = "mouse_UBR5_KD",
    root_folder = here("Mouse_UBR5_KD"),
    count_file = here("Data", "MouseKDCounts.xlsx"),
    control_group = "Control",
    treatment_group = "KD",
    group_levels = c("Control", "KD"),
    
    # Mouse file sample names:
    # KD:      sample.shrubr5_1_jw, sample.shrubr5_2_jw, sample.shrubr5_3_jw
    # Control: sample.jw23.3_1, sample.jw23.3_2, sample.jw23.3_3
    control_regex = "jw23\\.3",
    treatment_regex = "shrubr5",
    
    orgdb = org.Mm.eg.db,
    kegg_organism = "mmu",
    reactome_organism = "mouse",
    msig_db_species = "MM",
    msig_species = "Mus musculus",
    hallmark_collection = "MH"
  ),
  
  human_UBR5_KD = list(
    analysis_name = "human_UBR5_KD",
    root_folder = here("Human_UBR5_KD"),
    count_file = here("Data", "HumanKDCounts.xlsx"),
    control_group = "Control",
    treatment_group = "KD",
    group_levels = c("Control", "KD"),
    
    # Human file sample names:
    # KD:      sample.shrubr5_1_2002, sample.shrubr5_2_2002, sample.shrubr5_3_2002
    # Control: sample.jh_2_002_1, sample.jh_2_002_2, sample.jh_2_002_3
    control_regex = "jh_2_002",
    treatment_regex = "shrubr5",
    
    orgdb = org.Hs.eg.db,
    kegg_organism = "hsa",
    reactome_organism = "human",
    msig_db_species = "HS",
    msig_species = "Homo sapiens",
    hallmark_collection = "H"
  )
)

############################################################
# 2. Folder helper functions
############################################################

make_output_dirs <- function(config) {
  
  output_dirs <- c(
    config$root_folder,
    
    file.path(config$root_folder, "Results"),
    file.path(config$root_folder, "Results", "QC"),
    file.path(config$root_folder, "Results", "DESeq2"),
    file.path(config$root_folder, "Results", "GSEA"),
    
    file.path(config$root_folder, "Figures"),
    file.path(config$root_folder, "Figures", "QC"),
    file.path(config$root_folder, "Figures", "Volcano"),
    file.path(config$root_folder, "Figures", "GSEA")
  )
  
  walk(output_dirs, dir.create, showWarnings = FALSE, recursive = TRUE)
}

get_deseq_results_dir <- function(config) {
  file.path(config$root_folder, "Results", "DESeq2")
}

get_gsea_results_dir <- function(config) {
  file.path(config$root_folder, "Results", "GSEA")
}

get_qc_results_dir <- function(config) {
  file.path(config$root_folder, "Results", "QC")
}

get_volcano_fig_dir <- function(config) {
  file.path(config$root_folder, "Figures", "Volcano")
}

get_gsea_fig_dir <- function(config) {
  file.path(config$root_folder, "Figures", "GSEA")
}

get_qc_fig_dir <- function(config) {
  file.path(config$root_folder, "Figures", "QC")
}
############################################################
# 3. Input parsing and sample metadata
############################################################

read_count_file <- function(count_file) {
  
  counts_raw <- readxl::read_excel(
    count_file,
    sheet = 1,
    .name_repair = "minimal"
  ) %>%
    as.data.frame()
  
  return(counts_raw)
}

prepare_count_data <- function(counts_raw) {
  
  # New KD files use "entrezgene"; standardize to "entrezgene_id".
  if ("entrezgene" %in% colnames(counts_raw) && !"entrezgene_id" %in% colnames(counts_raw)) {
    counts_raw <- counts_raw %>%
      rename(entrezgene_id = entrezgene)
  }
  
  required_base_cols <- c(
    "ensembl_gene_id",
    "entrezgene_id",
    "external_gene_name",
    "gene_biotype",
    "description"
  )
  
  missing_required <- setdiff(required_base_cols, colnames(counts_raw))
  
  if (length(missing_required) > 0) {
    stop(
      "Missing required annotation columns: ",
      paste(missing_required, collapse = ", ")
    )
  }
  
  count_cols <- grep("^sample\\.", colnames(counts_raw), value = TRUE)
  
  if (length(count_cols) == 0) {
    stop("No sample count columns found. Expected columns beginning with 'sample.'.")
  }
  
  annotation_cols <- intersect(
    c(
      "ensembl_gene_id",
      "entrezgene_id",
      "external_gene_name",
      "gene_biotype",
      "external_gene_source",
      "transcript_count",
      "description",
      "chromosome_name",
      "start_position",
      "end_position"
    ),
    colnames(counts_raw)
  )
  
  gene_annot <- counts_raw %>%
    select(all_of(annotation_cols)) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)
  
  count_matrix <- counts_raw %>%
    select(ensembl_gene_id, all_of(count_cols)) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE) %>%
    column_to_rownames("ensembl_gene_id") %>%
    as.matrix()
  
  # Excel imports counts as numeric; DESeq2 needs integer-like counts.
  count_matrix <- round(count_matrix)
  storage.mode(count_matrix) <- "integer"
  
  if (any(is.na(count_matrix))) {
    stop("NA values detected in count matrix after import.")
  }
  
  if (any(count_matrix < 0)) {
    stop("Negative values detected in count matrix. Counts must be non-negative.")
  }
  
  list(
    gene_annot = gene_annot,
    count_matrix = count_matrix,
    count_cols = count_cols
  )
}

make_sample_info <- function(count_matrix, config) {
  
  sample_info <- tibble(
    sample = colnames(count_matrix),
    condition = case_when(
      str_detect(sample, regex(config$treatment_regex, ignore_case = TRUE)) ~ config$treatment_group,
      str_detect(sample, regex(config$control_regex, ignore_case = TRUE)) ~ config$control_group,
      TRUE ~ NA_character_
    )
  ) %>%
    mutate(
      condition = factor(condition, levels = config$group_levels)
    ) %>%
    column_to_rownames("sample")
  
  if (any(is.na(sample_info$condition))) {
    print(sample_info)
    stop("Some samples could not be assigned to a condition. Check control_regex and treatment_regex.")
  }
  
  stopifnot(all(rownames(sample_info) == colnames(count_matrix)))
  
  return(sample_info)
}

make_clean_sample_names <- function(sample_info) {
  
  clean_names_df <- sample_info %>%
    rownames_to_column("sample") %>%
    group_by(condition) %>%
    arrange(sample, .by_group = TRUE) %>%
    mutate(clean_name = paste0(condition, " ", row_number())) %>%
    ungroup()
  
  clean_sample_names <- setNames(
    clean_names_df$clean_name,
    clean_names_df$sample
  )
  
  return(clean_sample_names)
}

############################################################
# 4. QC diagnostics
############################################################

write_basic_qc_tables <- function(count_matrix, sample_info, config) {
  
  library_sizes <- colSums(count_matrix)
  detected_genes_10 <- colSums(count_matrix >= 10)
  detected_genes_1 <- colSums(count_matrix >= 1)
  
  qc_df <- tibble(
    sample = colnames(count_matrix),
    clean_condition = as.character(sample_info[colnames(count_matrix), "condition"]),
    library_size = library_sizes,
    detected_genes_count_ge_1 = detected_genes_1,
    detected_genes_count_ge_10 = detected_genes_10
  )
  
  write.csv(
    qc_df,
    file.path(get_qc_results_dir(config), "library_size_detected_genes.csv"),
    row.names = FALSE
  )
  
  return(qc_df)
}

plot_pca_custom <- function(vsd, sample_info, clean_sample_names, config, ntop = 500) {
  
  mat <- assay(vsd)
  
  # Use top variable genes, like DESeq2::plotPCA, but keep more precise variance labels.
  rv <- rowVars(mat)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  pca <- prcomp(t(mat[select, ]), center = TRUE, scale. = FALSE)
  
  percent_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  
  pca_data <- as.data.frame(pca$x[, 1:2]) %>%
    rownames_to_column("sample") %>%
    left_join(
      sample_info %>%
        rownames_to_column("sample"),
      by = "sample"
    ) %>%
    mutate(clean_name = clean_sample_names[sample])
  
  p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = clean_name)) +
    geom_point(size = 4) +
    geom_text(vjust = -1, size = 3, show.legend = FALSE) +
    xlab(paste0("PC1: ", round(percent_var[1], 2), "% variance")) +
    ylab(paste0("PC2: ", round(percent_var[2], 2), "% variance")) +
    ggtitle(paste0("PCA using VST counts: ", config$analysis_name)) +
    theme_bw()
  
  ggsave(
    file.path(get_qc_fig_dir(config), "PCA_vst_condition.png"),
    p_pca,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  pca_variance_df <- tibble(
    PC = paste0("PC", seq_along(percent_var)),
    percent_variance = percent_var
  )
  
  write.csv(
    pca_variance_df,
    file.path(get_qc_results_dir(config), "PCA_percent_variance.csv"),
    row.names = FALSE
  )
  
  if (percent_var[1] > 99) {
    warning(
      config$analysis_name,
      ": PC1 explains >99% variance. This may be real strong group separation, ",
      "but check sample labels, library sizes, and sample distance/correlation matrices."
    )
  }
  
  return(p_pca)
}

plot_sample_distance_heatmaps <- function(vsd, sample_info, clean_sample_names, config) {
  
  mat <- assay(vsd)
  
  sample_dists <- dist(t(mat))
  sample_dist_mat <- as.matrix(sample_dists)
  
  rownames(sample_dist_mat) <- colnames(vsd)
  colnames(sample_dist_mat) <- colnames(vsd)
  
  sample_corr_mat <- cor(mat, method = "pearson")
  
  sample_order_df <- sample_info %>%
    rownames_to_column("sample") %>%
    arrange(condition, sample)
  
  sample_order <- sample_order_df$sample
  
  sample_dist_mat_ordered <- sample_dist_mat[sample_order, sample_order]
  sample_corr_mat_ordered <- sample_corr_mat[sample_order, sample_order]
  
  rownames(sample_dist_mat_ordered) <- clean_sample_names[rownames(sample_dist_mat_ordered)]
  colnames(sample_dist_mat_ordered) <- clean_sample_names[colnames(sample_dist_mat_ordered)]
  
  rownames(sample_corr_mat_ordered) <- clean_sample_names[rownames(sample_corr_mat_ordered)]
  colnames(sample_corr_mat_ordered) <- clean_sample_names[colnames(sample_corr_mat_ordered)]
  
  group_counts <- table(sample_order_df$condition)
  gaps <- cumsum(as.integer(group_counts))
  gaps <- gaps[-length(gaps)]
  
  write.csv(
    sample_dist_mat_ordered,
    file.path(get_qc_results_dir(config), "sample_distance_matrix_vst.csv")
  )
  
  write.csv(
    sample_corr_mat_ordered,
    file.path(get_qc_results_dir(config), "sample_correlation_matrix_vst.csv")
  )
  
  # Grouped distance heatmap with numbers, so it is easier to see if the color scale looks uniform.
  png(
    file.path(get_qc_fig_dir(config), "sample_distance_heatmap_grouped.png"),
    width = 2000,
    height = 1800,
    res = 250
  )
  
  pheatmap(
    sample_dist_mat_ordered,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_row = gaps,
    gaps_col = gaps,
    display_numbers = round(sample_dist_mat_ordered, 1),
    number_color = "black",
    fontsize_number = 7,
    main = paste0("Sample distances using VST counts: ", config$analysis_name)
  )
  
  dev.off()
  
  # Clustered distance heatmap for QC.
  png(
    file.path(get_qc_fig_dir(config), "sample_distance_heatmap_clustered.png"),
    width = 2000,
    height = 1800,
    res = 250
  )
  
  pheatmap(
    sample_dist_mat_ordered,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = round(sample_dist_mat_ordered, 1),
    number_color = "black",
    fontsize_number = 7,
    main = paste0("Sample distances using VST counts, clustered: ", config$analysis_name)
  )
  
  dev.off()
  
  # Correlation heatmap is often easier to interpret when distance heatmap is dominated by group separation.
  png(
    file.path(get_qc_fig_dir(config), "sample_correlation_heatmap_grouped.png"),
    width = 2000,
    height = 1800,
    res = 250
  )
  
  pheatmap(
    sample_corr_mat_ordered,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_row = gaps,
    gaps_col = gaps,
    display_numbers = round(sample_corr_mat_ordered, 3),
    number_color = "black",
    fontsize_number = 7,
    main = paste0("Sample correlations using VST counts: ", config$analysis_name)
  )
  
  dev.off()
}

############################################################
# 5. DESeq2 helper functions
############################################################

export_deseq_result <- function(dds, contrast_vector, output_prefix, gene_annot, config) {
  
  res <- results(
    dds,
    contrast = contrast_vector,
    alpha = 0.05
  )
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(gene_annot, by = "ensembl_gene_id") %>%
    arrange(padj)
  
  write.csv(
    res_df,
    file.path(
      get_deseq_results_dir(config),
      paste0(output_prefix, "_DESeq2_results.csv")
    ),
    row.names = FALSE
  )
  
  return(res_df)
}

export_shrunk_lfc <- function(dds, coef_name, output_prefix, gene_annot, config) {
  
  if (!coef_name %in% resultsNames(dds)) {
    stop(
      "Coefficient '", coef_name, "' not found in resultsNames(dds). Available coefficients are: ",
      paste(resultsNames(dds), collapse = ", ")
    )
  }
  
  res_shrunk <- lfcShrink(
    dds,
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
      get_deseq_results_dir(config),
      paste0(output_prefix, "_DESeq2_shrunkLFC.csv")
    ),
    row.names = FALSE
  )
  
  return(res_shrunk_df)
}

summarize_de <- function(res_df, output_prefix, config, alpha = 0.05, lfc_cutoff = 1) {
  
  summary_df <- tibble(
    contrast = output_prefix,
    significant_padj = sum(res_df$padj < alpha, na.rm = TRUE),
    up_padj = sum(res_df$padj < alpha & res_df$log2FoldChange > 0, na.rm = TRUE),
    down_padj = sum(res_df$padj < alpha & res_df$log2FoldChange < 0, na.rm = TRUE),
    up_padj_lfc1 = sum(res_df$padj < alpha & res_df$log2FoldChange >= lfc_cutoff, na.rm = TRUE),
    down_padj_lfc1 = sum(res_df$padj < alpha & res_df$log2FoldChange <= -lfc_cutoff, na.rm = TRUE)
  )
  
  write.csv(
    summary_df,
    file.path(
      get_deseq_results_dir(config),
      paste0(output_prefix, "_DE_summary.csv")
    ),
    row.names = FALSE
  )
  
  return(summary_df)
}

plot_volcano <- function(res_df, output_prefix, config, title) {
  
  png(
    file.path(
      get_volcano_fig_dir(config),
      paste0("volcano_", output_prefix, ".png")
    ),
    width = 2000,
    height = 1800,
    res = 250
  )
  
  print(
    EnhancedVolcano(
      res_df,
      lab = clean_gene_labels(res_df),
      x = "log2FoldChange",
      y = "padj",
      title = title,
      subtitle = "DESeq2 with apeglm-shrunk log2FC",
      pCutoff = 0.05,
      FCcutoff = 1
    )
  )
  
  dev.off()
}

############################################################
# 6. GSEA helper functions
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

get_msig_gene_id_col <- function(msig_df) {
  
  if ("ncbi_gene" %in% colnames(msig_df)) {
    return("ncbi_gene")
  }
  
  if ("entrez_gene" %in% colnames(msig_df)) {
    return("entrez_gene")
  }
  
  stop(
    "Could not find an Entrez/NCBI gene ID column in msigdbr output. Available columns: ",
    paste(colnames(msig_df), collapse = ", ")
  )
}

run_hallmark_gsea <- function(gene_list, output_prefix, config) {
  
  hallmark_sets <- msigdbr(
    db_species = config$msig_db_species,
    species = config$msig_species,
    collection = config$hallmark_collection
  )
  
  gene_id_col <- get_msig_gene_id_col(hallmark_sets)
  
  hallmark_sets <- hallmark_sets %>%
    transmute(
      gs_name = gs_name,
      gene_id = as.character(.data[[gene_id_col]])
    ) %>%
    filter(!is.na(gene_id), gene_id != "")
  
  pathways <- split(
    hallmark_sets$gene_id,
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
      get_gsea_results_dir(config),
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
      get_gsea_results_dir(config),
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
      get_gsea_results_dir(config),
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
      get_gsea_results_dir(config),
      paste0(output_prefix, "_GSEA_GO_BP.csv")
    ),
    row.names = FALSE
  )
  
  return(go_res)
}

############################################################
# 7. GSEA plotting functions
############################################################

plot_fgsea_dotplot <- function(fgsea_df, output_prefix, config, title, top_n = 20) {
  
  plot_df <- fgsea_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    slice_head(n = top_n) %>%
    mutate(
      pathway = factor(pathway, levels = rev(pathway)),
      direction = ifelse(NES > 0, "Up", "Down")
    )
  
  if (nrow(plot_df) == 0) {
    message("No Hallmark GSEA rows to plot for ", output_prefix)
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
      get_gsea_fig_dir(config),
      paste0(output_prefix, "_Hallmark_dotplot.png")
    ),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

plot_top_fgsea_enrichment <- function(gsea_obj, gene_list, output_prefix, config) {
  
  fgsea_df <- gsea_obj$export
  pathways <- gsea_obj$pathways
  
  top_up <- fgsea_df %>%
    filter(!is.na(padj), !is.na(NES), NES > 0) %>%
    arrange(padj) %>%
    slice_head(n = 1) %>%
    pull(pathway)
  
  top_down <- fgsea_df %>%
    filter(!is.na(padj), !is.na(NES), NES < 0) %>%
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
        get_gsea_fig_dir(config),
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
        get_gsea_fig_dir(config),
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
    config,
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
  
  gsea_df <- gsea_df %>%
    filter(
      !is.na(ID),
      !is.na(p.adjust),
      !is.na(NES),
      is.finite(p.adjust),
      is.finite(NES)
    )
  
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
      get_gsea_fig_dir(config),
      paste0(output_prefix, "_", database_name, "_dotplot.png")
    ),
    p,
    width = fig_width,
    height = fig_height,
    dpi = 300
  )
  
  return(p)
}

save_top_gseaplot <- function(gsea_result, output_prefix, config, database_name, title_prefix) {
  
  gsea_df <- as.data.frame(gsea_result) %>%
    filter(
      !is.na(ID),
      !is.na(p.adjust),
      is.finite(p.adjust)
    )
  
  if (nrow(gsea_df) == 0) {
    message("No GSEA results to plot for: ", output_prefix, " ", database_name)
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
      get_gsea_fig_dir(config),
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
# 8. Main pipeline function
############################################################

run_pipeline <- function(config) {
  
  message("============================================================")
  message("Running analysis: ", config$analysis_name)
  message("Input file: ", config$count_file)
  message("Output folder: ", config$root_folder)
  message("============================================================")
  
  contrast_name <- paste0(config$treatment_group, "_vs_", config$control_group)
  
  contrast_list <- list()
  contrast_list[[contrast_name]] <- c("condition", config$treatment_group, config$control_group)
  
  make_output_dirs(config)
  
  ############################################################
  # Import and preprocess
  ############################################################
  
  counts_raw <- read_count_file(config$count_file)
  glimpse(counts_raw)
  
  prepared <- prepare_count_data(counts_raw)
  gene_annot <- prepared$gene_annot
  count_matrix <- prepared$count_matrix
  
  sample_info <- make_sample_info(count_matrix, config)
  print(sample_info)
  
  clean_sample_names <- make_clean_sample_names(sample_info)
  
  write_basic_qc_tables(
    count_matrix = count_matrix,
    sample_info = sample_info,
    config = config
  )
  
  ############################################################
  # DESeq2
  ############################################################
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_info,
    design = ~ condition
  )
  
  # Keep genes with at least 10 counts in at least 3 samples.
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  
  norm_counts <- counts(dds, normalized = TRUE)
  
  write.csv(
    as.data.frame(norm_counts) %>%
      rownames_to_column("ensembl_gene_id"),
    file.path(config$root_folder, "Results", "DESeq2", "normalized_counts.csv"),
    row.names = FALSE
  )
  
  ############################################################
  # QC plots
  ############################################################
  
  vsd <- vst(dds, blind = FALSE)
  
  plot_pca_custom(
    vsd = vsd,
    sample_info = sample_info,
    clean_sample_names = clean_sample_names,
    config = config,
    ntop = 500
  )
  
  plot_sample_distance_heatmaps(
    vsd = vsd,
    sample_info = sample_info,
    clean_sample_names = clean_sample_names,
    config = config
  )
  
  ############################################################
  # DESeq2 contrasts, shrinkage, summaries, volcanoes
  ############################################################
  
  message("Available DESeq2 coefficient names:")
  print(resultsNames(dds))
  
  writeLines(
    resultsNames(dds),
    file.path(config$root_folder, "Results", "DESeq2", "DESeq2_coefficient_names.txt")
  )
  
  deseq_results <- list()
  shrunk_results <- list()
  de_summaries <- list()
  
  for (contrast_name in names(contrast_list)) {
    
    message("Running contrast: ", contrast_name)
    
    deseq_results[[contrast_name]] <- export_deseq_result(
      dds = dds,
      contrast_vector = contrast_list[[contrast_name]],
      output_prefix = contrast_name,
      gene_annot = gene_annot,
      config = config
    )
    
    coef_name <- paste0("condition_", contrast_name)
    
    shrunk_results[[contrast_name]] <- export_shrunk_lfc(
      dds = dds,
      coef_name = coef_name,
      output_prefix = contrast_name,
      gene_annot = gene_annot,
      config = config
    )
    
    de_summaries[[contrast_name]] <- summarize_de(
      res_df = deseq_results[[contrast_name]],
      output_prefix = contrast_name,
      config = config
    )
    
    plot_volcano(
      res_df = shrunk_results[[contrast_name]],
      output_prefix = contrast_name,
      config = config,
      title = paste0(config$analysis_name, ": ", contrast_name)
    )
  }
  
  ############################################################
  # Ranked lists and GSEA
  ############################################################
  
  gene_lists <- list()
  gsea_hallmark <- list()
  gsea_kegg <- list()
  gsea_reactome <- list()
  gsea_go <- list()
  
  for (contrast_name in names(deseq_results)) {
    
    message("Preparing ranked gene list for: ", contrast_name)
    gene_lists[[contrast_name]] <- make_ranked_list(deseq_results[[contrast_name]])
    
    write.csv(
      tibble(
        entrezgene_id = names(gene_lists[[contrast_name]]),
        stat = as.numeric(gene_lists[[contrast_name]])
      ),
      file.path(
        get_gsea_results_dir(config),
        paste0(contrast_name, "_ranked_gene_list.csv")
      ),
      row.names = FALSE
    )
    
    message("Running Hallmark GSEA for: ", contrast_name)
    gsea_hallmark[[contrast_name]] <- run_hallmark_gsea(
      gene_list = gene_lists[[contrast_name]],
      output_prefix = contrast_name,
      config = config
    )
    
    message("Running KEGG GSEA for: ", contrast_name)
    gsea_kegg[[contrast_name]] <- run_kegg_gsea(
      gene_list = gene_lists[[contrast_name]],
      output_prefix = contrast_name,
      config = config
    )
    
    message("Running Reactome GSEA for: ", contrast_name)
    gsea_reactome[[contrast_name]] <- run_reactome_gsea(
      gene_list = gene_lists[[contrast_name]],
      output_prefix = contrast_name,
      config = config
    )
    
    message("Running GO BP GSEA for: ", contrast_name)
    gsea_go[[contrast_name]] <- run_go_bp_gsea(
      gene_list = gene_lists[[contrast_name]],
      output_prefix = contrast_name,
      config = config
    )
  }
  
  ############################################################
  # GSEA plots
  ############################################################
  
  for (contrast_name in names(gene_lists)) {
    
    plot_fgsea_dotplot(
      fgsea_df = gsea_hallmark[[contrast_name]]$export,
      output_prefix = contrast_name,
      config = config,
      title = paste0("Hallmark GSEA: ", config$analysis_name, " ", contrast_name),
      top_n = 20
    )
    
    plot_top_fgsea_enrichment(
      gsea_obj = gsea_hallmark[[contrast_name]],
      gene_list = gene_lists[[contrast_name]],
      output_prefix = contrast_name,
      config = config
    )
    
    save_gsea_dotplot(
      gsea_result = gsea_kegg[[contrast_name]],
      output_prefix = contrast_name,
      config = config,
      database_name = "KEGG",
      title = paste0("KEGG GSEA: ", config$analysis_name, " ", contrast_name),
      show_n = 10,
      label_width = 40,
      fig_width = 11,
      fig_height = 7,
      text_size = 8
    )
    
    save_gsea_dotplot(
      gsea_result = gsea_reactome[[contrast_name]],
      output_prefix = contrast_name,
      config = config,
      database_name = "Reactome",
      title = paste0("Reactome GSEA: ", config$analysis_name, " ", contrast_name),
      show_n = 10,
      label_width = 45,
      fig_width = 13,
      fig_height = 7,
      text_size = 8
    )
    
    save_gsea_dotplot(
      gsea_result = gsea_go[[contrast_name]],
      output_prefix = contrast_name,
      config = config,
      database_name = "GO_BP",
      title = paste0("GO Biological Process GSEA: ", config$analysis_name, " ", contrast_name),
      show_n = 8,
      label_width = 35,
      fig_width = 14,
      fig_height = 8,
      text_size = 7
    )
    
    save_top_gseaplot(
      gsea_result = gsea_kegg[[contrast_name]],
      output_prefix = contrast_name,
      config = config,
      database_name = "KEGG",
      title_prefix = paste0("Top KEGG GSEA ", config$analysis_name, " ", contrast_name)
    )
    
    save_top_gseaplot(
      gsea_result = gsea_reactome[[contrast_name]],
      output_prefix = contrast_name,
      config = config,
      database_name = "Reactome",
      title_prefix = paste0("Top Reactome GSEA ", config$analysis_name, " ", contrast_name)
    )
    
    save_top_gseaplot(
      gsea_result = gsea_go[[contrast_name]],
      output_prefix = contrast_name,
      config = config,
      database_name = "GO_BP",
      title_prefix = paste0("Top GO BP GSEA ", config$analysis_name, " ", contrast_name)
    )
  }
  
  ############################################################
  # Save session info
  ############################################################
  
  writeLines(
    capture.output(sessionInfo()),
    file.path(config$root_folder, "Results", "sessionInfo.txt")
  )
  
  message("Finished analysis: ", config$analysis_name)
  
  return(list(
    config = config,
    dds = dds,
    vsd = vsd,
    sample_info = sample_info,
    deseq_results = deseq_results,
    shrunk_results = shrunk_results,
    de_summaries = de_summaries,
    gene_lists = gene_lists,
    gsea_hallmark = gsea_hallmark,
    gsea_kegg = gsea_kegg,
    gsea_reactome = gsea_reactome,
    gsea_go = gsea_go
  ))
}

############################################################
# 9. Run both analyses
############################################################

mouse_UBR5_KD_results <- run_pipeline(analysis_configs$mouse_UBR5_KD)

human_UBR5_KD_results <- run_pipeline(analysis_configs$human_UBR5_KD)
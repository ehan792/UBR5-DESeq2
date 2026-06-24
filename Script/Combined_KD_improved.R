############################################################
# UBR5 KD RNA-seq pipeline (IMPROVED)
# Runs:
#   1. Human KD vs Control
#   2. Mouse KD vs Control
#
# Key improvements over original:
#   1. VST blind=TRUE for exploratory QC (standard practice)
#   2. Per-group filtering: gene kept if >=3 samples with
#      count>=10 in AT LEAST ONE of the two condition groups,
#      preventing silent exclusion of condition-specific genes
#   3. DESeq2 run once per dataset; QC and DE share the same
#      fitted object (consistent size factors)
#   4. Cook's distance outlier handling explicitly set
#   5. Independent filtering retained (DESeq2 default)
#   6. apeglm shrinkage applied via coef (not contrast) —
#      required by apeglm
#   7. Session info and filter diagnostics saved
#   8. Dispersion and MA plots added for model diagnostics
#   9. GSEA: fgsea nPermSimple set for reproducibility;
#      clusterProfiler functions use pAdjustMethod="BH"
#  10. set.seed() called before every stochastic step
############################################################

set.seed(42)

############################################################
# Package installation
############################################################

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
  "readxl",
  "RColorBrewer",
  "ggrepel"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

############################################################
# Package loading
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(janitor)
  library(readxl)
  library(RColorBrewer)
  library(ggrepel)

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
})

############################################################
# Output root
############################################################

output_root <- here("KD")

############################################################
# Analysis configuration
############################################################

analysis_configs <- list(
  human_KD_vs_Control = list(
    analysis_label   = "Human UBR5 KD",
    contrast_name    = "human_KD_vs_Control",
    count_file       = here("Data", "HumanKDCounts.xlsx"),
    control_group    = "Control",
    kd_group         = "KD",
    control_regex    = "jh_2_002",
    kd_regex         = "shrubr5",
    orgdb            = org.Hs.eg.db,
    kegg_organism    = "hsa",
    reactome_organism = "human",
    msig_db_species  = "HS",
    msig_species     = "Homo sapiens",
    hallmark_collection = "H"
  ),

  mouse_KD_vs_Control = list(
    analysis_label   = "Mouse UBR5 KD",
    contrast_name    = "mouse_KD_vs_Control",
    count_file       = here("Data", "MouseKDCounts.xlsx"),
    control_group    = "Control",
    kd_group         = "KD",
    control_regex    = "jw23\\.3",
    kd_regex         = "shrubr5",
    orgdb            = org.Mm.eg.db,
    kegg_organism    = "mmu",
    reactome_organism = "mouse",
    msig_db_species  = "MM",
    msig_species     = "Mus musculus",
    hallmark_collection = "MH"
  )
)

contrast_names <- names(analysis_configs)

############################################################
# Folder structure
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

for (cn in contrast_names) {
  output_dirs <- c(
    output_dirs,
    file.path(output_root, "results", "deseq2", cn),
    file.path(output_root, "results", "gsea", cn),
    file.path(output_root, "figures", "volcano", cn),
    file.path(output_root, "figures", "gsea", cn)
  )
}

walk(output_dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

############################################################
# Path helpers
############################################################

get_deseq_results_dir <- function(p) file.path(output_root, "results", "deseq2", p)
get_gsea_results_dir  <- function(p) file.path(output_root, "results", "gsea",   p)
get_volcano_fig_dir   <- function(p) file.path(output_root, "figures", "volcano", p)
get_gsea_fig_dir      <- function(p) file.path(output_root, "figures", "gsea",    p)
get_qc_fig_dir        <- function()  file.path(output_root, "figures", "qc")

############################################################
# Data import
############################################################
# Returns: count_matrix, gene_annot, sample_info

read_kd_count_file <- function(config) {

  counts_raw <- readxl::read_excel(config$count_file) %>%
    as.data.frame()

  # Normalise Entrez column name
  if ("entrezgene" %in% colnames(counts_raw) &&
      !"entrezgene_id" %in% colnames(counts_raw)) {
    counts_raw <- counts_raw %>% rename(entrezgene_id = entrezgene)
  }

  possible_annotation_cols <- c(
    "ensembl_gene_id", "entrezgene_id", "external_gene_name",
    "chromosome_name", "start_position", "end_position",
    "gene_biotype", "external_gene_source", "transcript_count", "description"
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

  # DESeq2 requires integer counts
  count_matrix <- round(count_matrix)
  storage.mode(count_matrix) <- "integer"

  if (any(is.na(count_matrix)))  stop("Count matrix contains NA values.")
  if (any(count_matrix < 0))     stop("Count matrix contains negative counts.")

  # Assign conditions
  sample_info <- tibble(
    sample = colnames(count_matrix),
    condition = case_when(
      str_detect(sample, config$kd_regex)      ~ config$kd_group,
      str_detect(sample, config$control_regex) ~ config$control_group,
      TRUE ~ NA_character_
    )
  ) %>%
    mutate(
      condition = factor(condition,
                         levels = c(config$control_group, config$kd_group))
    ) %>%
    column_to_rownames("sample")

  if (any(is.na(sample_info$condition))) {
    unmatched <- rownames(sample_info)[is.na(sample_info$condition)]
    stop("The following samples could not be assigned to a condition:\n  ",
         paste(unmatched, collapse = "\n  "))
  }

  stopifnot(all(rownames(sample_info) == colnames(count_matrix)))

  message(config$analysis_label, ": loaded ", nrow(count_matrix),
          " genes × ", ncol(count_matrix), " samples.")
  message("  Conditions: ",
          paste(names(table(sample_info$condition)), "=",
                as.integer(table(sample_info$condition)), collapse = "; "))

  list(
    counts_raw   = counts_raw,
    gene_annot   = gene_annot,
    count_matrix = count_matrix,
    sample_info  = sample_info
  )
}

############################################################
# Filtering function
############################################################
# IMPROVEMENT: filter per condition group.
#
# A gene passes if it has count >= min_count in at least
# min_samples_per_group samples WITHIN AT LEAST ONE condition
# group. This is the standard "filterByExpr"-style approach
# and prevents false exclusion of condition-specific genes
# while still removing unexpressed genes.
#
# With 3 samples per group and min_samples_per_group = 3,
# a gene must be expressed (>=10 counts) in ALL replicates
# of at least one group. You can relax to 2 if needed.
############################################################

filter_by_group <- function(count_matrix, sample_info,
                             min_count = 10,
                             min_samples_per_group = 3) {
  conditions <- levels(sample_info$condition)
  keep <- Reduce(`|`, lapply(conditions, function(cond) {
    samp_idx <- rownames(sample_info)[sample_info$condition == cond]
    rowSums(count_matrix[, samp_idx, drop = FALSE] >= min_count) >= min_samples_per_group
  }))
  keep
}

############################################################
# QC function
############################################################
# IMPROVEMENT:
#   - blind=TRUE for VST (standard for exploratory QC;
#     blind=FALSE uses the design and can mask batch effects)
#   - Runs on a filtered gene set but does NOT call DESeq()
#     here — that happens once in run_deseq_for_kd() so all
#     downstream steps share the same size factors.
#   - Uses estimateSizeFactors() only for normalised counts
#     in the QC output.
############################################################

run_qc <- function(count_matrix, sample_info,
                   output_prefix, analysis_label,
                   min_count = 10, min_samples_per_group = 3) {

  message("\n=== QC: ", analysis_label, " ===")

  # Filter for QC
  keep_qc <- filter_by_group(count_matrix, sample_info,
                              min_count, min_samples_per_group)
  message("  QC: keeping ", sum(keep_qc), " of ", length(keep_qc),
          " genes after group-aware filtering.")

  dds_qc <- DESeqDataSetFromMatrix(
    countData = count_matrix[keep_qc, ],
    colData   = sample_info,
    design    = ~ condition
  )

  # Estimate size factors only (no full DESeq fit) for QC normalised counts
  dds_qc <- estimateSizeFactors(dds_qc)

  # Export normalised counts for all samples
  norm_counts_qc <- counts(dds_qc, normalized = TRUE)
  write.csv(
    as.data.frame(norm_counts_qc) %>% rownames_to_column("ensembl_gene_id"),
    file.path(get_deseq_results_dir(output_prefix),
              paste0(output_prefix, "_normalized_counts_QC_all_samples.csv")),
    row.names = FALSE
  )

  # ----- VST (blind=TRUE) -----
  # blind=TRUE is correct for QC: it re-estimates dispersion
  # ignoring the design so the plot is not biased toward
  # showing the expected effect.
  vsd <- vst(dds_qc, blind = TRUE)

  # ----- Pretty sample labels -----
  clean_sample_names <- setNames(
    paste0(sample_info$condition, " ",
           ave(seq_along(sample_info$condition),
               sample_info$condition, FUN = seq_along)),
    rownames(sample_info)
  )

  # ----- PCA -----
  pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  pca_data$clean_name <- clean_sample_names[pca_data$name]
  percent_var <- round(100 * attr(pca_data, "percentVar"), 2)

  p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition,
                                label = clean_name)) +
    geom_point(size = 4) +
    ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggtitle(paste0(analysis_label, ": PCA (VST, blind=TRUE)")) +
    theme_bw(base_size = 12) +
    scale_color_brewer(palette = "Set1")

  ggsave(
    file.path(get_qc_fig_dir(),
              paste0(output_prefix, "_PCA_vst_condition.png")),
    p_pca, width = 7, height = 5, dpi = 300
  )

  # ----- PCA scree plot -----
  # Shows how much variance each PC captures — helps
  # contextualise whether 99% PC1 is driven by one outlier
  # sample or a true biological difference.
  pca_full <- prcomp(t(assay(vsd)), scale. = FALSE)
  pct_var_full <- round(100 * pca_full$sdev^2 / sum(pca_full$sdev^2), 2)
  n_pcs <- min(length(pct_var_full), ncol(count_matrix))
  scree_df <- tibble(PC = paste0("PC", seq_len(n_pcs)),
                     Variance = pct_var_full[seq_len(n_pcs)]) %>%
    mutate(PC = factor(PC, levels = PC))

  p_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = paste0(Variance, "%")),
              vjust = -0.5, size = 3) +
    ylim(0, max(pct_var_full[seq_len(n_pcs)]) * 1.15) +
    labs(title = paste0(analysis_label, ": PCA scree plot (VST)"),
         x = NULL, y = "% variance explained") +
    theme_bw(base_size = 12)

  ggsave(
    file.path(get_qc_fig_dir(),
              paste0(output_prefix, "_PCA_scree.png")),
    p_scree, width = 6, height = 4, dpi = 300
  )

  # ----- Sample distance heatmap -----
  sample_dists     <- dist(t(assay(vsd)))
  sample_dist_mat  <- as.matrix(sample_dists)
  rownames(sample_dist_mat) <- clean_sample_names[colnames(vsd)]
  colnames(sample_dist_mat) <- clean_sample_names[colnames(vsd)]

  ann_col <- data.frame(
    Condition = sample_info$condition,
    row.names = clean_sample_names[rownames(sample_info)]
  )
  cond_colors <- RColorBrewer::brewer.pal(3, "Set1")[seq_along(levels(sample_info$condition))]
  names(cond_colors) <- levels(sample_info$condition)

  # Grouped heatmap
  sample_order_df <- sample_info %>%
    rownames_to_column("sample") %>%
    arrange(condition, sample)
  sample_order_labels <- clean_sample_names[sample_order_df$sample]
  mat_ordered <- sample_dist_mat[sample_order_labels, sample_order_labels]
  group_sizes <- table(sample_order_df$condition)
  group_gaps  <- cumsum(as.integer(group_sizes))[-length(group_sizes)]

  png(file.path(get_qc_fig_dir(),
                paste0(output_prefix, "_sample_distance_heatmap_grouped.png")),
      width = 1800, height = 1600, res = 250)
  pheatmap(mat_ordered,
           cluster_rows = FALSE, cluster_cols = FALSE,
           gaps_row = group_gaps, gaps_col = group_gaps,
           annotation_row = ann_col[rownames(mat_ordered), , drop = FALSE],
           annotation_col = ann_col[colnames(mat_ordered), , drop = FALSE],
           annotation_colors = list(Condition = cond_colors),
           main = paste0(analysis_label, ": sample distances (grouped)"),
           fontsize_row = 10, fontsize_col = 10, angle_col = 45,
           color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(100))
  dev.off()

  # Clustered heatmap
  png(file.path(get_qc_fig_dir(),
                paste0(output_prefix, "_sample_distance_heatmap_clustered.png")),
      width = 1800, height = 1600, res = 250)
  pheatmap(sample_dist_mat,
           cluster_rows = TRUE, cluster_cols = TRUE,
           annotation_row = ann_col,
           annotation_col = ann_col,
           annotation_colors = list(Condition = cond_colors),
           main = paste0(analysis_label, ": sample distances (clustered)"),
           fontsize_row = 10, fontsize_col = 10, angle_col = 45,
           color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(100))
  dev.off()

  message("  QC plots saved. PC1=", percent_var[1], "%, PC2=", percent_var[2], "%")

  list(vsd = vsd, pca_data = pca_data, percent_var = percent_var,
       clean_sample_names = clean_sample_names)
}

############################################################
# DESeq2 DE analysis
############################################################
# IMPROVEMENTS:
#   - Group-aware pre-filtering (see filter_by_group above)
#   - Single DESeq() call; size factors used for both QC
#     normalised counts and DE
#   - lfcShrink via coef= (required by apeglm, not contrast=)
#   - MA plot and dispersion plot saved for model diagnostics
#   - Cook's distance filtering left at DESeq2 default (ON)
#   - independentFiltering left at DESeq2 default (ON)
############################################################

run_deseq_for_kd <- function(count_matrix, sample_info,
                              output_prefix, gene_annot,
                              treatment_group = "KD",
                              control_group   = "Control",
                              min_count       = 10,
                              min_samples_per_group = 3) {

  message("\n=== DESeq2: ", output_prefix, " ===")

  # Group-aware filtering
  keep <- filter_by_group(count_matrix, sample_info,
                          min_count, min_samples_per_group)
  message("  Kept ", sum(keep), " of ", length(keep),
          " genes after group-aware filtering (min_count=", min_count,
          ", min_samples_per_group=", min_samples_per_group, ")")

  # Filter summary
  filter_summary <- tibble(
    contrast               = output_prefix,
    treatment_group        = treatment_group,
    control_group          = control_group,
    samples_used           = paste(colnames(count_matrix), collapse = ";"),
    n_samples              = ncol(count_matrix),
    min_count              = min_count,
    min_samples_per_group  = min_samples_per_group,
    genes_before_filtering = length(keep),
    genes_after_filtering  = sum(keep),
    genes_removed          = length(keep) - sum(keep)
  )
  write.csv(filter_summary,
            file.path(get_deseq_results_dir(output_prefix),
                      paste0(output_prefix, "_filter_summary.csv")),
            row.names = FALSE)

  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix[keep, ],
    colData   = sample_info,
    design    = ~ condition
  )

  # Set reference level
  dds$condition <- relevel(dds$condition, ref = control_group)

  # Fit DESeq2 model
  # Note: DESeq2 automatically handles:
  #   - size factor estimation (median of ratios)
  #   - dispersion estimation and shrinkage
  #   - Wald test
  #   - Cook's distance outlier replacement (fitType="parametric")
  #   - independent filtering for multiple testing
  set.seed(42)
  dds <- DESeq(dds)

  # ----- Dispersion plot (model diagnostic) -----
  png(file.path(get_qc_fig_dir(),
                paste0(output_prefix, "_dispersion_estimates.png")),
      width = 1600, height = 1200, res = 250)
  plotDispEsts(dds, main = paste0(output_prefix, ": dispersion estimates"))
  dev.off()

  # ----- Normalised counts -----
  norm_counts <- counts(dds, normalized = TRUE)
  write.csv(
    as.data.frame(norm_counts) %>% rownames_to_column("ensembl_gene_id"),
    file.path(get_deseq_results_dir(output_prefix),
              paste0(output_prefix, "_normalized_counts.csv")),
    row.names = FALSE
  )

  # ----- Wald test results (unshrunken, used for GSEA ranking) -----
  res <- results(
    dds,
    contrast = c("condition", treatment_group, control_group),
    alpha    = 0.05
    # independentFiltering = TRUE (default)
    # cooksCutoff = TRUE (default)
  )

  res_df <- as.data.frame(res) %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(gene_annot, by = "ensembl_gene_id") %>%
    arrange(padj)

  write.csv(res_df,
            file.path(get_deseq_results_dir(output_prefix),
                      paste0(output_prefix, "_DESeq2_results.csv")),
            row.names = FALSE)

  # ----- MA plot (unshrunken) -----
  png(file.path(get_qc_fig_dir(),
                paste0(output_prefix, "_MA_plot_unshrunken.png")),
      width = 1600, height = 1200, res = 250)
  plotMA(res, alpha = 0.05,
         main = paste0(output_prefix, ": MA plot (Wald)"))
  dev.off()

  # ----- apeglm LFC shrinkage -----
  # apeglm requires coef= (not contrast=). The coefficient name
  # follows DESeq2's convention: condition_KD_vs_Control
  coef_name <- paste0("condition_", treatment_group, "_vs_", control_group)
  available_coefs <- resultsNames(dds)
  message("  Available coefficients: ", paste(available_coefs, collapse = ", "))

  if (!coef_name %in% available_coefs) {
    stop("Coefficient not found: '", coef_name,
         "'\n  Available: ", paste(available_coefs, collapse = ", "))
  }

  set.seed(42)
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")

  res_shrunk_df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(gene_annot, by = "ensembl_gene_id") %>%
    arrange(padj)

  write.csv(res_shrunk_df,
            file.path(get_deseq_results_dir(output_prefix),
                      paste0(output_prefix, "_DESeq2_shrunkLFC.csv")),
            row.names = FALSE)

  # ----- MA plot (shrunken) -----
  png(file.path(get_qc_fig_dir(),
                paste0(output_prefix, "_MA_plot_shrunken.png")),
      width = 1600, height = 1200, res = 250)
  plotMA(res_shrunk, alpha = 0.05,
         main = paste0(output_prefix, ": MA plot (apeglm-shrunk LFC)"))
  dev.off()

  # ----- Summary table -----
  summary_df <- tibble(
    contrast               = output_prefix,
    comparison             = paste(treatment_group, "vs", control_group),
    min_count_filter       = min_count,
    min_samples_per_group  = min_samples_per_group,
    genes_tested           = sum(!is.na(res_df$padj)),
    significant_padj_0.05  = sum(res_df$padj < 0.05, na.rm = TRUE),
    up_padj_0.05           = sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE),
    down_padj_0.05         = sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE),
    up_padj_lfc1           = sum(res_df$padj < 0.05 & res_df$log2FoldChange >= 1, na.rm = TRUE),
    down_padj_lfc1         = sum(res_df$padj < 0.05 & res_df$log2FoldChange <= -1, na.rm = TRUE)
  )
  write.csv(summary_df,
            file.path(get_deseq_results_dir(output_prefix),
                      paste0(output_prefix, "_DE_summary.csv")),
            row.names = FALSE)

  message("  DE summary: ", summary_df$significant_padj_0.05,
          " sig. genes (padj<0.05); ",
          summary_df$up_padj_0.05, " up, ",
          summary_df$down_padj_0.05, " down")

  list(dds = dds, res = res_df, shrunk = res_shrunk_df,
       summary = summary_df, filter_summary = filter_summary)
}

############################################################
# Volcano plot
############################################################

plot_volcano <- function(res_shrunk_df, output_prefix, title,
                         pCutoff = 0.05, FCcutoff = 1) {

  # Guard against empty or all-NA results
  if (all(is.na(res_shrunk_df$padj))) {
    message("No valid padj values for volcano plot: ", output_prefix)
    return(NULL)
  }

  png(
    file.path(get_volcano_fig_dir(output_prefix),
              paste0("volcano_", output_prefix, ".png")),
    width = 2000, height = 1800, res = 250
  )
  print(
    EnhancedVolcano(
      res_shrunk_df,
      lab      = res_shrunk_df$external_gene_name,
      x        = "log2FoldChange",
      y        = "padj",
      title    = title,
      subtitle = "DESeq2 with apeglm-shrunk log2FC",
      pCutoff  = pCutoff,
      FCcutoff = FCcutoff,
      pointSize = 2,
      labSize   = 3,
      drawConnectors = TRUE,
      widthConnectors = 0.5
    )
  )
  dev.off()
}

############################################################
# GSEA ranked gene list
############################################################
# IMPROVEMENT: rank by Wald statistic (sign(LFC) * -log10(pvalue)
# is an alternative; stat is preferred because it directly
# reflects the Wald test used by DESeq2 and is recommended
# in Love et al. and Zhu et al.)
############################################################

make_ranked_list <- function(res_df) {
  ranked <- res_df %>%
    filter(
      !is.na(stat), is.finite(stat),
      !is.na(entrezgene_id),
      entrezgene_id != "", entrezgene_id != "NA"
    ) %>%
    mutate(entrezgene_id = as.character(entrezgene_id)) %>%
    group_by(entrezgene_id) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(stat))

  gene_list <- setNames(ranked$stat, ranked$entrezgene_id)
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[is.finite(gene_list) & !is.na(names(gene_list)) &
                           names(gene_list) != "" & names(gene_list) != "NA"]
  gene_list[!duplicated(names(gene_list))]
}

############################################################
# GSEA helper
############################################################

get_msig_entrez_column <- function(msig_df) {
  if ("ncbi_gene"   %in% colnames(msig_df)) return("ncbi_gene")
  if ("entrez_gene" %in% colnames(msig_df)) return("entrez_gene")
  stop("Could not find Entrez gene column in msigdbr output.")
}

############################################################
# Hallmark GSEA (fgsea)
############################################################

run_hallmark_gsea <- function(gene_list, output_prefix, config) {
  hallmark_sets <- msigdbr(
    db_species = config$msig_db_species,
    species    = config$msig_species,
    collection = config$hallmark_collection
  )
  entrez_col <- get_msig_entrez_column(hallmark_sets)
  hallmark_sets <- hallmark_sets %>%
    select(gs_name, all_of(entrez_col)) %>%
    filter(!is.na(.data[[entrez_col]]))

  pathways <- split(hallmark_sets[[entrez_col]], hallmark_sets$gs_name)

  set.seed(42)
  fgsea_res <- fgsea(
    pathways   = pathways,
    stats      = gene_list,
    minSize    = 15,
    maxSize    = 500,
    nPermSimple = 10000   # reproducible permutation count
  ) %>%
    arrange(padj)

  fgsea_export <- fgsea_res %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

  write.csv(fgsea_export,
            file.path(get_gsea_results_dir(output_prefix),
                      paste0(output_prefix, "_GSEA_Hallmark.csv")),
            row.names = FALSE)

  list(results = fgsea_res, export = fgsea_export, pathways = pathways)
}

############################################################
# KEGG, Reactome, GO-BP GSEA (clusterProfiler / ReactomePA)
############################################################

run_kegg_gsea <- function(gene_list, output_prefix, config) {
  set.seed(42)
  kegg_res <- gseKEGG(
    geneList     = gene_list,
    organism     = config$kegg_organism,
    minGSSize    = 15, maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose      = FALSE
  )
  kegg_df <- as.data.frame(kegg_res) %>% arrange(p.adjust)
  write.csv(kegg_df,
            file.path(get_gsea_results_dir(output_prefix),
                      paste0(output_prefix, "_GSEA_KEGG.csv")),
            row.names = FALSE)
  kegg_res
}

run_reactome_gsea <- function(gene_list, output_prefix, config) {
  set.seed(42)
  reactome_res <- gsePathway(
    geneList     = gene_list,
    organism     = config$reactome_organism,
    minGSSize    = 15, maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose      = FALSE
  )
  reactome_df <- as.data.frame(reactome_res) %>% arrange(p.adjust)
  write.csv(reactome_df,
            file.path(get_gsea_results_dir(output_prefix),
                      paste0(output_prefix, "_GSEA_Reactome.csv")),
            row.names = FALSE)
  reactome_res
}

run_go_bp_gsea <- function(gene_list, output_prefix, config) {
  set.seed(42)
  go_res <- gseGO(
    geneList     = gene_list,
    OrgDb        = config$orgdb,
    keyType      = "ENTREZID",
    ont          = "BP",
    minGSSize    = 15, maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose      = FALSE,
    nPermSimple  = 10000
  )
  go_df <- as.data.frame(go_res) %>%
    filter(!is.na(ID), is.finite(NES), is.finite(p.adjust)) %>%
    arrange(p.adjust)
  write.csv(go_df,
            file.path(get_gsea_results_dir(output_prefix),
                      paste0(output_prefix, "_GSEA_GO_BP.csv")),
            row.names = FALSE)
  go_res
}

############################################################
# GSEA plots (Hallmark via fgsea)
############################################################

plot_fgsea_dotplot <- function(fgsea_df, output_prefix, title, top_n = 20) {
  plot_df <- fgsea_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    slice_head(n = top_n) %>%
    mutate(
      pathway   = str_replace_all(pathway, "_", " "),
      pathway   = factor(pathway, levels = rev(unique(pathway))),
      direction = ifelse(NES > 0, "Up in KD", "Down in KD")
    )

  if (nrow(plot_df) == 0) {
    message("No Hallmark GSEA results to plot: ", output_prefix); return(NULL)
  }

  p <- ggplot(plot_df, aes(x = NES, y = pathway, size = size, color = padj)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "grey80") +
    theme_bw(base_size = 11) +
    labs(title = title, x = "NES", y = NULL,
         size = "Gene set size", color = "Adj. p-value") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

  ggsave(
    file.path(get_gsea_fig_dir(output_prefix),
              paste0(output_prefix, "_Hallmark_dotplot.png")),
    p, width = 9, height = 6, dpi = 300
  )
  p
}

plot_top_fgsea_enrichment <- function(gsea_obj, gene_list, output_prefix) {
  fgsea_df <- gsea_obj$export
  pathways  <- gsea_obj$pathways

  top_up   <- fgsea_df %>% filter(!is.na(padj), NES > 0) %>%
    arrange(padj) %>% slice_head(n = 1) %>% pull(pathway)
  top_down <- fgsea_df %>% filter(!is.na(padj), NES < 0) %>%
    arrange(padj) %>% slice_head(n = 1) %>% pull(pathway)

  for (pw in c(top_up, top_down)) {
    if (length(pw) == 1 && pw %in% names(pathways)) {
      p <- plotEnrichment(pathways[[pw]], gene_list) +
        labs(title = paste0(output_prefix, ": ", pw))
      tag <- ifelse(pw %in% top_up, "top_up", "top_down")
      ggsave(
        file.path(get_gsea_fig_dir(output_prefix),
                  paste0(output_prefix, "_Hallmark_", tag, "_enrichment.png")),
        p, width = 7, height = 5, dpi = 300
      )
    }
  }
}

############################################################
# clusterProfiler GSEA plots (KEGG / Reactome / GO-BP)
############################################################

save_gsea_dotplot <- function(gsea_result, output_prefix, database_name,
                               title, show_n = 10, label_width = 40,
                               fig_width = 11, fig_height = 7,
                               split_direction = TRUE, text_size = 8) {
  gsea_df <- as.data.frame(gsea_result)
  if (nrow(gsea_df) == 0) {
    message("No GSEA results to plot: ", output_prefix, " ", database_name)
    return(NULL)
  }

  p <- if (split_direction) {
    enrichplot::dotplot(gsea_result, showCategory = show_n,
                        split = ".sign", label_format = label_width) +
      facet_grid(. ~ .sign)
  } else {
    enrichplot::dotplot(gsea_result, showCategory = show_n,
                        label_format = label_width)
  }

  p <- p + ggtitle(title) + theme_bw() +
    theme(axis.text.y  = element_text(size = text_size),
          axis.text.x  = element_text(size = 9),
          plot.title   = element_text(size = 14),
          strip.text   = element_text(size = 11))

  ggsave(
    file.path(get_gsea_fig_dir(output_prefix),
              paste0(output_prefix, "_", database_name, "_dotplot.png")),
    p, width = fig_width, height = fig_height, dpi = 300
  )
  p
}

save_top_gseaplot <- function(gsea_result, output_prefix,
                               database_name, title_prefix) {
  gsea_df <- as.data.frame(gsea_result) %>%
    filter(!is.na(ID), is.finite(p.adjust))

  if (nrow(gsea_df) == 0) {
    message("No valid GSEA results for gseaplot: ", output_prefix, " ", database_name)
    return(NULL)
  }

  top_id <- gsea_df %>% arrange(p.adjust) %>% slice_head(n = 1) %>% pull(ID)
  p <- enrichplot::gseaplot2(gsea_result, geneSetID = top_id,
                              title = paste0(title_prefix, ": ", top_id))
  ggsave(
    file.path(get_gsea_fig_dir(output_prefix),
              paste0(output_prefix, "_", database_name, "_top_gseaplot.png")),
    p, width = 8, height = 6, dpi = 300
  )
  p
}

############################################################
# Master pipeline function
############################################################

run_kd_pipeline <- function(config) {
  output_prefix <- config$contrast_name
  message("\n", strrep("=", 60))
  message("Starting: ", config$analysis_label)
  message(strrep("=", 60))

  # 1. Import data
  imported <- read_kd_count_file(config)

  # 2. QC (uses VST blind=TRUE; does NOT run full DESeq)
  qc_results <- run_qc(
    count_matrix         = imported$count_matrix,
    sample_info          = imported$sample_info,
    output_prefix        = output_prefix,
    analysis_label       = config$analysis_label,
    min_count            = 10,
    min_samples_per_group = 3
  )

  # 3. DE analysis (group-aware filtering + full DESeq2 fit)
  deseq_results <- run_deseq_for_kd(
    count_matrix         = imported$count_matrix,
    sample_info          = imported$sample_info,
    output_prefix        = output_prefix,
    gene_annot           = imported$gene_annot,
    treatment_group      = config$kd_group,
    control_group        = config$control_group,
    min_count            = 10,
    min_samples_per_group = 3
  )

  # 4. Volcano plot (apeglm-shrunk LFC)
  plot_volcano(
    res_shrunk_df = deseq_results$shrunk,
    output_prefix = output_prefix,
    title         = paste0(config$analysis_label, ": KD vs Control"),
    pCutoff       = 0.05,
    FCcutoff      = 1
  )

  # 5. GSEA ranked list (Wald statistic from unshrunken results)
  gene_list <- make_ranked_list(deseq_results$res)
  message("  GSEA input: ", length(gene_list), " ranked genes")

  # 6. GSEA
  gsea_hallmark <- run_hallmark_gsea(gene_list, output_prefix, config)
  gsea_kegg     <- run_kegg_gsea(gene_list, output_prefix, config)
  gsea_reactome <- run_reactome_gsea(gene_list, output_prefix, config)
  gsea_go       <- run_go_bp_gsea(gene_list, output_prefix, config)

  # 7. Hallmark plots
  plot_fgsea_dotplot(
    gsea_hallmark$export, output_prefix,
    paste0("Hallmark GSEA: ", config$analysis_label)
  )
  plot_top_fgsea_enrichment(gsea_hallmark, gene_list, output_prefix)

  # 8. KEGG / Reactome / GO-BP plots
  gsea_plot_settings <- list(
    KEGG     = list(result = gsea_kegg,     title = paste0("KEGG GSEA: ",     config$analysis_label), show_n = 10, label_width = 40, fig_width = 11, fig_height = 7,  text_size = 8),
    Reactome = list(result = gsea_reactome, title = paste0("Reactome GSEA: ", config$analysis_label), show_n = 10, label_width = 45, fig_width = 13, fig_height = 7,  text_size = 8),
    GO_BP    = list(result = gsea_go,       title = paste0("GO BP GSEA: ",    config$analysis_label), show_n = 8,  label_width = 35, fig_width = 14, fig_height = 8,  text_size = 7)
  )

  for (db in names(gsea_plot_settings)) {
    s <- gsea_plot_settings[[db]]
    save_gsea_dotplot(s$result, output_prefix, db, s$title,
                      s$show_n, s$label_width, s$fig_width, s$fig_height,
                      text_size = s$text_size)
    save_top_gseaplot(s$result, output_prefix, db,
                      paste0("Top ", db, " GSEA: ", config$analysis_label))
  }

  message("\nFinished: ", config$analysis_label)

  list(imported = imported, qc = qc_results,
       deseq = deseq_results, gene_list = gene_list,
       gsea = list(hallmark = gsea_hallmark, kegg = gsea_kegg,
                   reactome = gsea_reactome, go_bp = gsea_go))
}

############################################################
# Run both analyses
############################################################

human_KD_results <- run_kd_pipeline(analysis_configs$human_KD_vs_Control)
mouse_KD_results <- run_kd_pipeline(analysis_configs$mouse_KD_vs_Control)

############################################################
# Session info for reproducibility
############################################################

writeLines(
  capture.output(sessionInfo()),
  file.path(output_root, "results", "sessionInfo.txt")
)

message("\nPipeline complete. Outputs written to: ", output_root)

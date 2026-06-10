#mouse CRISPR KO mouse

#############################################
#Package download
#############################################
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "DESeq2",
  "apeglm",
  "org.Mm.eg.db",
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
  "janitor"
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
#Analysis & Visualization Tool Importing
############################################################
library(tidyverse)
library(here)
library(janitor)

library(DESeq2)
library(apeglm)
library(org.Mm.eg.db)
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
# All CRISPR KO/Het analysis outputs will be written inside:
#   CRISPR/results/
#   CRISPR/figures/
#
# The Data folder remains in the main project folder:
#   Data/MouseKOCounts.csv
############################################################

output_root <- here("CRISPR")

############################################################
# Folder organization
############################################################

output_dirs <- c(
  file.path(output_root, "results"),
  file.path(output_root, "results", "deseq2"),
  file.path(output_root, "results", "deseq2", "het_vs_WT"),
  file.path(output_root, "results", "deseq2", "KO_vs_WT"),
  file.path(output_root, "results", "gsea"),
  file.path(output_root, "results", "gsea", "het_vs_WT"),
  file.path(output_root, "results", "gsea", "KO_vs_WT"),
  
  file.path(output_root, "figures"),
  file.path(output_root, "figures", "qc"),
  file.path(output_root, "figures", "volcano"),
  file.path(output_root, "figures", "volcano", "het_vs_WT"),
  file.path(output_root, "figures", "volcano", "KO_vs_WT"),
  file.path(output_root, "figures", "gsea"),
  file.path(output_root, "figures", "gsea", "het_vs_WT"),
  file.path(output_root, "figures", "gsea", "KO_vs_WT")
)

walk(output_dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

############################################################
# Helper functions for output paths
############################################################

get_deseq_results_dir <- function(output_prefix) {
  if (output_prefix == "het_vs_WT") {
    return(file.path(output_root, "results", "deseq2", "het_vs_WT"))
  } else if (output_prefix == "KO_vs_WT") {
    return(file.path(output_root, "results", "deseq2", "KO_vs_WT"))
  } else {
    return(file.path(output_root, "results", "deseq2"))
  }
}

get_volcano_fig_dir <- function(output_prefix) {
  if (output_prefix == "het_vs_WT") {
    return(file.path(output_root, "figures", "volcano", "het_vs_WT"))
  } else if (output_prefix == "KO_vs_WT") {
    return(file.path(output_root, "figures", "volcano", "KO_vs_WT"))
  } else {
    return(file.path(output_root, "figures", "volcano"))
  }
}

get_gsea_results_dir <- function(output_prefix) {
  if (output_prefix == "het_vs_WT") {
    return(file.path(output_root, "results", "gsea", "het_vs_WT"))
  } else if (output_prefix == "KO_vs_WT") {
    return(file.path(output_root, "results", "gsea", "KO_vs_WT"))
  } else {
    return(file.path(output_root, "results", "gsea"))
  }
}

get_gsea_fig_dir <- function(output_prefix) {
  if (output_prefix == "het_vs_WT") {
    return(file.path(output_root, "figures", "gsea", "het_vs_WT"))
  } else if (output_prefix == "KO_vs_WT") {
    return(file.path(output_root, "figures", "gsea", "KO_vs_WT"))
  } else {
    return(file.path(output_root, "figures", "gsea"))
  }
}

############################################################
# Data Import
############################################################

counts_raw <- read.csv(
  here("Data", "MouseKOCounts.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

glimpse(counts_raw)

############################################################
# Pre-processing
############################################################

# Separate annotation columns and count columns
annotation_cols <- c(
  "ensembl_gene_id",
  "entrezgene_id",
  "external_gene_name",
  "chromosome_name",
  "start_position",
  "end_position",
  "gene_biotype",
  "description"
)

#value= true makes the matrix comprised of the sample data, not just index
count_cols <- grep("^sample\\.", colnames(counts_raw), value = TRUE)

gene_annot <- counts_raw %>%
  select(all_of(annotation_cols))


count_matrix <- counts_raw %>%
  select(ensembl_gene_id, all_of(count_cols)) %>%
  column_to_rownames("ensembl_gene_id") %>%
  as.matrix()

storage.mode(count_matrix) <- "integer"


# Make sample metadata
sample_info <- tibble(
  sample = colnames(count_matrix),
  genotype = case_when(
    str_detect(sample, "WT") ~ "WT",
    str_detect(sample, "pos.neg") ~ "het",
    str_detect(sample, "neg.neg") ~ "KO",
    TRUE ~ NA_character_
  )
) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "het", "KO"))
  ) %>%
  column_to_rownames("sample")

sample_info

#formatting check
stopifnot(all(rownames(sample_info) == colnames(count_matrix)))
stopifnot(!any(is.na(sample_info$genotype)))



############################################################
# 5. Create full DESeq2 object for QC only
############################################################

dds_all <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ genotype
)

############################################################
# 6. QC-filtered object for PCA / heatmaps only
############################################################
# This object is for global QC visualization.
# Differential expression contrasts will be filtered separately below.

keep_qc <- rowSums(counts(dds_all) >= 10) >= 3
dds_qc <- dds_all[keep_qc, ]

dds_qc <- DESeq(dds_qc)

# Save globally normalized counts for QC/reference
norm_counts_qc <- counts(dds_qc, normalized = TRUE)

write.csv(
  as.data.frame(norm_counts_qc) %>%
    rownames_to_column("ensembl_gene_id"),
  file.path(output_root, "results", "deseq2", "normalized_counts_QC_all_groups.csv"),
  row.names = FALSE
)




############################################################
# 8. QC: Initial visualizations (PCA, heatmap)
############################################################
clean_sample_names <- c(
  "sample.WT_1" = "WT 1",
  "sample.WT_2" = "WT 2",
  "sample.WT_3" = "WT 3",
  "sample.pos.neg_1" = "Het 1",
  "sample.pos.neg_2" = "Het 2",
  "sample.pos.neg_3" = "Het 3",
  "sample.neg.neg_1" = "KO 1",
  "sample.neg.neg_2" = "KO 2",
  "sample.neg.neg_3" = "KO 3"
)


# PCA using VST-transformed counts
vsd <- vst(dds_qc, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "genotype", returnData = TRUE)
pca_data$clean_name <- clean_sample_names[pca_data$name]

percent_var <- round(100 * attr(pca_data, "percentVar"), 2)

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = genotype, label = clean_name)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw()

ggsave(
  file.path(output_root, "figures", "qc", "PCA_vst_genotype.png"),
  p_pca,
  width = 7,
  height = 5,
  dpi = 300
)

############################################################
# 9. Sample distance heatmaps
############################################################

sample_dists <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dists)

rownames(sample_dist_mat) <- colnames(vsd)
colnames(sample_dist_mat) <- colnames(vsd)

# Order samples by genotype so the grouped heatmap is WT, het, KO
sample_order_df <- sample_info %>%
  rownames_to_column("sample") %>%
  arrange(genotype, sample)

sample_order <- sample_order_df$sample

sample_dist_mat_ordered <- sample_dist_mat[sample_order, sample_order]

rownames(sample_dist_mat_ordered) <- clean_sample_names[rownames(sample_dist_mat_ordered)]
colnames(sample_dist_mat_ordered) <- clean_sample_names[colnames(sample_dist_mat_ordered)]

# Automatically calculate group gaps instead of hardcoding c(3, 6)
group_sizes <- table(sample_order_df$genotype)
group_gaps <- cumsum(as.integer(group_sizes))
group_gaps <- group_gaps[-length(group_gaps)]

############################################################
# 9A. Grouped heatmap: preserves phenotype order
############################################################

png(
  file.path(output_root, "figures", "qc", "sample_distance_heatmap_grouped.png"),
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
  main = "Sample distances using VST counts, grouped by genotype",
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45
)

dev.off()

############################################################
# 9B. Clustered heatmap: lets samples cluster naturally
############################################################

png(
  file.path(output_root, "figures", "qc", "sample_distance_heatmap_clustered.png"),
  width = 1800,
  height = 1600,
  res = 250
)

pheatmap(
  sample_dist_mat_ordered,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Sample distances using VST counts, clustered",
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45
)

dev.off()


############################################################
# 10. Contrast-specific DESeq2 analysis function
############################################################
# This function runs a full DESeq2 workflow for ONE comparison at a time.
#
# Why this matters:
#   Instead of filtering genes globally across all 9 samples,
#   this function first subsets to only the two genotypes being compared.
#
# Example:
#   het_vs_WT uses only WT + het samples.
#   KO_vs_WT uses only WT + KO samples.
#
# Important parameters:
#   dds_all:
#     The unfiltered DESeq2 object containing all samples.
#
#   group_a:
#     The numerator group in the contrast.
#     For c("genotype", "KO", "WT"), group_a = "KO".
#     Positive log2FC means higher in group_a.
#
#   group_b:
#     The denominator / reference group.
#     For c("genotype", "KO", "WT"), group_b = "WT".
#
#   output_prefix:
#     Used for output folder routing and file names.
#     Examples: "het_vs_WT", "KO_vs_WT".
#
#   min_count:
#     Minimum raw count required for a sample to count as "expressed".
#     Here, min_count = 10.
#
#   min_samples:
#     Number of samples in the comparison that must pass min_count.
#     Here, min_samples = 3.
#     Since each contrast has 6 samples, this means:
#       keep genes with >=10 counts in at least 3 of the 6 contrast-specific samples.
############################################################

run_deseq_for_contrast <- function(
    dds_all,
    group_a,
    group_b,
    output_prefix,
    gene_annot,
    min_count = 10,
    min_samples = 3
) {
  
  ############################################################
  # 10A. Subset to only the samples in this comparison
  ############################################################
  
  samples_to_keep <- rownames(colData(dds_all))[
    dds_all$genotype %in% c(group_a, group_b)
  ]
  
  dds_contrast <- dds_all[, samples_to_keep]
  
  # Drop the unused genotype level.
  # Example: in KO_vs_WT, the "het" level is removed.
  dds_contrast$genotype <- droplevels(dds_contrast$genotype)
  
  # Set the reference group.
  # This ensures the coefficient is group_a vs group_b.
  dds_contrast$genotype <- relevel(dds_contrast$genotype, ref = group_b)
  
  ############################################################
  # 10B. Contrast-specific low-count filtering
  ############################################################
  # This is the key change from your old global filter.
  # Filtering is applied only after subsetting to the two groups being compared.
  #
  # For each gene:
  #   keep if raw count >= min_count in at least min_samples samples.
  #
  # With min_count = 10 and min_samples = 3:
  #   keep genes with >=10 counts in at least 3 of the 6 samples.
  ############################################################
  
  keep <- rowSums(counts(dds_contrast) >= min_count) >= min_samples
  dds_contrast <- dds_contrast[keep, ]
  
  message(output_prefix, ": kept ", sum(keep), " genes after filtering.")
  
  filter_summary <- tibble(
    contrast = output_prefix,
    group_a = group_a,
    group_b = group_b,
    samples_used = paste(samples_to_keep, collapse = ";"),
    number_of_samples_used = length(samples_to_keep),
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
  
  ############################################################
  # 10C. Run DESeq2 for this contrast-specific dataset
  ############################################################
  
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
  
  ############################################################
  # 10D. Extract standard DESeq2 results
  ############################################################
  # Positive log2FoldChange means higher in group_a than group_b.
  # Example:
  #   group_a = "KO", group_b = "WT"
  #   positive log2FC = higher in KO.
  ############################################################
  
  res <- results(
    dds_contrast,
    contrast = c("genotype", group_a, group_b),
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
  
  ############################################################
  # 10E. Shrink log2FC for visualization
  ############################################################
  # Shrunk LFCs are preferred for volcano plots because they reduce
  # unstable large fold-changes from low-count genes.
  #
  # Important:
  #   GSEA should still use the unshrunk Wald statistic from res_df$stat.
  #   Volcano plots can use the shrunk log2FoldChange.
  ############################################################
  
  message("Available coefficients for ", output_prefix, ":")
  print(resultsNames(dds_contrast))
  
  coef_name <- paste0("genotype_", group_a, "_vs_", group_b)
  
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
  
  ############################################################
  # 10F. Export compact DE summary
  ############################################################
  
  summary_df <- tibble(
    contrast = output_prefix,
    comparison = paste(group_a, "vs", group_b),
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
# 11. Run contrast-specific DESeq2 analyses
############################################################
# Each comparison gets its own DESeq2 object and its own filtering.
# This avoids using the het samples when filtering KO_vs_WT,
# and avoids using the KO samples when filtering het_vs_WT.
############################################################

contrast_specs <- list(
  het_vs_WT = list(
    group_a = "het",
    group_b = "WT",
    title = "+/- vs WT"
  ),
  KO_vs_WT = list(
    group_a = "KO",
    group_b = "WT",
    title = "-/- vs WT"
  )
)

contrast_results <- list()

for (contrast_name in names(contrast_specs)) {
  
  spec <- contrast_specs[[contrast_name]]
  
  contrast_results[[contrast_name]] <- run_deseq_for_contrast(
    dds_all = dds_all,
    group_a = spec$group_a,
    group_b = spec$group_b,
    output_prefix = contrast_name,
    gene_annot = gene_annot,
    min_count = 10,
    min_samples = 3
  )
}

# Preserve familiar object names for the rest of the script.
res_het_vs_wt <- contrast_results[["het_vs_WT"]]$res
res_KO_vs_wt <- contrast_results[["KO_vs_WT"]]$res

res_het_shrunk_df <- contrast_results[["het_vs_WT"]]$shrunk
res_KO_shrunk_df <- contrast_results[["KO_vs_WT"]]$shrunk

de_summary_het <- contrast_results[["het_vs_WT"]]$summary
de_summary_KO <- contrast_results[["KO_vs_WT"]]$summary

############################################################
# 12. Volcano plot function
############################################################
# This removes duplicate volcano code and prevents copy/paste errors.
#
# Important parameters:
#   res_shrunk_df:
#     DESeq2 result table using apeglm-shrunk log2FC.
#
#   output_prefix:
#     Controls the output folder and filename.
#
#   title:
#     Human-readable title shown on the plot.
#
#   pCutoff:
#     Adjusted p-value threshold used by EnhancedVolcano.
#
#   FCcutoff:
#     Absolute log2FC cutoff used by EnhancedVolcano.
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
# 13. Volcano plots
############################################################

for (contrast_name in names(contrast_specs)) {
  
  plot_volcano(
    res_shrunk_df = contrast_results[[contrast_name]]$shrunk,
    output_prefix = contrast_name,
    title = contrast_specs[[contrast_name]]$title,
    pCutoff = 0.05,
    FCcutoff = 1
  )
}


############################################################
# 15. Prepare ranked gene lists for GSEA
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

gene_list_het <- make_ranked_list(res_het_vs_wt)
gene_list_KO <- make_ranked_list(res_KO_vs_wt)


############################################################
# 16. GSEA Hallmark using mouse-native MSigDB + fgsea
############################################################

run_hallmark_gsea <- function(gene_list, output_prefix) {
  
  hallmark_sets <- msigdbr(
    db_species = "MM",
    species = "Mus musculus",
    collection = "MH"
  ) %>%
    dplyr::select(gs_name, ncbi_gene) %>%
    dplyr::filter(!is.na(ncbi_gene))
  
  pathways <- split(
    hallmark_sets$ncbi_gene,
    hallmark_sets$gs_name
  )
  
  fgsea_res <- fgsea(
    pathways = pathways,
    stats = gene_list,
    minSize = 15,
    maxSize = 500
  ) %>%
    dplyr::arrange(padj)
  
  # Exportable version: convert list-column to text
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

gsea_hallmark_het <- run_hallmark_gsea(gene_list_het, "het_vs_WT")
gsea_hallmark_KO <- run_hallmark_gsea(gene_list_KO, "KO_vs_WT")

############################################################
# 17. GSEA KEGG
############################################################

run_kegg_gsea <- function(gene_list, output_prefix) {
  
  kegg_res <- gseKEGG(
    geneList = gene_list,
    organism = "mmu",
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

gsea_kegg_het <- run_kegg_gsea(gene_list_het, "het_vs_WT")
gsea_kegg_KO <- run_kegg_gsea(gene_list_KO, "KO_vs_WT")

############################################################
# 18. GSEA Reactome
############################################################

run_reactome_gsea <- function(gene_list, output_prefix) {
  
  reactome_res <- gsePathway(
    geneList = gene_list,
    organism = "mouse",
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

gsea_reactome_het <- run_reactome_gsea(gene_list_het, "het_vs_WT")
gsea_reactome_KO <- run_reactome_gsea(gene_list_KO, "KO_vs_WT")

############################################################
# 19. GSEA GO Biological Process
############################################################

run_go_bp_gsea <- function(gene_list, output_prefix) {
  
  go_res <- gseGO(
    geneList = gene_list,
    OrgDb = org.Mm.eg.db,
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

gsea_go_het <- run_go_bp_gsea(gene_list_het, "het_vs_WT")
gsea_go_KO <- run_go_bp_gsea(gene_list_KO, "KO_vs_WT")


 ############################################################
# 20A. Helper function for fgsea Hallmark dotplots
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

p_hallmark_het <- plot_fgsea_dotplot(
  gsea_hallmark_het$export,
  "het_vs_WT",
  "Hallmark GSEA: +/- vs WT"
)

p_hallmark_KO <- plot_fgsea_dotplot(
  gsea_hallmark_KO$export,
  "KO_vs_WT",
  "Hallmark GSEA: -/- vs WT"
)

############################################################
# 20B. Hallmark enrichment curve plots
############################################################

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

plot_top_fgsea_enrichment(
  gsea_hallmark_het,
  gene_list_het,
  "het_vs_WT"
)

plot_top_fgsea_enrichment(
  gsea_hallmark_KO,
  gene_list_KO,
  "KO_vs_WT"
)
############################################################
# General dotplot function for clusterProfiler-style GSEA results
# Works for KEGG, Reactome, and GO
############################################################

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

############################################################
# 20C. Save KEGG, Reactome, and GO dotplots
############################################################
# This section uses one general dotplot function but applies
# different formatting settings for each gene-set database.
#
# Why settings differ:
#   KEGG terms are usually shorter.
#   Reactome terms can be moderately long.
#   GO Biological Process terms are often very long and redundant.
#
# Important parameters:
#   show_n:
#     Number of top categories shown.
#
#   label_width:
#     Controls line-wrapping of long pathway names.
#
#   fig_width / fig_height:
#     Controls saved PNG dimensions.
#
#   text_size:
#     Controls y-axis pathway label size.
############################################################

gsea_plot_settings <- list(
  KEGG = list(
    show_n = 10,
    label_width = 40,
    fig_width = 11,
    fig_height = 7,
    text_size = 8
  ),
  Reactome = list(
    show_n = 10,
    label_width = 45,
    fig_width = 13,
    fig_height = 7,
    text_size = 8
  ),
  GO_BP = list(
    show_n = 8,
    label_width = 35,
    fig_width = 14,
    fig_height = 8,
    text_size = 7
  )
)

gsea_objects <- list(
  het_vs_WT = list(
    KEGG = gsea_kegg_het,
    Reactome = gsea_reactome_het,
    GO_BP = gsea_go_het
  ),
  KO_vs_WT = list(
    KEGG = gsea_kegg_KO,
    Reactome = gsea_reactome_KO,
    GO_BP = gsea_go_KO
  )
)

gsea_titles <- list(
  het_vs_WT = list(
    KEGG = "KEGG GSEA: +/- vs WT",
    Reactome = "Reactome GSEA: +/- vs WT",
    GO_BP = "GO Biological Process GSEA: +/- vs WT"
  ),
  KO_vs_WT = list(
    KEGG = "KEGG GSEA: -/- vs WT",
    Reactome = "Reactome GSEA: -/- vs WT",
    GO_BP = "GO Biological Process GSEA: -/- vs WT"
  )
)

for (contrast_name in names(gsea_objects)) {
  
  for (database_name in names(gsea_objects[[contrast_name]])) {
    
    settings <- gsea_plot_settings[[database_name]]
    
    save_gsea_dotplot(
      gsea_result = gsea_objects[[contrast_name]][[database_name]],
      output_prefix = contrast_name,
      database_name = database_name,
      title = gsea_titles[[contrast_name]][[database_name]],
      show_n = settings$show_n,
      label_width = settings$label_width,
      fig_width = settings$fig_width,
      fig_height = settings$fig_height,
      text_size = settings$text_size
    )
  }
}

############################################################
# 20D. Running enrichment plots for clusterProfiler GSEA
############################################################
# These plots show the running enrichment score for the top pathway.
#
# The top pathway is selected by smallest adjusted p-value.
# Invalid rows where ID or p.adjust is NA are removed before choosing
# the top pathway. This avoids errors when GO/KEGG/Reactome return
# partially invalid GSEA rows.
############################################################

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

top_gseaplot_titles <- list(
  het_vs_WT = list(
    KEGG = "Top KEGG GSEA +/- vs WT",
    Reactome = "Top Reactome GSEA +/- vs WT",
    GO_BP = "Top GO BP GSEA +/- vs WT"
  ),
  KO_vs_WT = list(
    KEGG = "Top KEGG GSEA -/- vs WT",
    Reactome = "Top Reactome GSEA -/- vs WT",
    GO_BP = "Top GO BP GSEA -/- vs WT"
  )
)

for (contrast_name in names(gsea_objects)) {
  
  for (database_name in names(gsea_objects[[contrast_name]])) {
    
    save_top_gseaplot(
      gsea_result = gsea_objects[[contrast_name]][[database_name]],
      output_prefix = contrast_name,
      database_name = database_name,
      title_prefix = top_gseaplot_titles[[contrast_name]][[database_name]]
    )
  }
}

############################################################
# Save session info for reproducibility
############################################################

writeLines(
  capture.output(sessionInfo()),
  file.path(output_root, "results", "sessionInfo.txt")
)
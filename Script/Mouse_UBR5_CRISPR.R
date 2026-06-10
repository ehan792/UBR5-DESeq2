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


#Folder organization
output_dirs <- c(
  here("results"),
  here("results", "deseq2"),
  here("results", "deseq2", "het_vs_WT"),
  here("results", "deseq2", "KO_vs_WT"),
  here("results", "qc"),
  here("results", "gsea"),
  here("results", "gsea", "het_vs_WT"),
  here("results", "gsea", "KO_vs_WT"),
  here("figures"),
  here("figures", "qc"),
  here("figures", "volcano"),
  here("figures", "volcano", "het_vs_WT"),
  here("figures", "volcano", "KO_vs_WT"),
  here("figures", "gsea"),
  here("figures", "gsea", "het_vs_WT"),
  here("figures", "gsea", "KO_vs_WT")
)

walk(output_dirs, dir.create, showWarnings = FALSE, recursive = TRUE)
#helper fxns for folder organization
get_deseq_results_dir <- function(output_prefix) {
  if (output_prefix == "het_vs_WT") {
    return(here("results", "deseq2", "het_vs_WT"))
  } else if (output_prefix == "KO_vs_WT") {
    return(here("results", "deseq2", "KO_vs_WT"))
  } else {
    return(here("results", "deseq2"))
  }
}

get_volcano_fig_dir <- function(output_prefix) {
  if (output_prefix == "het_vs_WT") {
    return(here("figures", "volcano", "het_vs_WT"))
  } else if (output_prefix == "KO_vs_WT") {
    return(here("figures", "volcano", "KO_vs_WT"))
  } else {
    return(here("figures", "volcano"))
  }
}

get_gsea_results_dir <- function(output_prefix) {
  if (output_prefix == "het_vs_WT") {
    return(here("results", "gsea", "het_vs_WT"))
  } else if (output_prefix == "KO_vs_WT") {
    return(here("results", "gsea", "KO_vs_WT"))
  } else {
    return(here("results", "gsea"))
  }
}

get_gsea_fig_dir <- function(output_prefix) {
  if (output_prefix == "het_vs_WT") {
    return(here("figures", "gsea", "het_vs_WT"))
  } else if (output_prefix == "KO_vs_WT") {
    return(here("figures", "gsea", "KO_vs_WT"))
  } else {
    return(here("figures", "gsea"))
  }
}

############################################################
# Data Import
############################################################

counts_raw <- read.csv(
  here("Data/MouseKOCounts.csv"),
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



# 5. Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ genotype #design variable "formula"
)

############################################################
# 6. Filter low-count genes
############################################################
# Keep genes with at least 10 counts in at least 3 samples.
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds

############################################################
# 7. Run DESeq2
############################################################
dds <- DESeq(dds)

# Save normalized counts for visualization
norm_counts <- counts(dds, normalized = TRUE)

write.csv(
  as.data.frame(norm_counts) %>%
    rownames_to_column("ensembl_gene_id"),
  here("results", "deseq2", "normalized_counts.csv"),
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
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "genotype", returnData = TRUE)
pca_data$clean_name <- clean_sample_names[pca_data$name]

percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = genotype, label = clean_name)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw()

ggsave(
  here("figures", "qc", "PCA_vst_genotype.png"),
  p_pca,
  width = 7,
  height = 5,
  dpi = 300
)



# 9. Sample distance heatmap
sample_dists <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dists)

rownames(sample_dist_mat) <- colnames(vsd)
colnames(sample_dist_mat) <- colnames(vsd)

sample_order <- sample_info %>%
  rownames_to_column("sample") %>%
  arrange(genotype) %>%
  pull(sample)

sample_dist_mat_ordered <- sample_dist_mat[sample_order, sample_order]


rownames(sample_dist_mat_ordered) <- clean_sample_names[rownames(sample_dist_mat_ordered)]
colnames(sample_dist_mat_ordered) <- clean_sample_names[colnames(sample_dist_mat_ordered)]

png(
  here("figures", "qc", "sample_distance_heatmap.png"),
  width = 1800,
  height = 1600,
  res = 250
)

pheatmap(
  sample_dist_mat_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Sample distances using VST counts"
)

dev.off()

#QC heatmap, clustering
png(
  here("figures", "qc", "sample_distance_heatmap_clustered.png"),
  width = 1800,
  height = 1600,
  res = 250
)

pheatmap(
  sample_dist_mat_ordered,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Sample distances using VST counts, clustered"
)

dev.off()

############################################################
# 9A. Sample distance heatmap, grouped by genotype
############################################################

sample_dists <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dists)

rownames(sample_dist_mat) <- colnames(vsd)
colnames(sample_dist_mat) <- colnames(vsd)

sample_order <- sample_info %>%
  rownames_to_column("sample") %>%
  arrange(genotype, sample) %>%
  pull(sample)

sample_dist_mat_ordered <- sample_dist_mat[sample_order, sample_order]

rownames(sample_dist_mat_ordered) <- clean_sample_names[rownames(sample_dist_mat_ordered)]
colnames(sample_dist_mat_ordered) <- clean_sample_names[colnames(sample_dist_mat_ordered)]

png(
  here("figures", "qc", "sample_distance_heatmap_grouped.png"),
  width = 1800,
  height = 1600,
  res = 250
)

pheatmap(
  sample_dist_mat_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = c(3, 6),
  gaps_col = c(3, 6),
  main = "Sample distances using VST counts, grouped by genotype"
)

dev.off()

############################################################
# 9B. Sample distance heatmap, clustered QC version
############################################################

png(
  here("figures", "qc", "sample_distance_heatmap_clustered.png"),
  width = 1800,
  height = 1600,
  res = 250
)

pheatmap(
  sample_dist_mat_ordered,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Sample distances using VST counts, clustered"
)

dev.off()
# 10. Export DESeq2 results
export_deseq_result <- function(dds, contrast_vector, output_prefix) {
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
    here(
      get_deseq_results_dir(output_prefix),
      paste0(output_prefix, "_DESeq2_results.csv")
    ),
    row.names = FALSE
  )
  
  # LFC shrinkage for more stable effect-size visualization
  return(res_df)
}

resultsNames(dds)


# 11. Comparison to WT
res_het_vs_wt <- export_deseq_result(
  dds = dds,
  contrast_vector = c("genotype", "het", "WT"),
  output_prefix = "het_vs_WT"
)

res_KO_vs_wt <- export_deseq_result(
  dds = dds,
  contrast_vector = c("genotype", "KO", "WT"),
  output_prefix = "KO_vs_WT"
)

############################################################
# 12. Shrunk log2FC tables for visualization
############################################################

# Check exact coefficient names:
resultsNames(dds)

# These names should usually be:
# "Intercept"
# "genotype_het_vs_WT"
# "genotype_KO_vs_WT"

res_het_shrunk <- lfcShrink(
  dds,
  coef = "genotype_het_vs_WT",
  type = "apeglm"
)

res_KO_shrunk <- lfcShrink(
  dds,
  coef = "genotype_KO_vs_WT",
  type = "apeglm"
)

res_het_shrunk_df <- as.data.frame(res_het_shrunk) %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_annot, by = "ensembl_gene_id") %>%
  arrange(padj)

res_KO_shrunk_df <- as.data.frame(res_KO_shrunk) %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_annot, by = "ensembl_gene_id") %>%
  arrange(padj)

write.csv(
  res_het_shrunk_df,
  here(
    get_deseq_results_dir("het_vs_WT"),
    "het_vs_WT_DESeq2_shrunkLFC.csv"
  ),
  row.names = FALSE
)

write.csv(
  res_KO_shrunk_df,
  here(
    get_deseq_results_dir("KO_vs_WT"),
    "KO_vs_WT_DESeq2_shrunkLFC.csv"
  ),
  row.names = FALSE
)

summarize_de <- function(res_df, output_prefix, alpha = 0.05, lfc_cutoff = 1) {
  
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
    here(
      get_deseq_results_dir(output_prefix),
      paste0(output_prefix, "_DE_summary.csv")
    ),
    row.names = FALSE
  )
  
  return(summary_df)
}

de_summary_het <- summarize_de(res_het_vs_wt, "het_vs_WT")
de_summary_KO <- summarize_de(res_KO_vs_wt, "KO_vs_WT")

############################################################
# 13. Volcano plots
############################################################

png(
  here(get_volcano_fig_dir("het_vs_WT"), "volcano_het_vs_WT.png"),
  width = 2000,
  height = 1800,
  res = 250
)

EnhancedVolcano(
  res_het_shrunk_df,
  lab = res_het_shrunk_df$external_gene_name,
  x = "log2FoldChange",
  y = "padj",
  title = "+/- vs WT",
  subtitle = "DESeq2 with apeglm-shrunk log2FC",
  pCutoff = 0.05,
  FCcutoff = 1
)

dev.off()

png(
  here(get_volcano_fig_dir("KO_vs_WT"), "volcano_KO_vs_WT.png"),
  width = 2000,
  height = 1800,
  res = 250
)

EnhancedVolcano(
  res_KO_shrunk_df,
  lab = res_KO_shrunk_df$external_gene_name,
  x = "log2FoldChange",
  y = "padj",
  title = "-/- vs WT",
  subtitle = "DESeq2 with apeglm-shrunk log2FC",
  pCutoff = 0.05,
  FCcutoff = 1
)

dev.off()



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
    here(
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
    here(
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
    here(
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
    nPerm = 10000
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
    here(
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
    here(
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
      here(
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
      here(
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
    here(
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

#GSEA specific display settings
#KEGG
p_kegg_het <- save_gsea_dotplot(
  gsea_kegg_het,
  "het_vs_WT",
  "KEGG",
  "KEGG GSEA: +/- vs WT",
  show_n = 10,
  label_width = 40,
  fig_width = 11,
  fig_height = 7
)

p_kegg_KO <- save_gsea_dotplot(
  gsea_kegg_KO,
  "KO_vs_WT",
  "KEGG",
  "KEGG GSEA: -/- vs WT",
  show_n = 10,
  label_width = 40,
  fig_width = 11,
  fig_height = 7
)
#REACTOME
p_reactome_het <- save_gsea_dotplot(
  gsea_reactome_het,
  "het_vs_WT",
  "Reactome",
  "Reactome GSEA: +/- vs WT",
  show_n = 10,
  label_width = 45,
  fig_width = 13,
  fig_height = 7
)

p_reactome_KO <- save_gsea_dotplot(
  gsea_reactome_KO,
  "KO_vs_WT",
  "Reactome",
  "Reactome GSEA: -/- vs WT",
  show_n = 10,
  label_width = 45,
  fig_width = 13,
  fig_height = 7
)

#GO
p_go_het <- save_gsea_dotplot(
  gsea_go_het,
  "het_vs_WT",
  "GO_BP",
  "GO Biological Process GSEA: +/- vs WT",
  show_n = 8,
  label_width = 35,
  fig_width = 14,
  fig_height = 8,
  text_size = 7
)

p_go_KO <- save_gsea_dotplot(
  gsea_go_KO,
  "KO_vs_WT",
  "GO_BP",
  "GO Biological Process GSEA: -/- vs WT",
  show_n = 8,
  label_width = 35,
  fig_width = 14,
  fig_height = 8,
  text_size = 7
)


############################################################
# 20D. Running enrichment plots for clusterProfiler GSEA
############################################################

save_top_gseaplot <- function(gsea_result, output_prefix, database_name, title_prefix) {
  
  gsea_df <- as.data.frame(gsea_result)
  
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
    here(
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
save_top_gseaplot(
  gsea_kegg_het,
  "het_vs_WT",
  "KEGG",
  "Top KEGG GSEA +/- vs WT"
)

save_top_gseaplot(
  gsea_kegg_KO,
  "KO_vs_WT",
  "KEGG",
  "Top KEGG GSEA -/- vs WT"
)

save_top_gseaplot(
  gsea_reactome_het,
  "het_vs_WT",
  "Reactome",
  "Top Reactome GSEA +/- vs WT"
)

save_top_gseaplot(
  gsea_reactome_KO,
  "KO_vs_WT",
  "Reactome",
  "Top Reactome GSEA -/- vs WT"
)

save_top_gseaplot(
  gsea_go_het,
  "het_vs_WT",
  "GO_BP",
  "Top GO BP GSEA +/- vs WT"
)

save_top_gseaplot(
  gsea_go_KO,
  "KO_vs_WT",
  "GO_BP",
  "Top GO BP GSEA -/- vs WT"
)

############################################################
# Save session info for reproducibility
############################################################

writeLines(
  capture.output(sessionInfo()),
  here("results", "sessionInfo.txt")
)
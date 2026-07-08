################################################################################
# 02_cell_state_qc.R — QC assays: PCA, pooled PCA, ssGSEA, dispersion estimates
#
# Inputs:
#   results/data/deseq_all.rds from 01_import_deseq2.R
#
# Outputs:
#   results/tables/02_qc/<experiment>/
#   results/figures/02_qc/<experiment>/
#   results/tables/02_qc/pooledMouse/
#   results/figures/02_qc/pooledMouse/
################################################################################

if (!exists("PROJECT_ROOT")) {
  config_hits <- c(Sys.getenv("UBR5_CONFIG", unset = ""), file.path(getwd(), "00_config.R"), file.path(dirname(getwd()), "00_config.R"))
  config_hits <- config_hits[file.exists(config_hits)]
  if (!length(config_hits)) stop("Cannot find 00_config.R.")
  source(config_hits[1])
}

load_deseq_all <- function() {
  path <- p_data("deseq_all.rds")
  if (!file.exists(path)) stop("Missing ", path, ". Run 01_import_deseq2.R first.")
  readRDS(path)
}

sample_meta_df <- function(D) {
  # D$sample_info may carry row names after being saved/reloaded from RDS.
  # Remove them before any column_to_rownames() call to avoid tibble errors.
  D$sample_info %>%
    as.data.frame() %>%
    tibble::remove_rownames() %>%
    dplyr::select(sample, sample_label, condition, experiment)
}

plot_pca_one <- function(D, cfg, ntop = 2000) {
  mat <- assay(D$vsd)
  rv <- matrixStats::rowVars(mat)
  top <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pr <- prcomp(t(mat[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
  percent <- 100 * (pr$sdev^2 / sum(pr$sdev^2))

  pca_df <- tibble::as_tibble(pr$x[, seq_len(min(10, ncol(pr$x))), drop = FALSE], rownames = "sample") %>%
    left_join(sample_meta_df(D), by = "sample")

  write_tab(pca_df, "02_qc", cfg$exp_id, out_name(cfg, "PCA_coordinates.csv"))
  write_tab(tibble(PC = paste0("PC", seq_along(percent)), percent_variance = percent),
            "02_qc", cfg$exp_id, out_name(cfg, "PCA_percent_variance.csv"))

  p <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample_label)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
    scale_color_manual(values = COND_COLOURS, drop = FALSE) +
    labs(
      title = paste0(cfg$label, ": PCA"),
      subtitle = paste0("VST counts; top ", min(ntop, length(rv)), " variable genes"),
      x = paste0("PC1 (", round(percent[1], 1), "%)"),
      y = paste0("PC2 (", round(percent[2], 1), "%)")
    )
  save_fig(p, "02_qc", cfg$exp_id, out_name(cfg, "PCA.png"), width = 6.6, height = 5.2)

  p_scree <- tibble(PC = factor(paste0("PC", seq_len(min(10, length(percent)))), levels = paste0("PC", seq_len(min(10, length(percent))))),
                    percent_variance = percent[seq_len(min(10, length(percent)))]) %>%
    ggplot(aes(PC, percent_variance)) + geom_col() +
    labs(title = paste0(cfg$label, ": PCA scree"), x = NULL, y = "% variance")
  save_fig(p_scree, "02_qc", cfg$exp_id, out_name(cfg, "PCA_scree.png"), width = 5.6, height = 4.2)

  invisible(pca_df)
}

plot_sample_distance <- function(D, cfg) {
  mat <- assay(D$vsd)
  dist_mat <- as.matrix(dist(t(mat)))
  ann <- sample_meta_df(D) %>%
    dplyr::select(sample, condition) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("sample")
  dist_df <- as.data.frame(dist_mat) %>% rownames_to_column("sample")
  write_tab(dist_df, "02_qc", cfg$exp_id, out_name(cfg, "sample_distance_matrix.csv"))

  save_png("02_qc", cfg$exp_id, out_name(cfg, "sample_distance_heatmap.png"), width = 1800, height = 1600, res = 250, draw = function() {
    pheatmap::pheatmap(
      dist_mat,
      annotation_col = ann,
      annotation_row = ann,
      annotation_colors = list(condition = COND_COLOURS[names(COND_COLOURS) %in% levels(D$sample_info$condition)]),
      main = paste0(cfg$label, ": VST sample distances"),
      angle_col = "45"
    )
  })
}

export_dispersion_estimates <- function(D, cfg) {
  disp_tbl <- tibble(
    ensembl_gene_id = rownames(D$dds),
    baseMean = rowMeans(counts(D$dds, normalized = TRUE)),
    dispersion_gene_estimate = mcols(D$dds)$dispGeneEst,
    dispersion_fit = mcols(D$dds)$dispFit,
    dispersion_final = dispersions(D$dds)
  ) %>%
    left_join(D$annot, by = "ensembl_gene_id") %>%
    dplyr::select(ensembl_gene_id, any_of(c("entrezgene_id", "external_gene_name", "gene_biotype")), everything())
  write_tab(disp_tbl, "02_qc", cfg$exp_id, out_name(cfg, "dispersion_estimates.csv"))

  save_png("02_qc", cfg$exp_id, out_name(cfg, "dispersion_fit.png"), width = 1500, height = 1150, res = 240, draw = function() {
    plotDispEsts(D$dds, main = paste0(cfg$label, ": dispersion estimates"))
  })
}

run_ssgsea_one <- function(D, cfg) {
  message("  ssGSEA: ", cfg$exp_id)
  expr <- ens_to_entrez_mat(assay(D$vsd), cfg$orgdb_pkg)
  hallmark <- msig_to_pathways(fetch_msigdb(cfg$msig_species, "H"))
  scores <- safe_gsva(expr, hallmark, method = "ssgsea", min_size = 10, max_size = 500)

  score_long <- as.data.frame(scores) %>%
    rownames_to_column("gene_set") %>%
    pivot_longer(cols = -gene_set, names_to = "sample", values_to = "ssgsea_score") %>%
    left_join(sample_meta_df(D), by = "sample")

  write_tab(score_long, "02_qc", cfg$exp_id, out_name(cfg, "Hallmark_ssGSEA_scores_long.csv"))
  saveRDS(scores, p_data(out_name(cfg, "hallmark_ssgsea_scores.rds")))

  priority <- intersect(MPNST_PRIORITY_HALLMARK, rownames(scores))
  if (length(priority)) {
    plot_df <- score_long %>%
      filter(gene_set %in% priority) %>%
      mutate(gene_set = factor(gene_set, levels = priority))
    p <- ggplot(plot_df, aes(condition, ssgsea_score, color = condition)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.15) +
      geom_jitter(width = 0.12, height = 0, size = 2) +
      facet_wrap(~ gene_set, scales = "free_y") +
      scale_color_manual(values = COND_COLOURS, drop = FALSE) +
      labs(title = paste0(cfg$label, ": priority Hallmark ssGSEA"), x = NULL, y = "ssGSEA score") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    save_fig(p, "02_qc", cfg$exp_id, out_name(cfg, "priority_Hallmark_ssGSEA.png"), width = 12, height = 8)
  }
  scores
}

run_pooled_mouse_pca <- function(deseq_all) {
  mouse_ids <- intersect(c("CRISPR", "mouseKD"), names(deseq_all))
  if (length(mouse_ids) < 2) {
    message("Skipping pooled mouse PCA: need both CRISPR and mouseKD.")
    return(invisible(NULL))
  }
  mats <- lapply(mouse_ids, function(eid) assay(deseq_all[[eid]]$vsd))
  common <- Reduce(intersect, lapply(mats, rownames))
  if (length(common) < 1000) warning("Pooled mouse PCA has few common genes: ", length(common))
  mat <- do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE]))

  meta <- bind_rows(lapply(mouse_ids, function(eid) sample_meta_df(deseq_all[[eid]]))) %>%
    mutate(batch_experiment = experiment)
  rownames(meta) <- meta$sample
  mat <- mat[, meta$sample, drop = FALSE]

  mat_for_pca <- mat
  if (requireNamespace("limma", quietly = TRUE)) {
    mat_for_pca <- limma::removeBatchEffect(mat, batch = meta$batch_experiment)
  }

  rv <- matrixStats::rowVars(mat_for_pca)
  top <- order(rv, decreasing = TRUE)[seq_len(min(3000, length(rv)))]
  pr <- prcomp(t(mat_for_pca[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
  percent <- 100 * pr$sdev^2 / sum(pr$sdev^2)
  pca_df <- tibble::as_tibble(pr$x[, 1:5, drop = FALSE], rownames = "sample") %>% left_join(meta, by = "sample")
  write_tab(pca_df, "02_qc", "pooledMouse", file = "pooledMouse__PCA_coordinates.csv")

  p <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = experiment, label = sample_label)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
    scale_color_manual(values = COND_COLOURS, drop = FALSE) +
    labs(
      title = "Pooled mouse PCA",
      subtitle = "Common genes; batch-removed for visualization if limma is available",
      x = paste0("PC1 (", round(percent[1], 1), "%)"),
      y = paste0("PC2 (", round(percent[2], 1), "%)")
    )
  save_fig(p, "02_qc", "pooledMouse", file = "pooledMouse__PCA_batch_removed.png", width = 7.2, height = 5.4)
}

deseq_all <- load_deseq_all()
ssgsea_all <- list()
for (eid in names(deseq_all)) {
  cfg <- EXPERIMENTS[[eid]]
  D <- deseq_all[[eid]]
  message("\n── 02 QC: ", cfg$label, " ──")
  plot_pca_one(D, cfg)
  plot_sample_distance(D, cfg)
  export_dispersion_estimates(D, cfg)
  ssgsea_all[[eid]] <- run_ssgsea_one(D, cfg)
}
run_pooled_mouse_pca(deseq_all)
saveRDS(ssgsea_all, p_data("hallmark_ssgsea_all.rds"))

write_session("02_cell_state_qc")
message("\n02_cell_state_qc.R complete.")

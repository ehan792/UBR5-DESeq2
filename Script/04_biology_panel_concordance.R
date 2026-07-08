################################################################################
# 04_biology_panel_concordance.R — Biology gene panel and LFC concordance
#
# Outputs:
#   results/tables/04_biology_concordance/<experiment>/
#   results/figures/04_biology_concordance/<experiment>/
#   results/tables/04_biology_concordance/crossExperiment/
#   results/figures/04_biology_concordance/crossExperiment/
################################################################################

if (!exists("PROJECT_ROOT")) {
  config_hits <- c(Sys.getenv("UBR5_CONFIG", unset = ""), file.path(getwd(), "00_config.R"), file.path(dirname(getwd()), "00_config.R"))
  config_hits <- config_hits[file.exists(config_hits)]
  if (!length(config_hits)) stop("Cannot find 00_config.R.")
  source(config_hits[1])
}

path <- p_data("deseq_all.rds")
if (!file.exists(path)) stop("Missing ", path, ". Run 01_import_deseq2.R first.")
deseq_all <- readRDS(path)

# ── 1. Biology gene panel ────────────────────────────────────────────────────
default_mouse_panel <- tibble::tribble(
  ~class, ~gene_symbol,
  "UBR5", "Ubr5",
  "Interferon/viral mimicry", "Stat1",
  "Interferon/viral mimicry", "Stat2",
  "Interferon/viral mimicry", "Irf7",
  "Interferon/viral mimicry", "Isg15",
  "Interferon/viral mimicry", "Ifit1",
  "Interferon/viral mimicry", "Ifit2",
  "Interferon/viral mimicry", "Ifit3",
  "Interferon/viral mimicry", "Oas1a",
  "Interferon/viral mimicry", "Oas2",
  "Interferon/viral mimicry", "Mx1",
  "Interferon/viral mimicry", "Ddx58",
  "Interferon/viral mimicry", "Ifih1",
  "Interferon/viral mimicry", "Tbk1",
  "Interferon/viral mimicry", "Cxcl10",
  "NF-kB/inflammation", "Nfkb1",
  "NF-kB/inflammation", "Nfkbia",
  "NF-kB/inflammation", "Tnf",
  "NF-kB/inflammation", "Il6",
  "Cell cycle", "Mki67",
  "Cell cycle", "Ccnb1",
  "Cell cycle", "Cdk1",
  "Cell cycle", "Top2a",
  "Cell cycle", "Pcna",
  "MAPK/MPNST", "Nf1",
  "MAPK/MPNST", "Erbb2",
  "MAPK/MPNST", "Sox10",
  "MAPK/MPNST", "Pdgfra",
  "Chromatin/PRC2", "Eed",
  "Chromatin/PRC2", "Suz12",
  "Chromatin/PRC2", "Ezh2"
)

default_human_panel <- default_mouse_panel %>%
  mutate(gene_symbol = toupper(gene_symbol), gene_symbol = ifelse(gene_symbol == "UBR5", "UBR5", gene_symbol))

read_user_panel <- function() {
  if (!file.exists(BIOLOGY_PANEL_CSV)) return(NULL)
  x <- read.csv(BIOLOGY_PANEL_CSV, stringsAsFactors = FALSE, check.names = FALSE)
  names(x) <- str_replace_all(names(x), "\\s+", "_")
  if (all(c("Class", "Gene") %in% names(x))) {
    x <- x %>% transmute(class = Class, gene_symbol = Gene)
  } else if (all(c("class", "gene_symbol") %in% names(x))) {
    x <- x %>% transmute(class = class, gene_symbol = gene_symbol)
  } else {
    warning("Biology panel CSV found but does not have Class/Gene or class/gene_symbol columns. Using default panel.")
    return(NULL)
  }
  x %>% filter(!is.na(gene_symbol), gene_symbol != "") %>% distinct()
}

get_panel_for_experiment <- function(cfg) {
  user_panel <- read_user_panel()
  if (!is.null(user_panel)) return(user_panel)
  if (cfg$species == "Homo sapiens") default_human_panel else default_mouse_panel
}

export_panel_tables_and_heatmap <- function(D, cfg) {
  panel <- get_panel_for_experiment(cfg) %>% distinct()
  write.csv(panel, p_tab("04_biology_concordance", cfg$exp_id, file = paste0(cfg$exp_id, "__biology_gene_panel_used.csv")), row.names = FALSE)

  # Expression heatmap using VST values. Match by gene symbol.
  annot_panel <- D$annot %>%
    dplyr::select(ensembl_gene_id, external_gene_name, any_of(c("entrezgene_id", "gene_biotype"))) %>%
    inner_join(panel, by = c("external_gene_name" = "gene_symbol")) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)

  if (!nrow(annot_panel)) {
    warning("No panel genes matched annotation for ", cfg$exp_id)
    return(invisible(NULL))
  }

  # Some panel genes may have been removed by the low-count filter in file 01,
  # or may not exist in this experiment/species after annotation. Only subset
  # genes that are actually present in the VST matrix; otherwise R throws
  # "subscript out of bounds".
  vsd_mat <- assay(D$vsd)
  annot_panel_present <- annot_panel %>%
    dplyr::filter(ensembl_gene_id %in% rownames(vsd_mat))

  annot_panel_missing <- annot_panel %>%
    dplyr::filter(!ensembl_gene_id %in% rownames(vsd_mat))

  if (nrow(annot_panel_missing)) {
    write.csv(
      annot_panel_missing,
      p_tab(
        "04_biology_concordance",
        cfg$exp_id,
        file = paste0(cfg$exp_id, "__biology_gene_panel_dropped_not_in_filtered_VST_matrix.csv")
      ),
      row.names = FALSE
    )
    message(
      "  Note: ", nrow(annot_panel_missing),
      " panel gene rows were not present in the filtered VST matrix for ",
      cfg$exp_id, "; exported dropped-gene list."
    )
  }

  if (!nrow(annot_panel_present)) {
    warning("No panel genes remained in the filtered VST matrix for ", cfg$exp_id)
    return(invisible(NULL))
  }

  expr <- vsd_mat[annot_panel_present$ensembl_gene_id, , drop = FALSE]
  rownames(expr) <- make.unique(paste(annot_panel_present$external_gene_name, annot_panel_present$class, sep = " | "))

  # Keep displayed sample order explicit. For CRISPR: KO_1–3, Het_1–3, WT_1–3.
  # For KD experiments, keep the original column order.
  sample_order <- sample_order_for_experiment(colnames(expr), cfg$exp_id)
  expr <- expr[, sample_order, drop = FALSE]

  ann_col <- D$sample_info %>%
    as.data.frame() %>%
    tibble::remove_rownames() %>%
    dplyr::select(sample, condition) %>%
    tibble::column_to_rownames("sample")
  ann_col <- ann_col[colnames(expr), , drop = FALSE]

  expr_z <- t(scale(t(expr)))
  expr_z[is.na(expr_z)] <- 0
  png(p_fig("04_biology_concordance", cfg$exp_id, file = paste0(cfg$exp_id, "__biology_panel_VST_zscore_heatmap.png")), width = 1900, height = max(1400, 90 * nrow(expr_z)), res = 250)
  pheatmap::pheatmap(
    expr_z,
    annotation_col = ann_col,
    show_colnames = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 7,
    main = paste0(cfg$label, ": biology panel VST z-score"),
    angle_col = "45"
  )
  dev.off()

  # Per-contrast panel LFC tables.
  for (contrast_id in names(D$results)) {
    tab <- D$results[[contrast_id]]$wald %>%
      inner_join(panel, by = c("external_gene_name" = "gene_symbol")) %>%
      arrange(class, padj) %>%
      dplyr::select(class, external_gene_name, ensembl_gene_id, any_of(c("entrezgene_id")), baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
    write.csv(tab, p_tab("04_biology_concordance", cfg$exp_id, file = paste0(cfg$exp_id, "__", contrast_id, "__biology_panel_DESeq2_Wald.csv")), row.names = FALSE)
  }
  invisible(panel)
}

# ── 2. DEG set and overlap helpers ───────────────────────────────────────────
make_deg_table <- function(D, cfg, padj_cutoff = 0.05, lfc_cutoff = 1) {
  bind_rows(lapply(names(D$results), function(contrast_id) {
    D$results[[contrast_id]]$wald %>%
      mutate(
        experiment = cfg$exp_id,
        contrast = contrast_id,
        deg_direction = case_when(
          padj < padj_cutoff & log2FoldChange >= lfc_cutoff ~ "up",
          padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "down",
          TRUE ~ "not_DEG"
        ),
        is_DEG = deg_direction != "not_DEG"
      ) %>%
      filter(is_DEG) %>%
      dplyr::select(experiment, contrast, deg_direction, ensembl_gene_id, external_gene_name, any_of(c("entrezgene_id")), log2FoldChange, stat, padj)
  }))
}

export_deg_tables <- function(deseq_all) {
  all_deg <- bind_rows(lapply(names(deseq_all), function(eid) make_deg_table(deseq_all[[eid]], EXPERIMENTS[[eid]])))
  write.csv(all_deg, p_tab("04_biology_concordance", "crossExperiment", file = "all_experiments__DEG_lists_padj0.05_lfc1.csv"), row.names = FALSE)

  deg_summary <- all_deg %>%
    group_by(experiment, contrast, deg_direction) %>%
    summarise(n_genes = n_distinct(ensembl_gene_id), .groups = "drop")
  write.csv(deg_summary, p_tab("04_biology_concordance", "crossExperiment", file = "all_experiments__DEG_count_summary.csv"), row.names = FALSE)
  all_deg
}

compute_pairwise_overlap <- function(all_deg, universe_n = NULL) {
  sets <- all_deg %>%
    mutate(set_id = paste(experiment, contrast, deg_direction, sep = "__")) %>%
    group_by(set_id) %>%
    summarise(genes = list(unique(ensembl_gene_id)), .groups = "drop")
  if (nrow(sets) < 2) return(tibble())
  if (is.null(universe_n)) universe_n <- length(unique(unlist(sets$genes)))
  rows <- list()
  for (i in seq_len(nrow(sets) - 1)) {
    for (j in seq((i + 1), nrow(sets))) {
      a <- sets$genes[[i]]; b <- sets$genes[[j]]
      inter <- intersect(a, b); uni <- union(a, b)
      # Fisher table against all tested genes approximated by observed union universe.
      m <- matrix(c(length(inter), length(setdiff(a, b)), length(setdiff(b, a)), max(0, universe_n - length(uni))), nrow = 2)
      ft <- tryCatch(fisher.test(m), error = function(e) list(estimate = NA_real_, p.value = NA_real_))
      rows[[length(rows) + 1]] <- tibble(
        set_a = sets$set_id[i], set_b = sets$set_id[j],
        n_a = length(a), n_b = length(b), overlap = length(inter), union = length(uni),
        jaccard = ifelse(length(uni) > 0, length(inter) / length(uni), NA_real_),
        odds_ratio = unname(ft$estimate), pvalue = ft$p.value
      )
    }
  }
  bind_rows(rows) %>% mutate(padj = p.adjust(pvalue, method = "BH")) %>% arrange(desc(jaccard), padj)
}

# ── 3. LFC concordance ───────────────────────────────────────────────────────
get_contrast_lfc <- function(D, eid, contrast_id) {
  D$results[[contrast_id]]$wald %>%
    transmute(
      ensembl_gene_id,
      external_gene_name,
      contrast_key = paste(eid, contrast_id, sep = "__"),
      log2FoldChange,
      padj
    )
}

plot_lfc_concordance_pair <- function(df_a, df_b, label_a, label_b, file_prefix) {
  joined <- df_a %>%
    dplyr::select(ensembl_gene_id, external_gene_name, lfc_a = log2FoldChange, padj_a = padj) %>%
    inner_join(df_b %>% dplyr::select(ensembl_gene_id, lfc_b = log2FoldChange, padj_b = padj), by = "ensembl_gene_id") %>%
    filter(is.finite(lfc_a), is.finite(lfc_b)) %>%
    mutate(
      sig_class = case_when(
        padj_a < 0.05 & padj_b < 0.05 ~ "both padj<0.05",
        padj_a < 0.05 ~ paste0(label_a, " only"),
        padj_b < 0.05 ~ paste0(label_b, " only"),
        TRUE ~ "not significant"
      ),
      label_gene = ifelse(padj_a < 0.05 & padj_b < 0.05 & abs(lfc_a) + abs(lfc_b) > 3, external_gene_name, NA_character_)
    )
  rho <- suppressWarnings(cor(joined$lfc_a, joined$lfc_b, method = "spearman", use = "complete.obs"))
  write.csv(joined, p_tab("04_biology_concordance", "crossExperiment", file = paste0(file_prefix, "__LFC_concordance_table.csv")), row.names = FALSE)

  p <- ggplot(joined, aes(lfc_a, lfc_b)) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed") +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dotted") +
    geom_point(aes(color = sig_class), alpha = 0.55, size = 1.2) +
    ggrepel::geom_text_repel(data = joined %>% filter(!is.na(label_gene)) %>% slice_head(n = 25), aes(label = label_gene), size = 3, max.overlaps = Inf) +
    labs(
      title = paste0("LFC concordance: ", label_a, " vs ", label_b),
      subtitle = paste0("Spearman rho = ", round(rho, 3), "; each point is a gene"),
      x = paste0(label_a, " log2FC"),
      y = paste0(label_b, " log2FC"),
      color = NULL
    )
  ggsave(p_fig("04_biology_concordance", "crossExperiment", file = paste0(file_prefix, "__LFC_concordance.png")), p, width = 7, height = 6, dpi = 300)
  tibble(pair = file_prefix, label_a = label_a, label_b = label_b, spearman_rho = rho, n_genes = nrow(joined))
}

run_concordance <- function(deseq_all) {
  contrasts <- list()
  for (eid in names(deseq_all)) {
    for (contrast_id in names(deseq_all[[eid]]$results)) {
      key <- paste(eid, contrast_id, sep = "__")
      contrasts[[key]] <- get_contrast_lfc(deseq_all[[eid]], eid, contrast_id)
    }
  }
  pairs <- list(
    CRISPR_KO_vs_Het = c("CRISPR__KO_vs_WT", "CRISPR__Het_vs_WT"),
    mouseKD_vs_CRISPR_KO = c("mouseKD__KD_vs_Control", "CRISPR__KO_vs_WT"),
    mouseKD_vs_CRISPR_Het = c("mouseKD__KD_vs_Control", "CRISPR__Het_vs_WT")
  )
  summaries <- list()
  for (nm in names(pairs)) {
    a <- pairs[[nm]][1]; b <- pairs[[nm]][2]
    if (all(c(a, b) %in% names(contrasts))) {
      summaries[[nm]] <- plot_lfc_concordance_pair(contrasts[[a]], contrasts[[b]], a, b, nm)
    }
  }
  out <- bind_rows(summaries)
  if (nrow(out)) write.csv(out, p_tab("04_biology_concordance", "crossExperiment", file = "all_LFC_concordance_summary.csv"), row.names = FALSE)
}

# ── Run ──────────────────────────────────────────────────────────────────────
for (eid in names(deseq_all)) {
  message("\n── 04 Biology panel: ", eid, " ──")
  export_panel_tables_and_heatmap(deseq_all[[eid]], EXPERIMENTS[[eid]])
}
all_deg <- export_deg_tables(deseq_all)
if (nrow(all_deg)) {
  # Use number of tested genes from the largest DESeq2 result as approximate universe.
  universe_n <- max(unlist(lapply(deseq_all, function(D) nrow(D$results[[1]]$wald))))
  overlap <- compute_pairwise_overlap(all_deg, universe_n = universe_n)
  write.csv(overlap, p_tab("04_biology_concordance", "crossExperiment", file = "all_experiments__pairwise_DEG_overlap_Fisher_Jaccard.csv"), row.names = FALSE)
}
run_concordance(deseq_all)

write_session("04_biology_panel_concordance")
message("\n04_biology_panel_concordance.R complete.")

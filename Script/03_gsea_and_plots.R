################################################################################
# 03_gsea_and_plots.R — GSEA, volcano plots, dotplots, pathway comparisons
#
# Outputs:
#   results/tables/03_gsea/<experiment>/
#   results/figures/03_gsea/<experiment>/
#   results/tables/03_gsea/mouseDose/
#   results/figures/03_gsea/mouseDose/
#
# Notes:
#   Dotplots are ordered by significance and only show gene sets with padj <= 0.10.
#   Leading edge exports contain gene symbols, not only Entrez IDs.
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

make_ranked_vector <- function(wald_df) {
  ranked <- wald_df %>%
    filter(!is.na(stat), !is.na(entrezgene_id), entrezgene_id != "") %>%
    mutate(entrezgene_id = as.character(entrezgene_id)) %>%
    group_by(entrezgene_id) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(stat))
  ranks <- ranked$stat
  names(ranks) <- ranked$entrezgene_id
  ranks[is.finite(ranks)]
}

run_fgsea_collection <- function(pathways, ranks, cfg, collection_name, min_size = 15, max_size = 500) {
  sizes <- vapply(pathways, length, integer(1))
  pathways <- pathways[sizes >= min_size & sizes <= max_size]
  if (!length(pathways)) return(tibble())
  set.seed(GLOBAL_SEED)
  fg <- fgsea::fgsea(pathways = pathways, stats = ranks, minSize = min_size, maxSize = max_size, eps = 0) %>%
    as_tibble() %>%
    arrange(padj) %>%
    mutate(
      collection = collection_name,
      leadingEdge_entrez = vapply(leadingEdge, function(x) paste(as.character(x), collapse = ";"), character(1))
    )
  all_ids <- unique(unlist(fg$leadingEdge))
  sym_map <- entrez_to_symbol(all_ids, cfg$orgdb_pkg)
  fg %>%
    mutate(
      leadingEdge_symbols = vapply(leadingEdge, function(x) {
        syms <- unname(sym_map[as.character(x)])
        syms <- syms[!is.na(syms) & syms != ""]
        paste(unique(syms), collapse = ";")
      }, character(1)),
      leadingEdge_count = vapply(leadingEdge, length, integer(1))
    ) %>%
    dplyr::select(collection, pathway, pval, padj, ES, NES, size, leadingEdge_count, leadingEdge_entrez, leadingEdge_symbols)
}

plot_gsea_dotplot <- function(gsea_df, cfg, contrast_id, collection_name, padj_cutoff = 0.10, n_show = 20) {
  df <- gsea_df %>%
    filter(!is.na(padj), padj <= padj_cutoff) %>%
    arrange(padj) %>%
    slice_head(n = n_show) %>%
    mutate(
      pathway_label = pathway %>%
        str_replace_all("^HALLMARK_", "") %>%
        str_replace_all("^GOBP_", "") %>%
        str_replace_all("^REACTOME_", "") %>%
        str_replace_all("^KEGG_", "") %>%
        str_replace_all("_", " ") %>%
        str_to_sentence(),
      pathway_label = factor(pathway_label, levels = rev(pathway_label))
    )
  if (!nrow(df)) {
    message("  No ", collection_name, " dotplot entries for ", cfg$exp_id, " ", contrast_id, " at padj<=", padj_cutoff)
    return(invisible(NULL))
  }
  p <- ggplot(df, aes(x = NES, y = pathway_label, size = size, color = padj)) +
    geom_point(alpha = 0.95) +
    scale_color_gradient(low = GSEA_PADJ_LOW_COLOUR, high = GSEA_PADJ_HIGH_COLOUR, limits = c(0, padj_cutoff), oob = scales::squish) +
    labs(
      title = paste0(cfg$exp_id, " ", contrast_id, ": ", collection_name, " GSEA"),
      subtitle = paste0("Ordered by padj; showing padj <= ", padj_cutoff),
      x = "Normalized enrichment score (NES)",
      y = NULL,
      color = "padj",
      size = "Gene set size"
    )
  safe_collection <- str_replace_all(collection_name, "[^A-Za-z0-9]+", "_")
  save_fig(p, "03_gsea", cfg$exp_id, out_name(cfg, paste0(contrast_id, "__", safe_collection, "__dotplot_padj0.10.png")),
           width = 9.5, height = max(5, 0.33 * nrow(df) + 2.2))
}

plot_volcano <- function(shr_df, cfg, contrast_id) {
  df <- shr_df %>%
    mutate(
      neg_log10_padj = -log10(padj),
      sig_class = case_when(
        padj < 0.05 & log2FoldChange >= 1  ~ "Up, padj<0.05 & LFC>=1",
        padj < 0.05 & log2FoldChange <= -1 ~ "Down, padj<0.05 & LFC<=-1",
        padj < 0.05 ~ "padj<0.05 only",
        TRUE ~ "Not significant"
      ),
      label_gene = ifelse(padj < 0.05 & abs(log2FoldChange) >= 1.5, external_gene_name, NA_character_)
    )
  top_labels <- df %>% filter(!is.na(label_gene), label_gene != "") %>% arrange(padj) %>% slice_head(n = 20)
  p <- ggplot(df, aes(log2FoldChange, neg_log10_padj)) +
    geom_point(aes(color = sig_class), alpha = 0.55, size = 1.2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
    ggrepel::geom_text_repel(data = top_labels, aes(label = label_gene), size = 3, max.overlaps = Inf) +
    labs(
      title = paste0(cfg$exp_id, " ", contrast_id, ": volcano"),
      subtitle = "x-axis uses apeglm shrunken LFC; y-axis uses gene-level padj",
      x = "Shrunken log2 fold change",
      y = "-log10 adjusted p-value",
      color = NULL
    )
  save_fig(p, "03_gsea", cfg$exp_id, out_name(cfg, paste0(contrast_id, "__volcano_apeglm.png")), width = 7.2, height = 6.0)
}

export_important_gene_sets <- function(gsea_combined, cfg, contrast_id, padj_cutoff = 0.10, n_per_collection = 25) {
  important <- gsea_combined %>%
    filter(!is.na(padj), padj <= padj_cutoff) %>%
    group_by(collection) %>%
    arrange(padj, .by_group = TRUE) %>%
    slice_head(n = n_per_collection) %>%
    ungroup() %>%
    arrange(collection, padj) %>%
    mutate(experiment = cfg$exp_id, contrast = contrast_id) %>%
    dplyr::select(experiment, contrast, collection, pathway, NES, ES, pval, padj, size, leadingEdge_count, leadingEdge_symbols, leadingEdge_entrez)
  write_tab(important, "03_gsea", cfg$exp_id, out_name(cfg, paste0(contrast_id, "__important_gene_sets_padj0.10_with_leading_edge_symbols.csv")))
  important
}

run_gsea_one_contrast <- function(D, cfg, contrast_id, collections) {
  wald <- D$results[[contrast_id]]$wald
  shr  <- D$results[[contrast_id]]$shrunk
  ranks <- make_ranked_vector(wald)
  if (length(ranks) < 1000) warning(cfg$exp_id, " ", contrast_id, ": ranked vector has only ", length(ranks), " genes.")

  write_tab(tibble(entrezgene_id = names(ranks), stat = as.numeric(ranks)),
            "03_gsea", cfg$exp_id, out_name(cfg, paste0(contrast_id, "__ranked_Wald_stat.csv")))

  plot_volcano(shr, cfg, contrast_id)

  gsea_by_collection <- list()
  for (collection_name in names(collections)) {
    message("  GSEA ", cfg$exp_id, " ", contrast_id, " ", collection_name)
    res <- run_fgsea_collection(collections[[collection_name]], ranks, cfg, collection_name)
    gsea_by_collection[[collection_name]] <- res
    safe_collection <- str_replace_all(collection_name, "[^A-Za-z0-9]+", "_")
    write_tab(res, "03_gsea", cfg$exp_id, out_name(cfg, paste0(contrast_id, "__", safe_collection, "__GSEA.csv")))
    plot_gsea_dotplot(res, cfg, contrast_id, collection_name, padj_cutoff = 0.10, n_show = 20)
  }

  combined <- bind_rows(gsea_by_collection)
  write_tab(combined, "03_gsea", cfg$exp_id, out_name(cfg, paste0(contrast_id, "__all_collections_GSEA.csv")))
  important <- export_important_gene_sets(combined, cfg, contrast_id)
  list(collections = gsea_by_collection, combined = combined, important = important)
}

run_all_gsea <- function(deseq_all) {
  out <- list()
  for (eid in names(deseq_all)) {
    cfg <- EXPERIMENTS[[eid]]
    D <- deseq_all[[eid]]
    message("\n── 03 GSEA + plots: ", cfg$label, " ──")
    collections <- get_gsea_collections(cfg)
    out[[eid]] <- list()
    for (ct in cfg$contrasts) {
      out[[eid]][[ct$name]] <- run_gsea_one_contrast(D, cfg, ct$name, collections)
    }
  }
  out
}

# ── Mouse dose / GSEA comparison helpers ─────────────────────────────────────
get_hallmark_long <- function(gsea_all) {
  rows <- list()
  for (eid in names(gsea_all)) {
    for (contrast_id in names(gsea_all[[eid]])) {
      hm <- gsea_all[[eid]][[contrast_id]]$collections$Hallmark
      if (!is.null(hm) && nrow(hm)) {
        rows[[paste(eid, contrast_id, sep = "__")]] <- hm %>%
          mutate(experiment = eid, contrast = contrast_id) %>%
          dplyr::select(experiment, contrast, pathway, NES, padj, leadingEdge_symbols, leadingEdge_entrez, size)
      }
    }
  }
  bind_rows(rows)
}

plot_mouse_priority_nes <- function(hm_long) {
  mouse_df <- hm_long %>%
    filter(experiment %in% c("mouseKD", "CRISPR"), pathway %in% MPNST_PRIORITY_HALLMARK) %>%
    mutate(
      dose_label = case_when(
        experiment == "mouseKD" ~ "mouseKD",
        contrast == "Het_vs_WT" ~ "CRISPR Het",
        contrast == "KO_vs_WT" ~ "CRISPR KO",
        TRUE ~ contrast
      ),
      dose_label = factor(dose_label, levels = c("mouseKD", "CRISPR Het", "CRISPR KO")),
      pathway_label = pathway %>% str_replace("^HALLMARK_", "") %>% str_replace_all("_", " ") %>% str_to_sentence(),
      sig_tier = case_when(padj <= 0.05 ~ "padj <= 0.05", padj <= 0.10 ~ "padj <= 0.10", padj <= 0.25 ~ "padj <= 0.25", TRUE ~ "not significant")
    )
  if (!nrow(mouse_df)) return(invisible(NULL))
  write_tab(mouse_df, "03_gsea", "mouseDose", file = "mouseDose__MPNST_priority_Hallmark_NES_table.csv")

  p <- ggplot(mouse_df, aes(dose_label, NES, fill = dose_label)) +
    geom_col(width = 0.72) +
    facet_wrap(~ pathway_label, scales = "free_y") +
    scale_fill_manual(values = MOUSE_DOSE_COLOURS, drop = FALSE) +
    labs(
      title = "Mouse dose-priority Hallmark NES",
      subtitle = "Bar height is each contrast's direct NES, not a NES subtraction; fill shows mouse perturbation/comparison",
      x = NULL,
      y = "NES",
      fill = "Comparison"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p, "03_gsea", "mouseDose", file = "mouseDose__MPNST_priority_Hallmark_NES_bars.png", width = 13, height = 8.5)
}

compute_mouse_dose_ordering <- function(hm_long) {
  dose_df <- hm_long %>%
    filter(experiment %in% c("mouseKD", "CRISPR")) %>%
    mutate(dose_rank = case_when(
      experiment == "mouseKD" ~ 1,
      experiment == "CRISPR" & contrast == "Het_vs_WT" ~ 2,
      experiment == "CRISPR" & contrast == "KO_vs_WT" ~ 3,
      TRUE ~ NA_real_
    )) %>%
    filter(!is.na(dose_rank))
  tau_tbl <- dose_df %>%
    group_by(pathway) %>%
    summarise(
      n_conditions = n_distinct(dose_rank),
      tau = ifelse(n_conditions >= 3, suppressWarnings(cor(dose_rank, NES, method = "kendall", use = "complete.obs")), NA_real_),
      NES_mouseKD = NES[dose_rank == 1][1],
      NES_Het = NES[dose_rank == 2][1],
      NES_KO = NES[dose_rank == 3][1],
      min_padj = min(padj, na.rm = TRUE),
      .groups = "drop"
    ) %>% arrange(desc(abs(tau)), min_padj)
  write_tab(tau_tbl, "03_gsea", "mouseDose", file = "mouseDose__Hallmark_Kendall_NES_ordering.csv")
  tau_tbl
}

compute_leading_edge_jaccard <- function(hm_long, padj_cutoff = 0.25) {
  sig <- hm_long %>%
    filter(experiment %in% c("mouseKD", "CRISPR"), !is.na(padj), padj <= padj_cutoff, leadingEdge_entrez != "") %>%
    mutate(condition_label = case_when(
      experiment == "mouseKD" ~ "mouseKD",
      contrast == "Het_vs_WT" ~ "CRISPR_Het",
      contrast == "KO_vs_WT" ~ "CRISPR_KO",
      TRUE ~ paste(experiment, contrast, sep = "_")
    ))
  pathways <- sig %>%
    group_by(pathway) %>%
    summarise(n_sig_conditions = n_distinct(condition_label), .groups = "drop") %>%
    filter(n_sig_conditions >= 2) %>% pull(pathway)
  rows <- list()
  for (pw in pathways) {
    sub <- sig %>% filter(pathway == pw)
    labs <- unique(sub$condition_label)
    if (length(labs) < 2) next
    pairs <- combn(labs, 2, simplify = FALSE)
    for (pa in pairs) {
      a <- sub %>% filter(condition_label == pa[1]) %>% pull(leadingEdge_entrez) %>% str_split(";") %>% unlist() %>% unique()
      b <- sub %>% filter(condition_label == pa[2]) %>% pull(leadingEdge_entrez) %>% str_split(";") %>% unlist() %>% unique()
      inter <- intersect(a, b); uni <- union(a, b)
      rows[[length(rows) + 1]] <- tibble(
        pathway = pw,
        contrast_a = pa[1], contrast_b = pa[2],
        n_a = length(a), n_b = length(b), overlap = length(inter), union = length(uni),
        jaccard = ifelse(length(uni) > 0, length(inter) / length(uni), NA_real_),
        overlap_entrez = paste(inter, collapse = ";")
      )
    }
  }
  out <- bind_rows(rows) %>% arrange(desc(jaccard), pathway)
  if (nrow(out)) write_tab(out, "03_gsea", "mouseDose", file = "mouseDose__leading_edge_Jaccard_padj0.25.csv")
  out
}

# ── Run ──────────────────────────────────────────────────────────────────────
gsea_all <- run_all_gsea(deseq_all)
saveRDS(gsea_all, p_data("gsea_all.rds"))

hm_long <- get_hallmark_long(gsea_all)
if (nrow(hm_long)) {
  write_tab(hm_long, "03_gsea", "mouseDose", file = "all_experiments__Hallmark_NES_long.csv")
  plot_mouse_priority_nes(hm_long)
  compute_mouse_dose_ordering(hm_long)
  compute_leading_edge_jaccard(hm_long)
}

write_session("03_gsea_and_plots")
message("\n03_gsea_and_plots.R complete.")

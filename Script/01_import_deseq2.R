################################################################################
# 01_import_deseq2.R — Preprocess count tables, run DESeq2, export results
#
# Outputs:
#   results/tables/01_deseq2/<experiment>/
#     <experiment>__filter_summary.csv
#     <experiment>__DESeq2_normalized_counts_wide.csv
#     <experiment>__DESeq2_normalized_counts_long.csv
#     <experiment>__UBR5_percent_reference_per_sample.csv
#     <experiment>__UBR5_percent_reference_group_summary.csv
#     <experiment>__<contrast>__DESeq2_Wald.csv
#     <experiment>__<contrast>__DESeq2_apeglm_shrunk.csv
#     <experiment>__<contrast>__DE_summary.csv
#     <experiment>__<contrast>__DESeq2_Wald_by_symbol.csv
#
#   results/data/deseq_all.rds
#
# Important implementation note:
#   Saved objects do NOT include full cfg or OrgDb objects. This avoids the
#   C-stack error caused by serializing Bioconductor annotation databases.
#
# Downstream TF-inference (decoupler + CollecTRI, Python) note:
#   `<experiment>__<contrast>__DESeq2_Wald_by_symbol.csv` is Entrez-free,
#   gene-symbol-indexed, and deduplicated to one row per symbol — the format
#   decoupler's CollecTRI regulon matching expects. Load it with pandas,
#   e.g. `df.set_index("gene_symbol")["stat"]`, per experiment/contrast.
################################################################################

if (!exists("PROJECT_ROOT")) {
  config_hits <- c(
    Sys.getenv("UBR5_CONFIG", unset = ""),
    file.path(getwd(), "00_config.R"),
    file.path(dirname(getwd()), "00_config.R")
  )
  config_hits <- config_hits[file.exists(config_hits)]
  if (!length(config_hits)) stop("Cannot find 00_config.R. Set UBR5_CONFIG or run from the script/project folder.")
  source(config_hits[1])
}

read_count_file <- function(cfg) {
  if (!file.exists(cfg$file)) stop("Missing count file for ", cfg$exp_id, ": ", cfg$file)
  if (tolower(cfg$filetype) == "xlsx") {
    as.data.frame(readxl::read_excel(cfg$file))
  } else {
    read.csv(cfg$file, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

standardize_count_table <- function(raw, cfg) {
  if (!"ensembl_gene_id" %in% names(raw)) names(raw)[1] <- "ensembl_gene_id"
  if ("entrezgene" %in% names(raw) && !"entrezgene_id" %in% names(raw)) raw <- rename(raw, entrezgene_id = entrezgene)
  if (!"external_gene_name" %in% names(raw)) raw$external_gene_name <- NA_character_
  if (!"entrezgene_id" %in% names(raw)) raw$entrezgene_id <- NA_character_

  count_cols <- grep("^sample\\.", names(raw), value = TRUE)
  if (!length(count_cols)) stop(cfg$label, ": no sample count columns beginning with 'sample.' were found.")

  annot_cols <- intersect(c(
    "ensembl_gene_id", "entrezgene_id", "external_gene_name", "gene_biotype",
    "chromosome_name", "start_position", "end_position", "description"
  ), names(raw))

  counts_tbl <- raw %>%
    filter(!is.na(ensembl_gene_id), ensembl_gene_id != "") %>%
    mutate(across(all_of(count_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
    group_by(ensembl_gene_id) %>%
    summarise(across(all_of(count_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

  annot <- raw %>%
    filter(!is.na(ensembl_gene_id), ensembl_gene_id != "") %>%
    dplyr::select(all_of(annot_cols)) %>%
    group_by(ensembl_gene_id) %>%
    summarise(across(everything(), ~ dplyr::first(na.omit(.x), default = NA)), .groups = "drop") %>%
    mutate(entrezgene_id = as.character(entrezgene_id))

  mat <- counts_tbl %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("ensembl_gene_id") %>%
    as.matrix()
  mat <- round(mat)
  storage.mode(mat) <- "integer"
  if (anyNA(mat)) stop(cfg$label, ": NA values in count matrix after import.")
  if (any(mat < 0)) stop(cfg$label, ": negative counts found.")

  list(mat = mat, annot = annot, count_cols = count_cols)
}

make_sample_info <- function(mat, cfg) {
  condition <- rep(NA_character_, ncol(mat))
  for (nm in names(cfg$condition_map)) {
    hit <- stringr::str_detect(colnames(mat), cfg$condition_map[[nm]])
    condition[is.na(condition) & hit] <- nm
  }
  if (anyNA(condition)) {
    stop(cfg$label, ": unassigned samples: ", paste(colnames(mat)[is.na(condition)], collapse = ", "))
  }
  lvls <- c(cfg$ref_level, setdiff(names(cfg$condition_map), cfg$ref_level))
  si <- data.frame(
    sample = colnames(mat),
    sample_label = clean_sample_label(colnames(mat)),
    condition = factor(condition, levels = lvls),
    experiment = cfg$exp_id,
    stringsAsFactors = FALSE
  )
  rownames(si) <- si$sample
  si
}

export_normalized_counts <- function(norm, annot, sample_info, cfg) {
  norm_wide <- as.data.frame(norm) %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(annot, by = "ensembl_gene_id") %>%
    dplyr::select(
      ensembl_gene_id,
      any_of(c("entrezgene_id", "external_gene_name", "gene_biotype", "chromosome_name", "start_position", "end_position", "description")),
      everything()
    )

  write_tab(norm_wide, "01_deseq2", cfg$exp_id, out_name(cfg, "DESeq2_normalized_counts_wide.csv"))

  norm_long <- as.data.frame(norm) %>%
    rownames_to_column("ensembl_gene_id") %>%
    pivot_longer(cols = -ensembl_gene_id, names_to = "sample", values_to = "normalized_count") %>%
    left_join(annot %>% dplyr::select(ensembl_gene_id, any_of(c("entrezgene_id", "external_gene_name", "gene_biotype"))), by = "ensembl_gene_id") %>%
    left_join(sample_info %>% as.data.frame() %>% dplyr::select(sample, sample_label, condition), by = "sample") %>%
    mutate(experiment = cfg$exp_id) %>%
    dplyr::select(experiment, condition, sample, sample_label, ensembl_gene_id, any_of(c("entrezgene_id", "external_gene_name", "gene_biotype")), normalized_count)

  write_tab(norm_long, "01_deseq2", cfg$exp_id, out_name(cfg, "DESeq2_normalized_counts_long.csv"))

  invisible(list(wide = norm_wide, long = norm_long))
}

export_ubr5_percent_reference <- function(norm, annot, sample_info, cfg) {
  ids <- annot %>%
    filter(external_gene_name == cfg$ubr5_symbol) %>%
    pull(ensembl_gene_id) %>%
    unique()
  if (!length(ids)) {
    warning("Could not find ", cfg$ubr5_symbol, " in annotation for ", cfg$exp_id)
    return(invisible(NULL))
  }
  gid <- ids[1]
  per_sample <- tibble(
    experiment = cfg$exp_id,
    ensembl_gene_id = gid,
    gene_symbol = cfg$ubr5_symbol,
    sample = colnames(norm),
    normalized_count = as.numeric(norm[gid, ])
  ) %>%
    left_join(sample_info %>% as.data.frame() %>% dplyr::select(sample, sample_label, condition), by = "sample")

  ref_mean <- per_sample %>%
    filter(condition == cfg$ref_level) %>%
    summarise(reference_mean = mean(normalized_count, na.rm = TRUE)) %>%
    pull(reference_mean)

  per_sample <- per_sample %>%
    mutate(
      reference_condition = cfg$ref_level,
      reference_mean = ref_mean,
      percent_of_reference = 100 * normalized_count / ref_mean,
      remaining_fraction = normalized_count / ref_mean,
      loss_fraction = 1 - remaining_fraction
    )

  group_summary <- per_sample %>%
    group_by(experiment, gene_symbol, ensembl_gene_id, reference_condition, condition) %>%
    summarise(
      mean_normalized_count = mean(normalized_count, na.rm = TRUE),
      sd_normalized_count = sd(normalized_count, na.rm = TRUE),
      mean_percent_of_reference = mean(percent_of_reference, na.rm = TRUE),
      sd_percent_of_reference = sd(percent_of_reference, na.rm = TRUE),
      mean_remaining_fraction = mean(remaining_fraction, na.rm = TRUE),
      mean_loss_fraction = mean(loss_fraction, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )

  write_tab(per_sample, "01_deseq2", cfg$exp_id, out_name(cfg, "UBR5_percent_reference_per_sample.csv"))
  write_tab(group_summary, "01_deseq2", cfg$exp_id, out_name(cfg, "UBR5_percent_reference_group_summary.csv"))

  invisible(list(per_sample = per_sample, group_summary = group_summary))
}

# Deduplicated, gene-symbol-indexed Wald statistics for one contrast — the
# handoff format for the downstream Python TF-activity-inference step
# (decoupler + CollecTRI), which matches regulons by gene symbol rather than
# Ensembl/Entrez ID. Genes without a symbol are dropped; a symbol that maps
# to multiple Ensembl IDs keeps only the row with the largest |stat|,
# mirroring the Entrez-level collapsing used for GSEA ranking in
# 03_gsea_and_plots.R::make_ranked_vector().
export_tf_inference_stats <- function(wald_df, cfg, contrast_id) {
  by_symbol <- wald_df %>%
    filter(!is.na(external_gene_name), external_gene_name != "", !is.na(stat)) %>%
    group_by(external_gene_name) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(padj) %>%
    dplyr::rename(gene_symbol = external_gene_name) %>%
    dplyr::select(gene_symbol, ensembl_gene_id, any_of("entrezgene_id"), baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

  write_tab(by_symbol, "01_deseq2", cfg$exp_id, out_name(cfg, paste0(contrast_id, "__DESeq2_Wald_by_symbol.csv")))
  invisible(by_symbol)
}

run_deseq_one <- function(cfg, min_count = 10, min_samp = 3) {
  message("\n── 01 Import + DESeq2: ", cfg$label, " ──")
  raw <- read_count_file(cfg)
  std <- standardize_count_table(raw, cfg)
  si <- make_sample_info(std$mat, cfg)
  stopifnot(all(rownames(si) == colnames(std$mat)))

  message("  raw genes: ", nrow(std$mat), " | samples: ", ncol(std$mat))
  message("  groups: ", paste(names(table(si$condition)), as.integer(table(si$condition)), sep = "=", collapse = "; "))

  keep <- filter_by_group(std$mat, si, min_count = min_count, min_samp = min_samp)
  filter_tbl <- tibble(
    experiment = cfg$exp_id,
    label = cfg$label,
    genes_raw = nrow(std$mat),
    genes_kept = sum(keep),
    genes_removed = sum(!keep),
    min_count = min_count,
    min_samp = min_samp,
    filter_rule = paste0(">=", min_count, " counts in >=", min_samp, " samples within at least one condition")
  )
  write_tab(filter_tbl, "01_deseq2", cfg$exp_id, out_name(cfg, "filter_summary.csv"))

  dds <- DESeqDataSetFromMatrix(countData = std$mat[keep, , drop = FALSE], colData = si, design = ~ condition)
  dds$condition <- relevel(dds$condition, ref = cfg$ref_level)
  set.seed(GLOBAL_SEED)
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = TRUE)
  norm <- counts(dds, normalized = TRUE)

  export_normalized_counts(norm, std$annot, si, cfg)
  ubr5_reference <- export_ubr5_percent_reference(norm, std$annot, si, cfg)

  contrast_results <- list()
  for (ct in cfg$contrasts) {
    contrast_id <- ct$name
    coef_name <- paste0("condition_", ct$treat, "_vs_", ct$ctrl)
    if (!coef_name %in% resultsNames(dds)) {
      stop(cfg$exp_id, ": coefficient not found: ", coef_name, "\nAvailable: ", paste(resultsNames(dds), collapse = ", "))
    }

    wald_df <- as.data.frame(results(dds, contrast = c("condition", ct$treat, ct$ctrl), alpha = 0.05)) %>%
      rownames_to_column("ensembl_gene_id") %>%
      left_join(std$annot, by = "ensembl_gene_id") %>%
      arrange(padj)

    set.seed(GLOBAL_SEED)
    shr_df <- as.data.frame(lfcShrink(dds, coef = coef_name, type = "apeglm")) %>%
      rownames_to_column("ensembl_gene_id") %>%
      left_join(std$annot, by = "ensembl_gene_id") %>%
      arrange(padj)

    write_tab(wald_df, "01_deseq2", cfg$exp_id, out_name(cfg, paste0(contrast_id, "__DESeq2_Wald.csv")))
    write_tab(shr_df, "01_deseq2", cfg$exp_id, out_name(cfg, paste0(contrast_id, "__DESeq2_apeglm_shrunk.csv")))
    export_tf_inference_stats(wald_df, cfg, contrast_id)

    de_summary <- tibble(
      experiment = cfg$exp_id,
      contrast = contrast_id,
      tested_nonNA_padj = sum(!is.na(wald_df$padj)),
      padj_lt_0_10 = sum(wald_df$padj < 0.10, na.rm = TRUE),
      padj_lt_0_05 = sum(wald_df$padj < 0.05, na.rm = TRUE),
      padj_lt_0_01 = sum(wald_df$padj < 0.01, na.rm = TRUE),
      up_padj_0_05 = sum(wald_df$padj < 0.05 & wald_df$log2FoldChange > 0, na.rm = TRUE),
      down_padj_0_05 = sum(wald_df$padj < 0.05 & wald_df$log2FoldChange < 0, na.rm = TRUE),
      up_padj_0_05_lfc1 = sum(wald_df$padj < 0.05 & wald_df$log2FoldChange >= 1, na.rm = TRUE),
      down_padj_0_05_lfc1 = sum(wald_df$padj < 0.05 & wald_df$log2FoldChange <= -1, na.rm = TRUE)
    )
    write_tab(de_summary, "01_deseq2", cfg$exp_id, out_name(cfg, paste0(contrast_id, "__DE_summary.csv")))

    contrast_results[[contrast_id]] <- list(
      wald = wald_df,
      shrunk = shr_df,
      coef = coef_name,
      contrast = ct
    )
    message("  ", contrast_id, ": ", de_summary$padj_lt_0_05, " genes at padj<0.05")
  }

  # Lightweight bundle: no cfg, no OrgDb object.
  bundle <- list(
    experiment = cfg$exp_id,
    label = cfg$label,
    species = cfg$species,
    orgdb_pkg = cfg$orgdb_pkg,
    msig_species = cfg$msig_species,
    kegg_org = cfg$kegg_org,
    ref_level = cfg$ref_level,
    contrasts = cfg$contrasts,
    ubr5_symbol = cfg$ubr5_symbol,
    raw_dim = dim(std$mat),
    keep = keep,
    dds = dds,
    vsd = vsd,
    norm = norm,
    annot = std$annot,
    sample_info = si,
    ubr5_reference = ubr5_reference,
    results = contrast_results
  )
  saveRDS(bundle, p_data(out_name(cfg, "deseq_bundle.rds")))
  bundle
}


export_combined_ubr5_reference_tables <- function(deseq_all) {
  per_sample <- purrr::imap_dfr(deseq_all, function(D, eid) {
    if (is.null(D$ubr5_reference) || is.null(D$ubr5_reference$per_sample)) return(tibble())
    D$ubr5_reference$per_sample
  })

  phenotype_summary <- purrr::imap_dfr(deseq_all, function(D, eid) {
    if (is.null(D$ubr5_reference) || is.null(D$ubr5_reference$group_summary)) return(tibble())
    D$ubr5_reference$group_summary
  })

  if (nrow(per_sample)) {
    per_sample <- per_sample %>%
      mutate(
        ratio_to_reference = remaining_fraction,
        percent_of_reference = percent_of_reference
      ) %>%
      dplyr::select(
        experiment, condition, sample, sample_label,
        gene_symbol, ensembl_gene_id,
        normalized_count, reference_condition, reference_mean,
        ratio_to_reference, percent_of_reference, loss_fraction
      )

    write_tab(per_sample, "01_deseq2", file = "all_experiments__UBR5_normalized_expression_ratio_to_reference_per_sample.csv")
  }

  if (nrow(phenotype_summary)) {
    phenotype_summary <- phenotype_summary %>%
      mutate(
        mean_ratio_to_reference = mean_remaining_fraction,
        mean_percent_of_reference = mean_percent_of_reference
      ) %>%
      dplyr::select(
        experiment, condition, gene_symbol, ensembl_gene_id, reference_condition,
        mean_normalized_count, sd_normalized_count,
        mean_ratio_to_reference, mean_percent_of_reference,
        sd_percent_of_reference, mean_loss_fraction, n
      )

    write_tab(phenotype_summary, "01_deseq2", file = "all_experiments__UBR5_normalized_expression_ratio_to_reference_by_phenotype.csv")
  }

  invisible(list(per_sample = per_sample, phenotype_summary = phenotype_summary))
}

available <- EXPERIMENTS[vapply(EXPERIMENTS, function(x) file.exists(x$file), logical(1))]
missing <- setdiff(names(EXPERIMENTS), names(available))
if (length(missing)) warning("Skipping missing count files: ", paste(missing, collapse = ", "))
if (!length(available)) stop("No count files found. Check DATA_DIR and COUNT_FILES in 00_config.R: ", DATA_DIR)

deseq_all <- list()
for (eid in names(available)) {
  deseq_all[[eid]] <- run_deseq_one(available[[eid]])
}
export_combined_ubr5_reference_tables(deseq_all)
saveRDS(deseq_all, p_data("deseq_all.rds"))

write_session("01_import_deseq2")
message("\n01_import_deseq2.R complete: ", paste(names(deseq_all), collapse = ", "))

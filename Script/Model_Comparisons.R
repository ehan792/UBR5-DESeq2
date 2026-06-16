############################################################
# Compare GSEA results across models
############################################################
# Comparisons:
#   1. Mouse KO vs mouse KD vs mouse het
#   2. Human KD vs mouse KD
#
# This script reads completed GSEA result CSVs from:
#   CRISPR/results/gsea/
#   KD/results/gsea/
#
# and writes comparison outputs to:
#   Model_Comparisons/results/
#   Model_Comparisons/figures/
#
# Important:
#   This is pathway-level comparative profiling.
#   delta_NES values are descriptive effect-size differences,
#   not formal p-values for one model differing from another.
############################################################

library(tidyverse)
library(here)
library(pheatmap)

############################################################
# 1. Root-level output folders
############################################################

comparison_root <- here("Model_Comparisons")

comparison_results_dir <- file.path(comparison_root, "results")
comparison_figures_dir <- file.path(comparison_root, "figures")

comparison_dirs <- c(
  comparison_results_dir,
  file.path(comparison_results_dir, "mouse_KO_KD_het"),
  file.path(comparison_results_dir, "human_vs_mouse_KD"),
  
  comparison_figures_dir,
  file.path(comparison_figures_dir, "mouse_KO_KD_het"),
  file.path(comparison_figures_dir, "human_vs_mouse_KD")
)

walk(comparison_dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

############################################################
# 1B. Global significance and scatterplot filtering controls
############################################################
# primary_fdr_cutoff:
#   Used to create the general "significant" column in tables.
#
# scatter_fdr_cutoff:
#   Used to filter scatterplots and filtered scatter pathway lists.
#
# scatter_filter_mode options:
#
#   "none":
#     no FDR filtering
#
#   "either":
#     keep pathway if either comparison has padj < scatter_fdr_cutoff
#
#   "both":
#     keep pathway only if both comparisons have padj < scatter_fdr_cutoff
#
# Suggested:
#   primary_fdr_cutoff = 0.05 for strict table labels
#   scatter_fdr_cutoff = 0.25 for exploratory GSEA-style screening
#   scatter_fdr_cutoff = 0.05 for stricter plots/lists
############################################################

primary_fdr_cutoff <- 0.05
scatter_fdr_cutoff <- 0.05
scatter_filter_mode <- "either"

############################################################
# 2. Input file paths
############################################################
# These file names must match the output naming pattern from your
# CRISPR and KD scripts.
############################################################

gsea_files <- list(
  
  mouse_het = list(
    Hallmark = here("CRISPR", "results", "gsea", "het_vs_WT", "het_vs_WT_GSEA_Hallmark.csv"),
    KEGG     = here("CRISPR", "results", "gsea", "het_vs_WT", "het_vs_WT_GSEA_KEGG.csv"),
    Reactome = here("CRISPR", "results", "gsea", "het_vs_WT", "het_vs_WT_GSEA_Reactome.csv"),
    GO_BP    = here("CRISPR", "results", "gsea", "het_vs_WT", "het_vs_WT_GSEA_GO_BP.csv")
  ),
  
  mouse_KO = list(
    Hallmark = here("CRISPR", "results", "gsea", "KO_vs_WT", "KO_vs_WT_GSEA_Hallmark.csv"),
    KEGG     = here("CRISPR", "results", "gsea", "KO_vs_WT", "KO_vs_WT_GSEA_KEGG.csv"),
    Reactome = here("CRISPR", "results", "gsea", "KO_vs_WT", "KO_vs_WT_GSEA_Reactome.csv"),
    GO_BP    = here("CRISPR", "results", "gsea", "KO_vs_WT", "KO_vs_WT_GSEA_GO_BP.csv")
  ),
  
  mouse_KD = list(
    Hallmark = here("KD", "results", "gsea", "mouse_KD_vs_Control", "mouse_KD_vs_Control_GSEA_Hallmark.csv"),
    KEGG     = here("KD", "results", "gsea", "mouse_KD_vs_Control", "mouse_KD_vs_Control_GSEA_KEGG.csv"),
    Reactome = here("KD", "results", "gsea", "mouse_KD_vs_Control", "mouse_KD_vs_Control_GSEA_Reactome.csv"),
    GO_BP    = here("KD", "results", "gsea", "mouse_KD_vs_Control", "mouse_KD_vs_Control_GSEA_GO_BP.csv")
  ),
  
  human_KD = list(
    Hallmark = here("KD", "results", "gsea", "human_KD_vs_Control", "human_KD_vs_Control_GSEA_Hallmark.csv"),
    KEGG     = here("KD", "results", "gsea", "human_KD_vs_Control", "human_KD_vs_Control_GSEA_KEGG.csv"),
    Reactome = here("KD", "results", "gsea", "human_KD_vs_Control", "human_KD_vs_Control_GSEA_Reactome.csv"),
    GO_BP    = here("KD", "results", "gsea", "human_KD_vs_Control", "human_KD_vs_Control_GSEA_GO_BP.csv")
  )
)

############################################################
# 3. Safety check: confirm files exist
############################################################

all_input_files <- unlist(gsea_files)

missing_files <- all_input_files[!file.exists(all_input_files)]

if (length(missing_files) > 0) {
  stop(
    "These expected GSEA result files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

############################################################
# 4. Helper: clean cross-species pathway names
############################################################
# This is mainly for Reactome/KEGG cross-species matching.
# Human and mouse IDs often differ, so cleaned names are safer
# than direct ID matching for those databases.
############################################################

clean_pathway_name_for_cross_species <- function(pathway_name) {
  
  pathway_name %>%
    str_remove(" - Mus musculus \\(house mouse\\)$") %>%
    str_remove(" - Homo sapiens \\(human\\)$") %>%
    str_remove(" - mouse$") %>%
    str_remove(" - human$") %>%
    str_remove("\\s*\\(Mus musculus\\)$") %>%
    str_remove("\\s*\\(Homo sapiens\\)$") %>%
    str_squish()
}

############################################################
# 5. Standardize GSEA tables
############################################################
# Hallmark files come from fgsea:
#   pathway, NES, padj, leadingEdge
#
# KEGG/Reactome/GO files come from clusterProfiler:
#   ID, Description, NES, p.adjust, core_enrichment
#
# This function creates two keys:
#
#   within_species_key:
#     Used for mouse-only comparisons.
#     Since mouse het, KO, and KD all use mouse annotations,
#     pathway IDs are appropriate.
#
#   cross_species_key:
#     Used for human KD vs mouse KD.
#     Hallmark: use pathway name.
#     GO BP: use GO ID.
#     Reactome/KEGG: use cleaned pathway name.
############################################################

read_standardized_gsea <- function(file_path, model_name, database_name) {
  
  df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (database_name == "Hallmark") {
    
    standardized <- df %>%
      transmute(
        model = model_name,
        database = database_name,
        pathway_id = pathway,
        pathway_name = pathway,
        NES = NES,
        padj = padj,
        leading_edge = leadingEdge
      )
    
  } else {
    
    standardized <- df %>%
      transmute(
        model = model_name,
        database = database_name,
        pathway_id = ID,
        pathway_name = Description,
        NES = NES,
        padj = p.adjust,
        leading_edge = if ("core_enrichment" %in% colnames(df)) {
          core_enrichment
        } else {
          NA_character_
        }
      )
  }
  
  standardized <- standardized %>%
    filter(
      !is.na(pathway_id),
      !is.na(pathway_name),
      !is.na(NES),
      !is.na(padj),
      is.finite(NES),
      is.finite(padj)
    ) %>%
    mutate(
      within_species_key = pathway_id,
      
      cross_species_key = case_when(
        database == "Hallmark" ~ pathway_name,
        database == "GO_BP" ~ pathway_id,
        database %in% c("Reactome", "KEGG") ~ clean_pathway_name_for_cross_species(pathway_name),
        TRUE ~ pathway_name
      ),
      
      signed_log10padj = sign(NES) * -log10(pmax(padj, .Machine$double.xmin)),
      significant = padj < primary_fdr_cutoff
    )
  
  return(standardized)
}

read_database_for_all_models <- function(database_name) {
  
  bind_rows(
    read_standardized_gsea(
      gsea_files$mouse_het[[database_name]],
      "mouse_het",
      database_name
    ),
    read_standardized_gsea(
      gsea_files$mouse_KO[[database_name]],
      "mouse_KO",
      database_name
    ),
    read_standardized_gsea(
      gsea_files$mouse_KD[[database_name]],
      "mouse_KD",
      database_name
    ),
    read_standardized_gsea(
      gsea_files$human_KD[[database_name]],
      "human_KD",
      database_name
    )
  )
}

all_gsea_long <- bind_rows(
  read_database_for_all_models("Hallmark"),
  read_database_for_all_models("KEGG"),
  read_database_for_all_models("Reactome"),
  read_database_for_all_models("GO_BP")
)

write.csv(
  all_gsea_long,
  file.path(comparison_results_dir, "all_standardized_GSEA_results_long.csv"),
  row.names = FALSE
)

############################################################
# 6. Helper: make wide comparison table
############################################################
# comparison_type:
#   "within_species":
#     Uses within_species_key.
#     Best for mouse KO/KD/het comparison.
#
#   "cross_species":
#     Uses cross_species_key.
#     Best for human KD vs mouse KD.
############################################################

make_wide_comparison <- function(
    gsea_long,
    models_to_include,
    database_name,
    comparison_type = c("within_species", "cross_species")
) {
  
  comparison_type <- match.arg(comparison_type)
  
  key_col <- ifelse(
    comparison_type == "within_species",
    "within_species_key",
    "cross_species_key"
  )
  
  wide_df <- gsea_long %>%
    filter(
      model %in% models_to_include,
      database == database_name
    ) %>%
    mutate(pathway_key = .data[[key_col]]) %>%
    select(
      database,
      pathway_key,
      pathway_id,
      pathway_name,
      model,
      NES,
      padj,
      signed_log10padj,
      significant,
      leading_edge
    ) %>%
    group_by(database, pathway_key, model) %>%
    slice_max(order_by = abs(NES), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = c(database, pathway_key),
      names_from = model,
      values_from = c(
        pathway_id,
        pathway_name,
        NES,
        padj,
        signed_log10padj,
        significant,
        leading_edge
      ),
      names_sep = "_"
    )
  
  possible_name_cols <- paste0("pathway_name_", models_to_include)
  
  for (name_col in possible_name_cols) {
    if (!name_col %in% colnames(wide_df)) {
      wide_df[[name_col]] <- NA_character_
    }
  }
  
  wide_df <- wide_df %>%
    mutate(
      pathway_name = coalesce(!!!syms(possible_name_cols))
    )
  
  return(wide_df)
}

############################################################
# 7. Mouse-only comparison: KO vs KD vs het
############################################################
# This is the cleanest comparison because all three models are mouse.
#
# Heatmap order later will be:
#   KO, KD, het
############################################################

mouse_models <- c("mouse_het", "mouse_KO", "mouse_KD")

mouse_comparisons <- list()

for (database_name in c("Hallmark", "KEGG", "Reactome", "GO_BP")) {
  
  wide_df <- make_wide_comparison(
    gsea_long = all_gsea_long,
    models_to_include = mouse_models,
    database_name = database_name,
    comparison_type = "within_species"
  ) %>%
    mutate(
      delta_NES_KO_minus_het = NES_mouse_KO - NES_mouse_het,
      delta_NES_KD_minus_het = NES_mouse_KD - NES_mouse_het,
      delta_NES_KO_minus_KD = NES_mouse_KO - NES_mouse_KD,
      
      abs_delta_NES_KO_minus_het = abs(delta_NES_KO_minus_het),
      abs_delta_NES_KD_minus_het = abs(delta_NES_KD_minus_het),
      abs_delta_NES_KO_minus_KD = abs(delta_NES_KO_minus_KD),
      
      direction_pattern = case_when(
        NES_mouse_het > 0 & NES_mouse_KO > 0 & NES_mouse_KD > 0 ~ "all positive",
        NES_mouse_het < 0 & NES_mouse_KO < 0 & NES_mouse_KD < 0 ~ "all negative",
        NES_mouse_het > 0 & NES_mouse_KO > 0 & NES_mouse_KD < 0 ~ "KD opposite",
        NES_mouse_het > 0 & NES_mouse_KO < 0 & NES_mouse_KD > 0 ~ "KO opposite",
        NES_mouse_het < 0 & NES_mouse_KO > 0 & NES_mouse_KD > 0 ~ "het opposite",
        TRUE ~ "mixed or incomplete"
      )
    ) %>%
    arrange(desc(abs_delta_NES_KO_minus_het))
  
  mouse_comparisons[[database_name]] <- wide_df
  
  write.csv(
    wide_df,
    file.path(
      comparison_results_dir,
      "mouse_KO_KD_het",
      paste0("mouse_KO_KD_het_GSEA_comparison_", database_name, ".csv")
    ),
    row.names = FALSE
  )
}

mouse_comparison_all <- bind_rows(mouse_comparisons)

write.csv(
  mouse_comparison_all,
  file.path(
    comparison_results_dir,
    "mouse_KO_KD_het",
    "mouse_KO_KD_het_GSEA_comparison_ALL.csv"
  ),
  row.names = FALSE
)

############################################################
# 8. Human KD vs mouse KD comparison
############################################################
# Recommended primary comparison:
#   Hallmark and GO_BP
#
# Reactome and KEGG are exported as exploratory comparisons.
# Their IDs and names can differ across species, so they may have
# fewer matched rows or no matched rows.
############################################################

kd_models <- c("mouse_KD", "human_KD")

kd_cross_species_comparisons <- list()

for (database_name in c("Hallmark", "GO_BP", "Reactome", "KEGG")) {
  
  wide_df <- make_wide_comparison(
    gsea_long = all_gsea_long,
    models_to_include = kd_models,
    database_name = database_name,
    comparison_type = "cross_species"
  ) %>%
    mutate(
      delta_NES_human_minus_mouse = NES_human_KD - NES_mouse_KD,
      abs_delta_NES_human_minus_mouse = abs(delta_NES_human_minus_mouse),
      direction_pattern = case_when(
        NES_human_KD > 0 & NES_mouse_KD > 0 ~ "shared positive",
        NES_human_KD < 0 & NES_mouse_KD < 0 ~ "shared negative",
        NES_human_KD > 0 & NES_mouse_KD < 0 ~ "opposite direction",
        NES_human_KD < 0 & NES_mouse_KD > 0 ~ "opposite direction",
        TRUE ~ "mixed or incomplete"
      )
    ) %>%
    arrange(desc(abs_delta_NES_human_minus_mouse))
  
  kd_cross_species_comparisons[[database_name]] <- wide_df
  
  write.csv(
    wide_df,
    file.path(
      comparison_results_dir,
      "human_vs_mouse_KD",
      paste0("human_vs_mouse_KD_GSEA_comparison_", database_name, ".csv")
    ),
    row.names = FALSE
  )
}

kd_cross_species_all <- bind_rows(kd_cross_species_comparisons)

write.csv(
  kd_cross_species_all,
  file.path(
    comparison_results_dir,
    "human_vs_mouse_KD",
    "human_vs_mouse_KD_GSEA_comparison_ALL.csv"
  ),
  row.names = FALSE
)

############################################################
# 9. Plot helper: NES heatmap
############################################################

plot_nes_heatmap <- function(
    comparison_df,
    database_name,
    model_columns,
    output_file,
    title,
    top_n = 30,
    ranking_column = NULL,
    clean_colnames = NULL
) {
  
  if (is.null(ranking_column)) {
    ranking_column <- model_columns[1]
  }
  
  plot_df <- comparison_df %>%
    filter(database == database_name)
  
  if (ranking_column %in% colnames(plot_df)) {
    plot_df <- plot_df %>%
      filter(!is.na(.data[[ranking_column]])) %>%
      arrange(desc(abs(.data[[ranking_column]])))
  }
  
  plot_df <- plot_df %>%
    filter(
      if_all(
        all_of(model_columns),
        ~ !is.na(.) & is.finite(.)
      )
    ) %>%
    slice_head(n = top_n) %>%
    select(pathway_name, all_of(model_columns))
  
  if (nrow(plot_df) == 0) {
    message("No rows to plot for: ", database_name)
    return(NULL)
  }
  
  heatmap_mat <- plot_df %>%
    column_to_rownames("pathway_name") %>%
    as.matrix()
  
  if (!is.null(clean_colnames)) {
    if (length(clean_colnames) == ncol(heatmap_mat)) {
      colnames(heatmap_mat) <- clean_colnames
    }
  }
  
  rownames(heatmap_mat) <- str_wrap(rownames(heatmap_mat), width = 45)
  
  png(
    output_file,
    width = 2000,
    height = 2400,
    res = 250
  )
  
  pheatmap(
    heatmap_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = title,
    fontsize_row = 7,
    fontsize_col = 10
  )
  
  dev.off()
  
  return(heatmap_mat)
}

############################################################
# 10. Plot helper: pairwise NES scatter
############################################################
# This function can make either unfiltered or FDR-filtered
# scatterplots, controlled by:
#
#   padj_cutoff
#   filter_mode
#
# filter_mode:
#   "none":
#     no FDR filtering
#
#   "either":
#     keep pathways if either x or y comparison has padj < cutoff
#
#   "both":
#     keep pathways only if both x and y comparisons have padj < cutoff
############################################################

plot_pairwise_nes_scatter <- function(
    comparison_df,
    database_name,
    x_col,
    y_col,
    x_padj_col,
    y_padj_col,
    x_label,
    y_label,
    output_file,
    title,
    padj_cutoff = scatter_fdr_cutoff,
    filter_mode = scatter_filter_mode
) {
  
  filter_mode <- match.arg(filter_mode, choices = c("none", "either", "both"))
  
  plot_df <- comparison_df %>%
    filter(
      database == database_name,
      !is.na(.data[[x_col]]),
      !is.na(.data[[y_col]]),
      is.finite(.data[[x_col]]),
      is.finite(.data[[y_col]]),
      !is.na(.data[[x_padj_col]]),
      !is.na(.data[[y_padj_col]]),
      is.finite(.data[[x_padj_col]]),
      is.finite(.data[[y_padj_col]])
    )
  
  if (filter_mode == "either") {
    
    plot_df <- plot_df %>%
      filter(
        .data[[x_padj_col]] < padj_cutoff |
          .data[[y_padj_col]] < padj_cutoff
      )
    
  } else if (filter_mode == "both") {
    
    plot_df <- plot_df %>%
      filter(
        .data[[x_padj_col]] < padj_cutoff &
          .data[[y_padj_col]] < padj_cutoff
      )
  }
  
  plot_df <- plot_df %>%
    mutate(
      direction_class = case_when(
        .data[[x_col]] > 0 & .data[[y_col]] > 0 ~ "shared positive",
        .data[[x_col]] < 0 & .data[[y_col]] < 0 ~ "shared negative",
        .data[[x_col]] > 0 & .data[[y_col]] < 0 ~ "opposite direction",
        .data[[x_col]] < 0 & .data[[y_col]] > 0 ~ "opposite direction",
        TRUE ~ "unclear"
      )
    )
  
  if (nrow(plot_df) == 0) {
    message(
      "No valid rows for scatterplot: ",
      title,
      " using filter_mode = ",
      filter_mode,
      " and FDR cutoff = ",
      padj_cutoff
    )
    return(NULL)
  }
  
  axis_limit <- max(abs(c(plot_df[[x_col]], plot_df[[y_col]])), na.rm = TRUE)
  
  subtitle_text <- case_when(
    filter_mode == "none" ~ "No FDR filtering",
    filter_mode == "either" ~ paste0("Showing pathways with adjusted p-value < ", padj_cutoff, " in at least one comparison"),
    filter_mode == "both" ~ paste0("Showing pathways with adjusted p-value < ", padj_cutoff, " in both comparisons")
  )
  
  p <- ggplot(
    plot_df,
    aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      color = direction_class
    )
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(size = 3, alpha = 0.8) +
    coord_equal(
      xlim = c(-axis_limit, axis_limit),
      ylim = c(-axis_limit, axis_limit)
    ) +
    theme_bw() +
    labs(
      title = title,
      subtitle = subtitle_text,
      x = x_label,
      y = y_label,
      color = "Direction"
    )
  
  ggsave(
    output_file,
    p,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

############################################################
# 11. Mouse comparison plots
############################################################
# Heatmap column order:
#   KO, KD, het
############################################################

for (database_name in c("Hallmark", "KEGG", "Reactome", "GO_BP")) {
  
  plot_nes_heatmap(
    comparison_df = mouse_comparison_all,
    database_name = database_name,
    model_columns = c("NES_mouse_KO", "NES_mouse_KD", "NES_mouse_het"),
    output_file = file.path(
      comparison_figures_dir,
      "mouse_KO_KD_het",
      paste0("mouse_KO_KD_het_NES_heatmap_", database_name, ".png")
    ),
    title = paste0("Mouse KO vs KD vs het NES: ", database_name),
    top_n = 30,
    ranking_column = "delta_NES_KO_minus_het",
    clean_colnames = c("KO", "KD", "Het")
  )
  
  plot_pairwise_nes_scatter(
    comparison_df = mouse_comparison_all,
    database_name = database_name,
    x_col = "NES_mouse_het",
    y_col = "NES_mouse_KO",
    x_padj_col = "padj_mouse_het",
    y_padj_col = "padj_mouse_KO",
    x_label = "NES: mouse het vs WT",
    y_label = "NES: mouse KO vs WT",
    output_file = file.path(
      comparison_figures_dir,
      "mouse_KO_KD_het",
      paste0("mouse_KO_vs_het_NES_scatter_", database_name, "_FDR", scatter_fdr_cutoff, "_", scatter_filter_mode, ".png")
    ),
    title = paste0("Mouse KO vs het GSEA NES: ", database_name)
  )
  
  plot_pairwise_nes_scatter(
    comparison_df = mouse_comparison_all,
    database_name = database_name,
    x_col = "NES_mouse_KD",
    y_col = "NES_mouse_KO",
    x_padj_col = "padj_mouse_KD",
    y_padj_col = "padj_mouse_KO",
    x_label = "NES: mouse KD vs control",
    y_label = "NES: mouse KO vs WT",
    output_file = file.path(
      comparison_figures_dir,
      "mouse_KO_KD_het",
      paste0("mouse_KO_vs_KD_NES_scatter_", database_name, "_FDR", scatter_fdr_cutoff, "_", scatter_filter_mode, ".png")
    ),
    title = paste0("Mouse KO vs KD GSEA NES: ", database_name)
  )
  
  plot_pairwise_nes_scatter(
    comparison_df = mouse_comparison_all,
    database_name = database_name,
    x_col = "NES_mouse_het",
    y_col = "NES_mouse_KD",
    x_padj_col = "padj_mouse_het",
    y_padj_col = "padj_mouse_KD",
    x_label = "NES: mouse het vs WT",
    y_label = "NES: mouse KD vs control",
    output_file = file.path(
      comparison_figures_dir,
      "mouse_KO_KD_het",
      paste0("mouse_KD_vs_het_NES_scatter_", database_name, "_FDR", scatter_fdr_cutoff, "_", scatter_filter_mode, ".png")
    ),
    title = paste0("Mouse KD vs het GSEA NES: ", database_name)
  )
}

############################################################
# 12. Human KD vs mouse KD plots
############################################################
# Primary cross-species plots:
#   Hallmark
#   GO_BP
#
# Reactome and KEGG are still exported as CSVs.
# Plot them only if matched rows exist.
############################################################

primary_cross_species_databases <- c("Hallmark", "GO_BP")
exploratory_cross_species_databases <- c("Reactome", "KEGG")

for (database_name in primary_cross_species_databases) {
  
  plot_nes_heatmap(
    comparison_df = kd_cross_species_all,
    database_name = database_name,
    model_columns = c("NES_mouse_KD", "NES_human_KD"),
    output_file = file.path(
      comparison_figures_dir,
      "human_vs_mouse_KD",
      paste0("human_vs_mouse_KD_NES_heatmap_", database_name, ".png")
    ),
    title = paste0("Human KD vs mouse KD NES: ", database_name),
    top_n = 30,
    ranking_column = "delta_NES_human_minus_mouse",
    clean_colnames = c("Mouse KD", "Human KD")
  )
  
  plot_pairwise_nes_scatter(
    comparison_df = kd_cross_species_all,
    database_name = database_name,
    x_col = "NES_mouse_KD",
    y_col = "NES_human_KD",
    x_padj_col = "padj_mouse_KD",
    y_padj_col = "padj_human_KD",
    x_label = "NES: mouse KD vs control",
    y_label = "NES: human KD vs control",
    output_file = file.path(
      comparison_figures_dir,
      "human_vs_mouse_KD",
      paste0("human_vs_mouse_KD_NES_scatter_", database_name, "_FDR", scatter_fdr_cutoff, "_", scatter_filter_mode, ".png")
    ),
    title = paste0("Human KD vs mouse KD GSEA NES: ", database_name)
  )
}

for (database_name in exploratory_cross_species_databases) {
  
  matched_rows <- kd_cross_species_all %>%
    filter(
      database == database_name,
      !is.na(NES_mouse_KD),
      !is.na(NES_human_KD),
      is.finite(NES_mouse_KD),
      is.finite(NES_human_KD)
    )
  
  if (nrow(matched_rows) == 0) {
    message(
      "Skipping human-vs-mouse KD plots for ",
      database_name,
      ": no matched cross-species rows found."
    )
    next
  }
  
  plot_nes_heatmap(
    comparison_df = kd_cross_species_all,
    database_name = database_name,
    model_columns = c("NES_mouse_KD", "NES_human_KD"),
    output_file = file.path(
      comparison_figures_dir,
      "human_vs_mouse_KD",
      paste0("human_vs_mouse_KD_NES_heatmap_", database_name, ".png")
    ),
    title = paste0("Human KD vs mouse KD NES: ", database_name),
    top_n = 30,
    ranking_column = "delta_NES_human_minus_mouse",
    clean_colnames = c("Mouse KD", "Human KD")
  )
  
  plot_pairwise_nes_scatter(
    comparison_df = kd_cross_species_all,
    database_name = database_name,
    x_col = "NES_mouse_KD",
    y_col = "NES_human_KD",
    x_padj_col = "padj_mouse_KD",
    y_padj_col = "padj_human_KD",
    x_label = "NES: mouse KD vs control",
    y_label = "NES: human KD vs control",
    output_file = file.path(
      comparison_figures_dir,
      "human_vs_mouse_KD",
      paste0("human_vs_mouse_KD_NES_scatter_", database_name, "_FDR", scatter_fdr_cutoff, "_", scatter_filter_mode, ".png")
    ),
    title = paste0("Human KD vs mouse KD GSEA NES: ", database_name)
  )
}

############################################################
# 13. Cross-species concordance summary
############################################################
# This summarizes all matched Hallmark/GO_BP pathways with valid
# human and mouse NES values, and labels each pathway as significant
# in either/both species using scatter_fdr_cutoff.
############################################################

cross_species_concordance_summary <- kd_cross_species_all %>%
  filter(
    database %in% c("Hallmark", "GO_BP"),
    !is.na(NES_mouse_KD),
    !is.na(NES_human_KD),
    is.finite(NES_mouse_KD),
    is.finite(NES_human_KD)
  ) %>%
  mutate(
    significant_either = padj_mouse_KD < scatter_fdr_cutoff | padj_human_KD < scatter_fdr_cutoff,
    significant_both = padj_mouse_KD < scatter_fdr_cutoff & padj_human_KD < scatter_fdr_cutoff,
    direction_class = case_when(
      NES_mouse_KD > 0 & NES_human_KD > 0 ~ "shared positive",
      NES_mouse_KD < 0 & NES_human_KD < 0 ~ "shared negative",
      NES_mouse_KD > 0 & NES_human_KD < 0 ~ "opposite direction",
      NES_mouse_KD < 0 & NES_human_KD > 0 ~ "opposite direction",
      TRUE ~ "unclear"
    )
  ) %>%
  group_by(database, significant_either, significant_both, direction_class) %>%
  summarise(
    n_pathways = n(),
    mean_delta_NES = mean(delta_NES_human_minus_mouse, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  cross_species_concordance_summary,
  file.path(
    comparison_results_dir,
    "human_vs_mouse_KD",
    "human_vs_mouse_KD_concordance_summary.csv"
  ),
  row.names = FALSE
)

############################################################
# 14. Helper: export filtered pathway lists used for pairwise scatterplots
############################################################

make_filtered_pairwise_pathway_list <- function(
    comparison_df,
    database_name,
    x_model,
    y_model,
    x_label,
    y_label,
    padj_cutoff = scatter_fdr_cutoff,
    filter_mode = scatter_filter_mode
) {
  
  filter_mode <- match.arg(filter_mode, choices = c("none", "either", "both"))
  
  x_nes_col <- paste0("NES_", x_model)
  y_nes_col <- paste0("NES_", y_model)
  x_padj_col <- paste0("padj_", x_model)
  y_padj_col <- paste0("padj_", y_model)
  x_le_col <- paste0("leading_edge_", x_model)
  y_le_col <- paste0("leading_edge_", y_model)
  
  out_df <- comparison_df %>%
    filter(
      database == database_name,
      !is.na(.data[[x_nes_col]]),
      !is.na(.data[[y_nes_col]]),
      is.finite(.data[[x_nes_col]]),
      is.finite(.data[[y_nes_col]]),
      !is.na(.data[[x_padj_col]]),
      !is.na(.data[[y_padj_col]]),
      is.finite(.data[[x_padj_col]]),
      is.finite(.data[[y_padj_col]])
    )
  
  if (filter_mode == "either") {
    out_df <- out_df %>%
      filter(
        .data[[x_padj_col]] < padj_cutoff |
          .data[[y_padj_col]] < padj_cutoff
      )
  } else if (filter_mode == "both") {
    out_df <- out_df %>%
      filter(
        .data[[x_padj_col]] < padj_cutoff &
          .data[[y_padj_col]] < padj_cutoff
      )
  }
  
  out_df <- out_df %>%
    mutate(
      comparison = paste0(y_label, " vs ", x_label),
      filter_mode_used = filter_mode,
      fdr_cutoff_used = padj_cutoff,
      delta_NES = .data[[y_nes_col]] - .data[[x_nes_col]],
      abs_delta_NES = abs(delta_NES),
      direction_class = case_when(
        .data[[x_nes_col]] > 0 & .data[[y_nes_col]] > 0 ~ "shared positive",
        .data[[x_nes_col]] < 0 & .data[[y_nes_col]] < 0 ~ "shared negative",
        .data[[x_nes_col]] > 0 & .data[[y_nes_col]] < 0 ~ "opposite direction",
        .data[[x_nes_col]] < 0 & .data[[y_nes_col]] > 0 ~ "opposite direction",
        TRUE ~ "unclear"
      )
    ) %>%
    select(
      comparison,
      database,
      pathway_key,
      pathway_name,
      filter_mode_used,
      fdr_cutoff_used,
      direction_class,
      delta_NES,
      abs_delta_NES,
      NES_x = all_of(x_nes_col),
      NES_y = all_of(y_nes_col),
      padj_x = all_of(x_padj_col),
      padj_y = all_of(y_padj_col),
      leading_edge_x = any_of(x_le_col),
      leading_edge_y = any_of(y_le_col),
      everything()
    ) %>%
    arrange(desc(abs_delta_NES))
  
  return(out_df)
}

############################################################
# 15. Export filtered pathway lists for all scatterplots
############################################################
# Files are saved into:
#
#   Model_Comparisons/results/mouse_KO_KD_het/filtered_scatter_pathway_lists/
#   Model_Comparisons/results/human_vs_mouse_KD/filtered_scatter_pathway_lists/
#
# using the same FDR rule as the scatterplots.
############################################################

mouse_scatter_list_dir <- file.path(
  comparison_results_dir,
  "mouse_KO_KD_het",
  "filtered_scatter_pathway_lists"
)

human_mouse_scatter_list_dir <- file.path(
  comparison_results_dir,
  "human_vs_mouse_KD",
  "filtered_scatter_pathway_lists"
)

dir.create(mouse_scatter_list_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(human_mouse_scatter_list_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
# 15A. Mouse KO/KD/het filtered scatter pathway lists
############################################################

mouse_filtered_scatter_lists <- list()

for (database_name in c("Hallmark", "KEGG", "Reactome", "GO_BP")) {
  
  mouse_filtered_scatter_lists[[paste0("mouse_KO_vs_het_", database_name)]] <-
    make_filtered_pairwise_pathway_list(
      comparison_df = mouse_comparison_all,
      database_name = database_name,
      x_model = "mouse_het",
      y_model = "mouse_KO",
      x_label = "mouse het vs WT",
      y_label = "mouse KO vs WT"
    )
  
  mouse_filtered_scatter_lists[[paste0("mouse_KO_vs_KD_", database_name)]] <-
    make_filtered_pairwise_pathway_list(
      comparison_df = mouse_comparison_all,
      database_name = database_name,
      x_model = "mouse_KD",
      y_model = "mouse_KO",
      x_label = "mouse KD vs control",
      y_label = "mouse KO vs WT"
    )
  
  mouse_filtered_scatter_lists[[paste0("mouse_KD_vs_het_", database_name)]] <-
    make_filtered_pairwise_pathway_list(
      comparison_df = mouse_comparison_all,
      database_name = database_name,
      x_model = "mouse_het",
      y_model = "mouse_KD",
      x_label = "mouse het vs WT",
      y_label = "mouse KD vs control"
    )
}

for (list_name in names(mouse_filtered_scatter_lists)) {
  
  write.csv(
    mouse_filtered_scatter_lists[[list_name]],
    file.path(
      mouse_scatter_list_dir,
      paste0(list_name, "_FDR", scatter_fdr_cutoff, "_", scatter_filter_mode, ".csv")
    ),
    row.names = FALSE
  )
}

mouse_filtered_scatter_lists_all <- bind_rows(
  mouse_filtered_scatter_lists,
  .id = "scatter_comparison"
)

write.csv(
  mouse_filtered_scatter_lists_all,
  file.path(
    mouse_scatter_list_dir,
    paste0(
      "ALL_mouse_KO_KD_het_filtered_scatter_pathway_lists_FDR",
      scatter_fdr_cutoff,
      "_",
      scatter_filter_mode,
      ".csv"
    )
  ),
  row.names = FALSE
)

############################################################
# 15B. Human KD vs mouse KD filtered scatter pathway lists
############################################################

human_mouse_filtered_scatter_lists <- list()

for (database_name in c("Hallmark", "GO_BP", "Reactome", "KEGG")) {
  
  human_mouse_filtered_scatter_lists[[paste0("human_vs_mouse_KD_", database_name)]] <-
    make_filtered_pairwise_pathway_list(
      comparison_df = kd_cross_species_all,
      database_name = database_name,
      x_model = "mouse_KD",
      y_model = "human_KD",
      x_label = "mouse KD vs control",
      y_label = "human KD vs control"
    )
}

for (list_name in names(human_mouse_filtered_scatter_lists)) {
  
  write.csv(
    human_mouse_filtered_scatter_lists[[list_name]],
    file.path(
      human_mouse_scatter_list_dir,
      paste0(list_name, "_FDR", scatter_fdr_cutoff, "_", scatter_filter_mode, ".csv")
    ),
    row.names = FALSE
  )
}

human_mouse_filtered_scatter_lists_all <- bind_rows(
  human_mouse_filtered_scatter_lists,
  .id = "scatter_comparison"
)

write.csv(
  human_mouse_filtered_scatter_lists_all,
  file.path(
    human_mouse_scatter_list_dir,
    paste0(
      "ALL_human_vs_mouse_KD_filtered_scatter_pathway_lists_FDR",
      scatter_fdr_cutoff,
      "_",
      scatter_filter_mode,
      ".csv"
    )
  ),
  row.names = FALSE
)

############################################################
# 16. Save session info
############################################################

writeLines(
  capture.output(sessionInfo()),
  file.path(comparison_results_dir, "sessionInfo_model_comparisons.txt")
)
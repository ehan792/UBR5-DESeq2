################################################################################
# 00_config.R — UBR5 / MPNST RNA-seq pipeline configuration
#
# Purpose:
#   Shared paths, experiment definitions, package loading, gene-set helpers,
#   plotting theme, and output path helpers for scripts 01–04.
#
# Design choice:
#   This file defines config globally, but downstream .rds objects do NOT save
#   the full config or OrgDb objects. Later scripts source 00_config.R again.
#   This prevents the C-stack / heavy serialization issue.
################################################################################

# ── 0. Project paths ─────────────────────────────────────────────────────────
PROJECT_ROOT <- Sys.getenv("UBR5_PROJECT_ROOT", unset = NA_character_)
if (is.na(PROJECT_ROOT) || PROJECT_ROOT == "") {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (dir.exists(file.path(wd, "Data"))) {
    PROJECT_ROOT <- wd
  } else if (dir.exists(file.path(dirname(wd), "Data"))) {
    PROJECT_ROOT <- dirname(wd)
  } else {
    PROJECT_ROOT <- wd
  }
}
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = FALSE)
DATA_DIR <- file.path(PROJECT_ROOT, "Data")
OUT <- file.path(PROJECT_ROOT, "results")

COUNT_FILES <- list(
  CRISPR  = file.path(DATA_DIR, "MouseCRISPRCounts.csv"),
  mouseKD = file.path(DATA_DIR, "MouseKDCounts.xlsx"),
  humanKD = file.path(DATA_DIR, "HumanKDCounts.xlsx")
)

# Optional user-curated gene panel. Accepted columns are either:
#   Gene, Class
# or
#   gene_symbol, class
BIOLOGY_PANEL_CSV <- file.path(DATA_DIR, "UBR5_biology_gene_table_ALL.csv")

GLOBAL_SEED <- 42
set.seed(GLOBAL_SEED)

# ── 1. Packages ──────────────────────────────────────────────────────────────
# CRAN and Bioconductor dependencies are listed explicitly (rather than
# discovered implicitly) so the required environment is documented in one
# place. For fully locked, publication-grade reproducibility, pair this with
# an `renv` lockfile (`renv::init()` / `renv::snapshot()` at the project root);
# the install-if-missing logic below is a lightweight fallback for that.
cran_pkgs <- c(
  "tidyverse", "readxl", "ggrepel", "pheatmap", "RColorBrewer",
  "scales", "matrixStats", "viridisLite"
)
bioc_pkgs <- c(
  "DESeq2", "apeglm", "org.Mm.eg.db", "org.Hs.eg.db", "clusterProfiler",
  "fgsea", "msigdbr", "GSVA", "AnnotationDbi"
)
optional_pkgs <- c("limma")

install_if_missing <- function(pkgs, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing missing package: ", p)
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(p, update = FALSE, ask = FALSE)
      } else {
        install.packages(p)
      }
    }
  }
}
install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)
for (p in optional_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    try(install_if_missing(p, bioc = TRUE), silent = TRUE)
  }
}

still_missing <- Filter(function(p) !requireNamespace(p, quietly = TRUE), c(cran_pkgs, bioc_pkgs))
if (length(still_missing)) {
  stop(
    "Required package(s) could not be installed automatically: ",
    paste(still_missing, collapse = ", "),
    ". Install manually, e.g. BiocManager::install(c(",
    paste(sprintf("\"%s\"", still_missing), collapse = ", "),
    ")), then re-source 00_config.R."
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(matrixStats)
  library(DESeq2)
  library(apeglm)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(fgsea)
  library(msigdbr)
  library(GSVA)
  library(AnnotationDbi)
})

# ── 1b. Namespace safety ─────────────────────────────────────────────────────
# AnnotationDbi (and, via DESeq2, BiocGenerics) re-export many base/dplyr verb
# names as S4 generics — e.g. AnnotationDbi::select() vs dplyr::select(), or
# BiocGenerics::setdiff()/intersect()/table() vs the base versions — which
# otherwise fails with errors like "unable to find an inherited method for
# function 'select' for signature x = 'data.frame'".
#
# {conflicted} was tried here and rejected: it requires an explicit
# preference for *every* ambiguous name on the search path, and BiocGenerics
# alone re-exports several dozen base functions, so it surfaces one masking
# error at a time deep into a run rather than up front. Explicit dplyr
# aliases in the global session are more predictable for this particular
# tidyverse + Bioconductor combination. Sites that need the non-dplyr version
# (e.g. AnnotationDbi::select()) call it via an explicit `::`.
select     <- dplyr::select
filter     <- dplyr::filter
rename     <- dplyr::rename
mutate     <- dplyr::mutate
summarise  <- dplyr::summarise
summarize  <- dplyr::summarize
arrange    <- dplyr::arrange
distinct   <- dplyr::distinct
group_by   <- dplyr::group_by
ungroup    <- dplyr::ungroup
left_join  <- dplyr::left_join
inner_join <- dplyr::inner_join
full_join  <- dplyr::full_join
bind_rows  <- dplyr::bind_rows
case_when  <- dplyr::case_when
across     <- dplyr::across
everything <- dplyr::everything
any_of     <- dplyr::any_of
all_of     <- dplyr::all_of
n_distinct <- dplyr::n_distinct

# ── 2. Output helpers ────────────────────────────────────────────────────────
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

p_data <- function(...) file.path(ensure_dir(file.path(OUT, "data")), ...)
p_log  <- function(...) file.path(ensure_dir(file.path(OUT, "logs")), ...)

# No contrast subfolders. Put contrast names in filenames for clarity.
p_tab <- function(section, experiment = NULL, file) {
  parts <- c(OUT, "tables", section, experiment)
  dir <- do.call(file.path, as.list(parts[!is.na(parts) & nzchar(parts)]))
  file.path(ensure_dir(dir), file)
}
p_fig <- function(section, experiment = NULL, file) {
  parts <- c(OUT, "figures", section, experiment)
  dir <- do.call(file.path, as.list(parts[!is.na(parts) & nzchar(parts)]))
  file.path(ensure_dir(dir), file)
}

# Build "<experiment>__<suffix>" filenames, the convention used throughout
# scripts 01-04. `id` may be an experiment config list (uses id$exp_id) or a
# bare experiment/group id string (e.g. "pooledMouse", "mouseDose").
out_name <- function(id, suffix) {
  exp_id <- if (is.list(id)) id$exp_id else id
  paste0(exp_id, "__", suffix)
}

# Thin wrappers around write.csv()/ggsave() that route through p_tab()/p_fig()
# so every table/figure export in 01-04 shares one place for path building
# and CSV/plot defaults (e.g. row.names = FALSE, dpi = 300).
write_tab <- function(df, section, experiment = NULL, file, ...) {
  write.csv(df, p_tab(section, experiment, file = file), row.names = FALSE, ...)
  invisible(df)
}
save_fig <- function(plot, section, experiment = NULL, file, width, height, dpi = 300, ...) {
  ggsave(p_fig(section, experiment, file = file), plot = plot, width = width, height = height, dpi = dpi, ...)
  invisible(plot)
}
# For base-graphics/pheatmap output, which must be drawn inside an open
# png() device rather than returned as a ggplot object. `draw` is a
# zero-argument function that performs the plotting call(s).
save_png <- function(section, experiment = NULL, file, width, height, res, draw) {
  png(p_fig(section, experiment, file = file), width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  draw()
}

invisible(lapply(c(
  file.path(OUT, "data"),
  file.path(OUT, "logs"),
  file.path(OUT, "tables", "01_deseq2"),
  file.path(OUT, "tables", "02_qc"),
  file.path(OUT, "tables", "03_gsea"),
  file.path(OUT, "tables", "04_biology_concordance"),
  file.path(OUT, "figures", "02_qc"),
  file.path(OUT, "figures", "03_gsea"),
  file.path(OUT, "figures", "04_biology_concordance")
), ensure_dir))

# ── 3. Experiment metadata ──────────────────────────────────────────────────
# Folder names intentionally match requested organization: CRISPR, mouseKD, humanKD.
EXPERIMENTS <- list(
  CRISPR = list(
    exp_id = "CRISPR",
    label = "Mouse CRISPR Ubr5 series",
    file = COUNT_FILES$CRISPR,
    filetype = "csv",
    condition_map = list(WT = "\\.WT_", Het = "pos\\.neg", KO = "neg\\.neg"),
    ref_level = "WT",
    contrasts = list(
      list(name = "Het_vs_WT", treat = "Het", ctrl = "WT"),
      list(name = "KO_vs_WT",  treat = "KO",  ctrl = "WT")
    ),
    species = "Mus musculus",
    orgdb_pkg = "org.Mm.eg.db",
    msig_species = "Mus musculus",
    kegg_org = "mmu",
    ubr5_symbol = "Ubr5"
  ),
  mouseKD = list(
    exp_id = "mouseKD",
    label = "Mouse Ubr5 knockdown",
    file = COUNT_FILES$mouseKD,
    filetype = "xlsx",
    condition_map = list(Control = "jw23\\.3", KD = "shrubr5"),
    ref_level = "Control",
    contrasts = list(list(name = "KD_vs_Control", treat = "KD", ctrl = "Control")),
    species = "Mus musculus",
    orgdb_pkg = "org.Mm.eg.db",
    msig_species = "Mus musculus",
    kegg_org = "mmu",
    ubr5_symbol = "Ubr5"
  ),
  humanKD = list(
    exp_id = "humanKD",
    label = "Human UBR5 knockdown",
    file = COUNT_FILES$humanKD,
    filetype = "xlsx",
    condition_map = list(Control = "jh_2_002", KD = "shrubr5"),
    ref_level = "Control",
    contrasts = list(list(name = "KD_vs_Control", treat = "KD", ctrl = "Control")),
    species = "Homo sapiens",
    orgdb_pkg = "org.Hs.eg.db",
    msig_species = "Homo sapiens",
    kegg_org = "hsa",
    ubr5_symbol = "UBR5"
  )
)

# Priority Hallmark sets for MPNST/UBR5 biology summaries.
MPNST_PRIORITY_HALLMARK <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_APOPTOSIS"
)

COND_COLOURS <- c(
  Control = "#2166AC", WT = "#2166AC",
  KD = "#F4A582", Het = "#D6604D", KO = "#A50026"
)
CONTRAST_COLOURS <- c(
  KD_vs_Control = "#F4A582", Het_vs_WT = "#D6604D", KO_vs_WT = "#A50026"
)

GSEA_PADJ_LOW_COLOUR <- "#001F4D"  # dark navy for most significant padj values
GSEA_PADJ_HIGH_COLOUR <- "#F4A6A6" # light red for larger padj values in single-contrast GSEA plots

# Colours for plots that compare multiple mouse perturbations/experiments.
# These encode the comparison identity rather than statistical significance.
MOUSE_DOSE_COLOURS <- c(
  mouseKD = "#2166AC",
  `CRISPR Het` = "#F4A582",
  `CRISPR KO` = "#B2182B"
)

# ── 4. General helpers ──────────────────────────────────────────────────────
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace% theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "grey96", colour = "grey45", linewidth = 0.4),
    plot.title = element_text(face = "bold"),
    legend.key = element_blank()
  )
}
theme_set(theme_pub())

clean_sample_label <- function(x) {
  x %>%
    str_replace("^sample\\.", "") %>%
    str_replace_all("pos\\.neg", "Het") %>%
    str_replace_all("neg\\.neg", "KO") %>%
    str_replace_all("WT_", "WT ") %>%
    str_replace_all("shrubr5", "KD") %>%
    str_replace_all("jw23\\.3", "Control") %>%
    str_replace_all("jh_2_002", "Control") %>%
    str_replace_all("_", " ")
}


# Explicit sample ordering for figures/tables where sample order is displayed.
# CRISPR order requested: KO_1–3, Het_1–3, WT_1–3. KD order is kept as-is.
sample_order_for_experiment <- function(samples, exp_id) {
  samples <- as.character(samples)
  if (identical(exp_id, "CRISPR")) {
    rank_one <- function(x) {
      # sample.neg.neg_1/2/3 = KO; sample.pos.neg_1/2/3 = Het; sample.WT_1/2/3 = WT
      group_rank <- dplyr::case_when(
        stringr::str_detect(x, "neg\\.neg") ~ 1L,
        stringr::str_detect(x, "pos\\.neg") ~ 2L,
        stringr::str_detect(x, "WT_") ~ 3L,
        TRUE ~ 9L
      )
      rep_num <- suppressWarnings(as.integer(stringr::str_match(x, "_(\\d+)(?:$|[^0-9])")[, 2]))
      ifelse(is.na(rep_num), 99L, rep_num) + group_rank * 100L
    }
    samples[order(vapply(samples, rank_one, numeric(1)), samples)]
  } else {
    samples
  }
}

write_session <- function(script_id) {
  writeLines(capture.output(sessionInfo()), p_log(paste0("sessionInfo_", script_id, ".txt")))
}

get_orgdb <- function(orgdb_pkg) {
  switch(orgdb_pkg,
         org.Mm.eg.db = org.Mm.eg.db,
         org.Hs.eg.db = org.Hs.eg.db,
         stop("Unsupported OrgDb package: ", orgdb_pkg))
}

filter_by_group <- function(mat, sample_info, min_count = 10, min_samp = 3) {
  conds <- levels(droplevels(sample_info$condition))
  keep_list <- lapply(conds, function(cd) {
    idx <- rownames(sample_info)[sample_info$condition == cd]
    needed <- min(min_samp, length(idx))
    rowSums(mat[, idx, drop = FALSE] >= min_count) >= needed
  })
  Reduce(`|`, keep_list)
}

# Map Ensembl rows to Entrez rows and average duplicated Entrez IDs. Useful for GSVA/ssGSEA.
ens_to_entrez_mat <- function(mat, orgdb_pkg) {
  orgdb <- get_orgdb(orgdb_pkg)
  map <- AnnotationDbi::select(orgdb, keys = rownames(mat), keytype = "ENSEMBL", columns = "ENTREZID") %>%
    as_tibble() %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENSEMBL, ENTREZID)
  df <- as.data.frame(mat) %>%
    rownames_to_column("ENSEMBL") %>%
    inner_join(map, by = "ENSEMBL") %>%
    dplyr::select(-ENSEMBL) %>%
    group_by(ENTREZID) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")
  out <- as.matrix(df[, -1, drop = FALSE])
  rownames(out) <- df$ENTREZID
  out
}

entrez_to_symbol <- function(ids, orgdb_pkg) {
  ids <- unique(as.character(ids[!is.na(ids) & ids != ""]))
  if (!length(ids)) return(setNames(character(0), character(0)))
  orgdb <- get_orgdb(orgdb_pkg)
  sym <- AnnotationDbi::mapIds(orgdb, keys = ids, keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")
  sym[is.na(sym)] <- names(sym)[is.na(sym)]
  sym
}

# Robust GSVA wrapper for old and new GSVA APIs.
safe_gsva <- function(expr, gene_sets, method = c("ssgsea", "gsva"), min_size = 5, max_size = 500) {
  method <- match.arg(method)
  sizes <- vapply(gene_sets, length, integer(1))
  gene_sets <- gene_sets[sizes >= min_size & sizes <= max_size]
  if (!length(gene_sets)) stop("No gene sets remain after size filtering.")
  if (method == "ssgsea") {
    tryCatch(
      GSVA::gsva(GSVA::ssgseaParam(expr, gene_sets, minSize = min_size, maxSize = max_size), verbose = FALSE),
      error = function(e) GSVA::gsva(expr, gene_sets, method = "ssgsea", min.sz = min_size, max.sz = max_size, kcdf = "Gaussian", verbose = FALSE)
    )
  } else {
    tryCatch(
      GSVA::gsva(GSVA::gsvaParam(expr, gene_sets, minSize = min_size, maxSize = max_size), verbose = FALSE),
      error = function(e) GSVA::gsva(expr, gene_sets, method = "gsva", min.sz = min_size, max.sz = max_size, kcdf = "Gaussian", verbose = FALSE)
    )
  }
}

# ── 5. MSigDB / GSEA gene-set helpers ────────────────────────────────────────
msig_gene_col <- function(df) {
  cols <- c("ncbi_gene", "entrez_gene", "entrezgene", "gene_symbol", "human_gene_symbol")
  hit <- cols[cols %in% colnames(df)][1]
  if (is.na(hit)) stop("No supported gene-ID column found in msigdbr output.")
  hit
}

fetch_msigdb <- function(species, collection, subcollection = NULL) {
  # Handles msigdbr API changes: collection/subcollection vs category/subcategory,
  # and mouse-specific collection names in newer MSigDB releases.
  is_mouse <- identical(species, "Mus musculus")
  coll_try <- collection
  if (is_mouse) {
    coll_try <- switch(collection, H = "MH", C2 = "M2", C5 = "M5", collection)
  }
  attempts <- list(
    list(db_species = if (is_mouse) "MM" else "HS", species = species, collection = coll_try, subcollection = subcollection),
    list(species = species, collection = collection, subcollection = subcollection),
    list(species = species, category = collection, subcategory = subcollection),
    list(species = species, collection = collection),
    list(species = species, category = collection)
  )
  for (args in attempts) {
    args <- args[!vapply(args, is.null, logical(1))]
    out <- tryCatch(suppressWarnings(do.call(msigdbr::msigdbr, args)), error = function(e) NULL)
    if (!is.null(out) && nrow(out) > 0) return(as_tibble(out))
  }
  stop("Could not fetch MSigDB collection ", collection, " / ", subcollection, " for ", species)
}

msig_to_pathways <- function(msig_df) {
  gene_col <- msig_gene_col(msig_df)
  split(as.character(msig_df[[gene_col]]), msig_df$gs_name) %>%
    lapply(function(x) unique(na.omit(x)))
}

get_gsea_collections <- function(cfg) {
  # All returned gene sets use Entrez IDs, so they match the ranked vector.
  collections <- list()
  collections$Hallmark <- msig_to_pathways(fetch_msigdb(cfg$msig_species, "H"))
  collections$GO_BP    <- msig_to_pathways(fetch_msigdb(cfg$msig_species, "C5", "GO:BP"))
  collections$Reactome <- msig_to_pathways(fetch_msigdb(cfg$msig_species, "C2", "CP:REACTOME"))
  collections$KEGG     <- msig_to_pathways(fetch_msigdb(cfg$msig_species, "C2", "CP:KEGG"))
  collections
}

message("00_config.R loaded")
message("  PROJECT_ROOT: ", PROJECT_ROOT)
message("  DATA_DIR:     ", DATA_DIR)
message("  OUT:          ", OUT)

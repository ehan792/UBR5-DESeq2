# UBR5 / MPNST RNA-seq pipeline

Four analysis scripts, run in order, covering DESeq2 differential expression, cell-state QC, GSEA, and biology-panel/cross-experiment concordance for the CRISPR, mouseKD, and humanKD experiments.

## Required data files

Put these in `PROJECT_ROOT/Data/`:

- `MouseCRISPRCounts.csv`
- `MouseKDCounts.xlsx`
- `HumanKDCounts.xlsx`

Optional biology panel:

- `UBR5_biology_gene_table_ALL.csv`
  - accepted columns: `Class,Gene` or `class,gene_symbol`

## Run

Open RStudio in the script folder or project folder, then run:

```r
Sys.setenv(UBR5_PROJECT_ROOT = "/path/to/UBR5-DESeq2") # optional but recommended
source("run_all.R")
```

## Scripts

### 00_config.R

Shared configuration only: paths, experiment metadata, packages, output helpers, gene-set helpers, and plotting theme.

The DESeq2 saved objects do **not** store full config or OrgDb objects. This prevents the previous C-stack / heavy serialization issue.

### 01_import_deseq2.R

Preprocesses counts, runs DESeq2, exports normalized counts, UBR5 percent-reference tables, Wald tables, shrunken LFC tables, and DE summaries.

Outputs go to:

```text
results/tables/01_deseq2/CRISPR/
results/tables/01_deseq2/mouseKD/
results/tables/01_deseq2/humanKD/
results/data/deseq_all.rds
```

Each contrast also gets a `<experiment>__<contrast>__DESeq2_Wald_by_symbol.csv` — the same Wald statistics as `..._DESeq2_Wald.csv`, but indexed by gene symbol and deduplicated to one row per symbol (largest `|stat|` wins). This is the handoff format for the planned Python TF-activity-inference step (`decoupler` + CollecTRI), which resolves regulons by gene symbol rather than Ensembl/Entrez ID. In pandas: `pd.read_csv(path).set_index("gene_symbol")["stat"]`.

### 02_cell_state_qc.R

Runs QC assays:

- PCA
- sample distance heatmaps
- pooled mouse PCA
- Hallmark ssGSEA
- dispersion estimates

Outputs go to:

```text
results/tables/02_qc/<experiment>/
results/figures/02_qc/<experiment>/
```

### 03_gsea_and_plots.R

Runs GSEA for:

- Hallmark
- GO Biological Process
- Reactome
- KEGG

Also exports:

- volcano plots
- GSEA dotplots ordered by `padj`
- dotplots filtered to `padj <= 0.10`
- important gene-set tables with leading-edge gene symbols
- mouse dose-priority MPNST Hallmark NES plot
- mouse Hallmark Kendall NES ordering
- leading-edge Jaccard comparisons

Outputs go to:

```text
results/tables/03_gsea/<experiment>/
results/figures/03_gsea/<experiment>/
results/tables/03_gsea/mouseDose/
results/figures/03_gsea/mouseDose/
```

### 04_biology_panel_concordance.R

Creates/uses the biology gene panel and exports:

- biology panel heatmaps
- panel DESeq2 tables
- DEG lists
- pairwise DEG overlap Fisher/Jaccard tables
- LFC concordance plots

Outputs go to:

```text
results/tables/04_biology_concordance/<experiment>/
results/figures/04_biology_concordance/<experiment>/
results/tables/04_biology_concordance/crossExperiment/
results/figures/04_biology_concordance/crossExperiment/
```

## Folder organization

Contrasts are not separated into extra subfolders. Instead, contrast names are embedded directly in filenames, for example:

```text
results/tables/03_gsea/CRISPR/CRISPR__KO_vs_WT__Hallmark__GSEA.csv
results/figures/03_gsea/mouseKD/mouseKD__KD_vs_Control__Hallmark__dotplot_padj0.10.png
```

This keeps the hierarchy shallow:

```text
results/
  tables/
    01_deseq2/
      CRISPR/
      mouseKD/
      humanKD/
    02_qc/
    03_gsea/
    04_biology_concordance/
  figures/
    02_qc/
    03_gsea/
    04_biology_concordance/
  data/
  logs/
```

## Code organization

`00_config.R` defines a few small IO helpers used throughout `01`-`04` to avoid repeating path/filename construction at each of the ~60 table and figure exports:

- `out_name(cfg_or_id, suffix)` — builds the `<experiment>__<suffix>` filename convention (`cfg_or_id` can be an experiment config list or a bare id like `"pooledMouse"`).
- `write_tab(df, section, experiment, file)` — `write.csv()` through `p_tab()`, always `row.names = FALSE`.
- `save_fig(plot, section, experiment, file, width, height)` — `ggsave()` through `p_fig()`, always `dpi = 300`.
- `save_png(section, experiment, file, width, height, res, draw)` — opens a `png()` device, runs `draw()` (a base-graphics/`pheatmap()` call), and always closes the device, for the few plots that aren't `ggplot` objects.

## Package management

Packages are loaded explicitly (CRAN vs. Bioconductor vs. optional) in `00_config.R`, with a fallback installer for anything missing. After install attempts, the script verifies every required package actually loaded and stops with an explicit list if not, rather than failing later with a cryptic `library()` error.

Namespace collisions (e.g. `AnnotationDbi::select()`/`BiocGenerics::setdiff()` vs. the `dplyr`/base versions) are resolved by explicitly re-binding the dplyr verbs in the global session (`00_config.R`, "Namespace safety"), so unqualified calls like `filter()` or `mutate()` always resolve to `dplyr`. The [`conflicted`](https://conflicted.r-lib.org/) package was evaluated as a more idiomatic alternative but rejected: `BiocGenerics` (pulled in by `DESeq2`) re-exports several dozen base functions as S4 generics (`setdiff`, `intersect`, `table`, `order`, `unique`, `Reduce`, ...), so `conflicted` surfaces one masking error at a time, deep into a pipeline run, rather than up front. Sites that need a non-dplyr version of a masked verb call it via an explicit `::` (e.g. `AnnotationDbi::select()`).

For fully locked, publication-grade reproducibility beyond this install-if-missing fallback, pair the project with an [`renv`](https://rstudio.github.io/renv/) lockfile (`renv::init()` / `renv::snapshot()` at the project root).

## Implementation notes

- `deseq_all.rds` and the per-experiment `<experiment>__deseq_bundle.rds` files never store the full `cfg` list or OrgDb objects, only the fields scripts 02-04 need. Serializing Bioconductor annotation databases into an `.rds` triggers a C-stack error; this avoids it.
- Sample metadata is attached to expression/GSEA/QC tables via explicit joins on `sample` or `ensembl_gene_id`, not row order, to avoid duplicate/misaligned `sample` columns.
- Data frames reloaded from `.rds` can carry stale row names, so `tibble::remove_rownames()` runs before every `column_to_rownames()` in the QC and biology-panel heatmap annotation code (otherwise: `.data must be a data frame without row names`).
- `pheatmap()` is always called with `angle_col = "45"` (a string) — the numeric `45` is silently mishandled by pheatmap.
- Bare `count()` is avoided everywhere in favor of explicit `dplyr::n()`/`n_distinct()`, since `count()` is masked by more than one loaded package (see "Package management" above).
- `run_all.R` wraps each script in `withCallingHandlers()` so warnings are logged without an `invokeRestart("muffleWarning")` crash, and stops the pipeline (rather than printing a false success message) if any script errors.
- CRISPR biology-panel heatmaps (`04_biology_panel_concordance.R`) show samples in a fixed order — KO_1-3, Het_1-3, WT_1-3 — via `sample_order_for_experiment()`; KD experiments keep their original column order.
- Figure colors encode either statistical significance (GSEA dotplots and the mouse-dose priority plot: dark navy = low padj, light red = high padj) or group/comparison identity (PCA and ssGSEA: `COND_COLOURS`; mouse-dose comparison bars: `MOUSE_DOSE_COLOURS`) — never both in the same plot.
- `01_import_deseq2.R` additionally exports combined UBR5 normalized-expression ratio-to-reference tables across all three experiments:
  - `results/tables/01_deseq2/all_experiments__UBR5_normalized_expression_ratio_to_reference_per_sample.csv`
  - `results/tables/01_deseq2/all_experiments__UBR5_normalized_expression_ratio_to_reference_by_phenotype.csv`


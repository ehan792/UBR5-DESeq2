# UBR5 / MPNST RNA-seq pipeline

Five analysis scripts, run in order: four in R (00-04) covering DESeq2 differential expression, cell-state QC, GSEA, and biology-panel/cross-experiment concordance, then one in Python (05) covering transcription-factor activity inference — for the CRISPR, mouseKD, and humanKD experiments.

## Required data files

Put these in `PROJECT_ROOT/Data/`:

- `MouseCRISPRCounts.csv`
- `MouseKDCounts.xlsx`
- `HumanKDCounts.xlsx`

Optional biology panel:

- `UBR5_biology_gene_table_ALL.csv`
  - accepted columns: `Class,Gene` or `class,gene_symbol`

## Run

### R stages (00-04)

Open RStudio in the script folder or project folder, then run:

```r
Sys.setenv(UBR5_PROJECT_ROOT = "/path/to/UBR5-DESeq2") # optional but recommended
source("run_all.R")
```

### Python stage (05)

Requires 01-04 to have already run (it reads their `*_DESeq2_Wald_by_symbol.csv` outputs). Uses a virtual environment, isolated from the R setup:

```bash
cd Script
python3 -m venv ../.venv
source ../.venv/bin/activate   # Windows: ..\.venv\Scripts\activate
pip install -r requirements.txt
python 05_tf_inference.py
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

### 05_tf_inference.py

Infers transcription factor (TF) activity per experiment/contrast using [decoupler](https://decoupler.readthedocs.io)'s Univariate Linear Model (ULM) against the [CollecTRI](https://github.com/saezlab/CollecTRI) TF-target regulon network (fetched via OmniPath, cached locally per organism). Input is each contrast's `DESeq2_Wald_by_symbol.csv` from `01_import_deseq2.R` — no additional normalization needed, since ULM is designed to run directly on contrast-level statistics (Wald stat, logFC, etc.), and the Wald statistic (rather than raw log2FoldChange) is used specifically because it's variance-aware.

Outputs go to:

```text
results/tables/05_tf_inference/<experiment>/
  <experiment>__TF_activity_ulm_wide.csv   (rows=contrast, columns=TF)
  <experiment>__TF_activity_ulm_long.csv   (one row per experiment/contrast/TF, sorted by padj, with n_targets)
results/figures/05_tf_inference/<experiment>/
  <experiment>__<contrast>__top_TF_barplot.png
results/tables/05_tf_inference/mouseDose/
  mouseDose__all_experiments_TF_activity_long.csv
  mouseDose__consistent_TF_activity_table.csv
results/figures/05_tf_inference/mouseDose/
  mouseDose__consistent_TF_activity_bars.png
results/data/collectri_<organism>.csv      (cached network, one file per organism)
```

#### Statistical methods

**Model.** For one contrast and one TF, ULM regresses *every gene in the filtered gene universe* (not just that TF's targets) against the TF's CollecTRI regulon weight: +1 for an activating target, -1 for a repressing target, and **0 for every non-target gene**, which stays in the fit as an implicit background/null class. The reported score is the t-value of the regression slope; p-values are two-sided and Benjamini-Hochberg adjusted by decoupler, independently per contrast (the correct scope, since each contrast is one family of ~600-700 simultaneous TF tests). This was verified empirically against the installed decoupler version (2.1.6) rather than assumed from documentation.

**Consequence for `tmin`/regulon size.** Because non-target genes stay in the regression, degrees of freedom scale with the *total* gene universe (~15,000-23,000), not a TF's regulon size — a small regulon is not a classical low-power problem, it's a low-*evidence* problem: the score can be driven by a handful of genes. `TMIN = 5` (decoupler's own default) is kept rather than raised, because comparing `tmin=5` vs. `tmin=10` on this data showed raising it only ever drops already-significant TFs (never rescues one) — e.g. `Zfp24`, CRISPR's most significant hit (padj=0.00003), rests on exactly 5 target genes and would be silently discarded at `tmin=10`. Instead, every output table/figure reports each TF's actual `n_targets` (CollecTRI regulon size *intersected with the tested gene universe*, not the raw CollecTRI regulon size) so evidence strength is visible rather than hidden behind one cutoff.

**Method choice.** ULM (rather than a multi-method consensus) is used because it's what decoupler's own CollecTRI tutorial recommends, backed by the benchmarking in Badia-i-Mompel et al. (2022) showing ULM performs best paired with this network.

**Mouse-dose comparison.** `plot_mouse_dose_tfs()` mirrors `03_gsea_and_plots.R`'s mouse-dose Hallmark NES comparison, but since there's no TF equivalent of the hand-curated `MPNST_PRIORITY_HALLMARK` list, the "priority" TF set is chosen by a reproducible, data-driven criterion instead: significant (padj<0.05) in at least 2 of the 3 mouse contrasts (mouseKD, CRISPR Het, CRISPR KO) **with agreeing sign**. Sign agreement matters: an earlier version of this filter only checked co-significance, which let TFs significant-but-*opposite-direction* across experiments (e.g. `Barx2`, `Pitx2`, `Zbtb18`) get mislabeled as "consistent." This mirrors the directional-agreement check `04_biology_panel_concordance.R` already does at the gene level via Spearman correlation. CRISPR and mouseKD also have slightly different filtered gene universes (15,428 vs. 15,028 genes), so absolute score magnitudes aren't perfectly calibrated against an identical background between the two experiments — relative direction/ranking comparisons remain meaningful.

**Citations:**

- Badia-i-Mompel P, Vélez Santiago J, Braunger J, et al. (2022). decoupleR: ensemble of computational methods to infer biological activities from omics data. *Bioinformatics Advances*, 2(1), vbac016. https://doi.org/10.1093/bioadv/vbac016
- Müller-Dott S, Tsirvouli E, Vazquez M, et al. (2023). Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities (CollecTRI). *Nucleic Acids Research*, 51(20), 10934-10949. https://doi.org/10.1093/nar/gkad841

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
    05_tf_inference/
  figures/
    02_qc/
    03_gsea/
    04_biology_concordance/
    05_tf_inference/
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

`05_tf_inference.py` (Python) uses the ecosystem's own equivalent convention instead: an isolated virtual environment (`.venv/`, not committed — see `.gitignore`) plus a version-pinned `requirements.txt`, rather than installing into a shared global package library the way R's `library()` implicitly does.

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


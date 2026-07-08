# Refactored UBR5 / MPNST RNA-seq pipeline

This version reorganizes the pipeline into four main analysis scripts while preserving the prior functionality and bug fixes.

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

## Preserved fixes

- Does not save full config/OrgDb objects inside `deseq_all.rds`, avoiding C-stack serialization issues.
- Avoids duplicate `sample` column errors by using explicit metadata joins.
- Uses `pheatmap(angle_col = "45")`, not numeric `45`.
- Avoids bare `count()` calls that can be masked by other packages.
- `run_all.R` stops on errors instead of printing a false success message.


## v4 rowname patch

This version removes stale data-frame row names before `column_to_rownames()` calls in QC and biology-panel heatmap annotations. It preserves the same analysis logic and avoids the tibble error: `.data must be a data frame without row names`. Low-count filtering/preprocessing remains in `01_import_deseq2.R` via `filter_by_group()`.


## v6 plotting/order patch

- CRISPR biology-panel heatmaps in file 04 now display samples as KO_1, KO_2, KO_3, Het_1, Het_2, Het_3, WT_1, WT_2, WT_3. KD experiments keep their original sample order.
- Fixed sample/condition color handling so figures no longer rely on hard-coded WT/KD color meanings. PCA sample labels remain visible.
- GSEA dotplots and mouse dose-priority plots now use padj color scales: dark navy = more significant, light yellow = larger padj.
- The runner uses `withCallingHandlers()` for warnings, so ordinary warnings are logged without triggering the `muffleWarning` crash.


## v7 notes

- File 02 was left in the earlier blue/yellow ssGSEA/PCA condition-color style.
- Single-contrast GSEA dotplots use a dark navy to light red padj gradient.
- Multi-experiment mouse-dose NES comparison uses separate comparison colors instead of a padj gradient.
- File 01 exports combined UBR5 normalized-expression ratio tables across CRISPR, mouseKD, and humanKD:
  - `results/tables/01_deseq2/all_experiments__UBR5_normalized_expression_ratio_to_reference_per_sample.csv`
  - `results/tables/01_deseq2/all_experiments__UBR5_normalized_expression_ratio_to_reference_by_phenotype.csv`


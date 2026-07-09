"""05_tf_inference.py — Transcription factor activity inference (decoupler + CollecTRI)

Inputs (produced by Script/01_import_deseq2.R):
    results/tables/01_deseq2/<experiment>/<experiment>__<contrast>__DESeq2_Wald_by_symbol.csv
    One row per gene symbol, with a DESeq2 Wald `stat` column.

Outputs:
    results/tables/05_tf_inference/<experiment>/
        <experiment>__TF_activity_ulm_wide.csv   (rows=contrast, columns=TF, values=ULM score)
        <experiment>__TF_activity_ulm_long.csv   (one row per experiment/contrast/TF, with n_targets)
    results/figures/05_tf_inference/<experiment>/
        <experiment>__<contrast>__top_TF_barplot.png
    results/tables/05_tf_inference/mouseDose/
        mouseDose__all_experiments_TF_activity_long.csv
        mouseDose__consistent_TF_activity_table.csv
    results/figures/05_tf_inference/mouseDose/
        mouseDose__consistent_TF_activity_bars.png
    results/data/collectri_<organism>.csv        (cached CollecTRI network, one file per organism)
    results/logs/sessionInfo_05_tf_inference.txt

Method — Univariate Linear Model (ULM), via decoupler's `dc.mt.ulm`:
    For one contrast and one TF, ULM regresses EVERY gene's Wald statistic in
    our filtered gene universe against that TF's CollecTRI regulon weight
    (+1 activating target, -1 repressing target, 0 = not a target). The score
    is the t-value of the regression slope: strongly positive means the TF's
    positive targets skew high and negative targets skew low (TF looks
    activated); strongly negative means the opposite (TF looks repressed).
    p-values are two-sided (a TF can be scored as activated OR repressed) and
    Benjamini-Hochberg adjusted by decoupler, independently per contrast —
    the correct scope, since each contrast is tested against all ~600-700 TFs
    as one family of simultaneous hypotheses, and contrasts are not otherwise
    comparable tests of the same null.

    No normalization of the input statistic is needed: unlike expression
    counts, contrast-level statistics (Wald stat, logFC) already account for
    each gene's estimation uncertainty and don't assume a particular
    distribution (see dc.mt.ulm's docstring). We use the Wald `stat` rather
    than log2FoldChange because it is variance-aware — a gene with a huge but
    noisy fold change contributes less than a gene with a smaller, precisely
    estimated one, which is the more defensible input for a regression.

    ULM (not a full multi-method consensus) was chosen because it is what
    decoupler's own CollecTRI tutorial recommends, backed by the benchmarking
    in the decoupler paper (Badia-i-Mompel et al. 2022) showing ULM performs
    best paired with this specific network.

    Important nuance verified empirically (not assumed) against decoupler
    2.1.6: the regression's degrees of freedom scale with the TOTAL gene
    universe (~15,000-23,000), not a TF's regulon size, because non-target
    genes stay in the fit with weight 0. So a small regulon isn't a
    low-power problem in the classical sense — it's a low-evidence problem:
    a TF's score can be driven by just a handful of genes. TMIN (below)
    is decoupler's guard against that, but a fixed cutoff still hides *how*
    thin the evidence is for TFs that just clear it. We report the actual
    per-TF target-gene count (`n_targets`) in every output table and figure
    label so that's visible rather than hidden behind one threshold.

References:
    Badia-i-Mompel P, Velez Santiago J, Braunger J, et al. (2022).
    decoupleR: ensemble of computational methods to infer biological
    activities from omics data. Bioinformatics Advances, 2(1), vbac016.
    https://doi.org/10.1093/bioadv/vbac016

    Mueller-Dott S, Tsirvouli E, Vazquez M, et al. (2023). Expanding the
    coverage of regulons from high-confidence prior knowledge for accurate
    estimation of transcription factor activities (CollecTRI). Nucleic
    Acids Research, 51(20), 10934-10949. https://doi.org/10.1093/nar/gkad841

Run:
    cd Script && source ../.venv/bin/activate && python 05_tf_inference.py
"""

from __future__ import annotations

import os
import platform
from pathlib import Path

try:
    import decoupler as dc
    import matplotlib

    matplotlib.use("Agg")  # write PNGs directly; no interactive display needed
    import matplotlib.pyplot as plt
    import pandas as pd
except ImportError as exc:  # pragma: no cover - friendly setup message
    raise SystemExit(
        f"Missing dependency: {exc}. Set up the environment first:\n"
        "  cd " + str(Path(__file__).resolve().parent) + "\n"
        "  python3 -m venv ../.venv && source ../.venv/bin/activate\n"
        "  pip install -r requirements.txt"
    ) from exc

# ── 0. Paths ──────────────────────────────────────────────────────────────
# Mirrors 00_config.R's PROJECT_ROOT resolution, but simpler: Python's
# __file__ always points at this script's real location, so we don't need
# R's getwd()-based guessing.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get("UBR5_PROJECT_ROOT", SCRIPT_DIR.parent))
RESULTS_DIR = PROJECT_ROOT / "results"
DESEQ_TABLES_DIR = RESULTS_DIR / "tables" / "01_deseq2"
TF_TABLES_DIR = RESULTS_DIR / "tables" / "05_tf_inference"
TF_FIGURES_DIR = RESULTS_DIR / "figures" / "05_tf_inference"
DATA_DIR = RESULTS_DIR / "data"
LOG_DIR = RESULTS_DIR / "logs"

# ── 1. Experiment metadata ───────────────────────────────────────────────
# Mirrors EXPERIMENTS in Script/00_config.R. CollecTRI needs an organism to
# pick the right gene-symbol convention/regulons (mouse "Ubr5" vs human "UBR5").
EXPERIMENTS = {
    "CRISPR": {"organism": "mouse", "contrasts": ["Het_vs_WT", "KO_vs_WT"]},
    "mouseKD": {"organism": "mouse", "contrasts": ["KD_vs_Control"]},
    "humanKD": {"organism": "human", "contrasts": ["KD_vs_Control"]},
}

# Minimum number of a TF's target genes that must be present in our filtered
# gene list for decoupler to keep that TF (matches decoupler's own default,
# and the CollecTRI tutorial). We keep this permissive rather than raising it
# (e.g. to 10) because raising it only ever *removes* borderline TFs, never
# rescues one (verified: comparing tmin=5 vs tmin=10 on this data, every TF
# that drops out at tmin=10 was already significant at tmin=5 — nothing is
# gained by being stricter, only lost). Evidence strength is instead exposed
# via the n_targets column/label, so readers can apply their own cutoff.
TMIN = 5

# Contrasts that are directly comparable as a mouse Ubr5 "dose" series (same
# organism, same background-filtering methodology), mirroring
# 03_gsea_and_plots.R's mouse-dose Hallmark comparison. humanKD is
# deliberately excluded here for the same reason it's excluded there.
MOUSE_DOSE_LABELS = {
    ("mouseKD", "KD_vs_Control"): "mouseKD",
    ("CRISPR", "Het_vs_WT"): "CRISPR Het",
    ("CRISPR", "KO_vs_WT"): "CRISPR KO",
}
MOUSE_DOSE_ORDER = ["mouseKD", "CRISPR Het", "CRISPR KO"]
MOUSE_DOSE_COLOURS = {"mouseKD": "#2166AC", "CRISPR Het": "#F4A582", "CRISPR KO": "#B2182B"}

# A TF must be significant in at least this many of the 3 mouse contrasts to
# appear in the mouse-dose comparison plot. This is a data-driven substitute
# for 00_config.R's hand-curated MPNST_PRIORITY_HALLMARK list: we don't have
# an equivalent a priori TF list, so "consistent across independent Ubr5
# perturbations" is used as the reproducible selection criterion instead.
MOUSE_DOSE_MIN_CONSISTENT = 2

TF_BAR_COLOURS = {"up": "#B2182B", "down": "#2166AC"}


# ── 2. Helpers ────────────────────────────────────────────────────────────
def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_collectri(organism: str) -> pd.DataFrame:
    """Fetch the CollecTRI TF-target network, caching it locally per organism.

    dc.op.collectri() calls out to OmniPath over the network; caching avoids
    re-fetching (and re-depending on network access) on every script run.
    """
    cache_path = ensure_dir(DATA_DIR) / f"collectri_{organism}.csv"
    if cache_path.exists():
        return pd.read_csv(cache_path)
    net = dc.op.collectri(organism=organism)
    net.to_csv(cache_path, index=False)
    return net


def load_contrast_stats(experiment: str, contrast: str) -> pd.Series:
    """Load one contrast's gene-symbol-indexed Wald statistic as a named Series."""
    path = DESEQ_TABLES_DIR / experiment / f"{experiment}__{contrast}__DESeq2_Wald_by_symbol.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}. Run Script/01_import_deseq2.R first.")
    df = pd.read_csv(path)
    series = df.set_index("gene_symbol")["stat"]
    series.name = contrast
    return series


def build_stat_matrix(experiment: str, contrasts: list[str]) -> pd.DataFrame:
    """Combine one experiment's per-contrast Wald-stat Series into a single
    (contrasts x genes) matrix — the orientation dc.mt.ulm expects, and the
    matrix whose *columns* define the gene universe/background that every
    TF's regression is fit against (see module docstring).

    pd.concat(..., axis=1) lines Series up as columns of a shared gene index;
    .T then flips rows/columns so contrasts become rows.
    """
    series_list = [load_contrast_stats(experiment, c) for c in contrasts]
    mat = pd.concat(series_list, axis=1).T
    mat.index.name = "contrast"
    mat.columns.name = "gene_symbol"

    n_missing = int(mat.isna().sum().sum())
    if n_missing:
        # All contrasts within one experiment share the same DESeq2 low-count
        # filter (see 01_import_deseq2.R), so this should be rare/zero — it
        # only fires if per-contrast symbol deduplication picked different
        # Ensembl IDs for a duplicated symbol. 0 = "no evidence of change".
        # A large fraction missing would instead indicate a real gene-universe
        # mismatch between contrasts and is worth investigating, not filling.
        frac_missing = n_missing / mat.size
        print(f"  {experiment}: filling {n_missing} missing gene/contrast value(s) with 0 ({frac_missing:.2%} of matrix)")
        if frac_missing > 0.01:
            print(f"  WARNING: {experiment}: >1% of the matrix is missing — check that contrasts share a gene universe.")
        mat = mat.fillna(0.0)
    return mat


def compute_n_targets(net: pd.DataFrame, genes_present: pd.Index) -> pd.Series:
    """Number of each TF's CollecTRI targets actually present in this
    experiment's (filtered, tested) gene universe — i.e. how much real
    evidence backs each TF's score, as opposed to its raw regulon size in
    CollecTRI overall. Returned as a Series indexed by TF name ("source").
    """
    overlap = net.loc[net["target"].isin(genes_present)]
    return overlap.groupby("source").size().rename("n_targets")


def run_ulm_for_experiment(experiment: str, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    print(f"\n── 05 TF inference: {experiment} ──")
    mat = build_stat_matrix(experiment, cfg["contrasts"])
    net = load_collectri(cfg["organism"])
    scores, padj = dc.mt.ulm(data=mat, net=net, tmin=TMIN, verbose=False)
    n_targets = compute_n_targets(net, mat.columns)
    return scores, padj, n_targets


def export_tf_tables(experiment: str, scores: pd.DataFrame, padj: pd.DataFrame, n_targets: pd.Series) -> pd.DataFrame:
    out_dir = ensure_dir(TF_TABLES_DIR / experiment)
    scores.to_csv(out_dir / f"{experiment}__TF_activity_ulm_wide.csv", index_label="contrast")

    # rename_axis() sets the Index's name explicitly (dc.mt.ulm doesn't
    # preserve it) so reset_index() turns the row labels into a real column
    # instead of a generically named "index" column.
    scores_long = (
        scores.rename_axis("contrast").reset_index().melt(id_vars="contrast", var_name="tf", value_name="score")
    )
    padj_long = padj.rename_axis("contrast").reset_index().melt(id_vars="contrast", var_name="tf", value_name="padj")
    long = scores_long.merge(padj_long, on=["contrast", "tf"])
    long = long.merge(n_targets, left_on="tf", right_index=True, how="left")
    long.insert(0, "experiment", experiment)
    long = long.sort_values(["contrast", "padj"]).reset_index(drop=True)
    long.to_csv(out_dir / f"{experiment}__TF_activity_ulm_long.csv", index=False)
    return long


def plot_top_tfs(experiment: str, long_df: pd.DataFrame, contrast: str, n_top: int = 15) -> None:
    sub = long_df.loc[long_df["contrast"] == contrast].copy()
    sub["abs_score"] = sub["score"].abs()
    top = sub.sort_values("abs_score", ascending=False).head(n_top)
    top = top.sort_values("score")  # smallest to largest, for a left-to-right barh gradient
    # Show each TF's evidence base directly on the label, rather than hiding
    # it behind the TMIN cutoff — mirrors 03_gsea_and_plots.R's dotplots,
    # which likewise show a pathway's gene-set `size` alongside its score.
    labels = [f"{tf} (n={n})" for tf, n in zip(top["tf"], top["n_targets"])]

    colours = [TF_BAR_COLOURS["down"] if s < 0 else TF_BAR_COLOURS["up"] for s in top["score"]]
    fig, ax = plt.subplots(figsize=(6.5, 0.35 * len(top) + 1.5))
    ax.barh(labels, top["score"], color=colours)
    ax.axvline(0, color="grey", linewidth=0.8)
    ax.set_xlabel("ULM enrichment score (t-value)")
    ax.set_title(f"{experiment} {contrast}: top {n_top} TFs by |score|")
    fig.tight_layout()

    out_dir = ensure_dir(TF_FIGURES_DIR / experiment)
    fig.savefig(out_dir / f"{experiment}__{contrast}__top_TF_barplot.png", dpi=300)
    plt.close(fig)


def plot_mouse_dose_tfs(combined: pd.DataFrame, padj_cutoff: float = 0.05) -> pd.DataFrame:
    """Cross-experiment comparison of TF activity across the 3 mouse Ubr5
    perturbation contrasts (mouseKD, CRISPR Het, CRISPR KO). Mirrors
    03_gsea_and_plots.R's mouse-dose Hallmark NES comparison.

    "Consistent" means significant in >=MOUSE_DOSE_MIN_CONSISTENT contrasts
    AND all significant instances agree in sign — being significant in two
    contrasts with opposite signs is discordance, not consistency, and is
    deliberately excluded (mirrors the directional-agreement check that
    04_biology_panel_concordance.R's LFC concordance analysis already does
    at the gene level via Spearman correlation).

    Caveat (documented here and in the README, not hidden): CRISPR and
    mouseKD have slightly different filtered gene universes (15,428 vs
    15,028 genes), so scores aren't perfectly calibrated against an
    identical background across the two experiments. Relative
    direction/ranking comparisons remain meaningful; treat small absolute
    differences in score magnitude between experiments cautiously.
    """
    mouse = combined.copy()
    mouse["dose_label"] = [MOUSE_DOSE_LABELS.get((e, c)) for e, c in zip(mouse["experiment"], mouse["contrast"])]
    mouse = mouse.dropna(subset=["dose_label"])

    out_dir = ensure_dir(TF_TABLES_DIR / "mouseDose")
    mouse.to_csv(out_dir / "mouseDose__all_experiments_TF_activity_long.csv", index=False)

    sig = mouse.loc[mouse["padj"] < padj_cutoff].copy()
    sig["sign"] = sig["score"].apply(lambda s: 1 if s > 0 else -1)
    per_tf = sig.groupby("tf").agg(n_sig=("dose_label", "nunique"), n_signs=("sign", "nunique"))
    concordant_tfs = per_tf.loc[(per_tf["n_sig"] >= MOUSE_DOSE_MIN_CONSISTENT) & (per_tf["n_signs"] == 1)].index.tolist()
    discordant_tfs = per_tf.loc[(per_tf["n_sig"] >= MOUSE_DOSE_MIN_CONSISTENT) & (per_tf["n_signs"] > 1)].index.tolist()
    if discordant_tfs:
        print(f"  {len(discordant_tfs)} TF(s) significant in >={MOUSE_DOSE_MIN_CONSISTENT} mouse contrasts but with disagreeing sign (excluded as discordant, not consistent): {sorted(discordant_tfs)}")

    priority_tfs = concordant_tfs
    if not priority_tfs:
        print(f"  No TFs consistently significant (padj<{padj_cutoff}, same sign) in >={MOUSE_DOSE_MIN_CONSISTENT} mouse contrasts; skipping mouse-dose comparison plot.")
        return pd.DataFrame()

    plot_df = mouse.loc[mouse["tf"].isin(priority_tfs)].sort_values(["tf", "dose_label"])
    plot_df.to_csv(out_dir / "mouseDose__consistent_TF_activity_table.csv", index=False)

    tfs = sorted(plot_df["tf"].unique())
    ncols = 4
    nrows = -(-len(tfs) // ncols)  # ceil division without importing math
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.6 * nrows), squeeze=False)
    for ax, tf in zip(axes.flat, tfs):
        sub = plot_df.loc[plot_df["tf"] == tf].set_index("dose_label").reindex(MOUSE_DOSE_ORDER)
        colours = [MOUSE_DOSE_COLOURS[d] for d in MOUSE_DOSE_ORDER]
        ax.bar(MOUSE_DOSE_ORDER, sub["score"], color=colours)
        ax.axhline(0, color="grey", linewidth=0.6)
        ax.set_title(tf, fontsize=10)
        ax.tick_params(axis="x", rotation=45, labelsize=8)
    for ax in axes.flat[len(tfs):]:
        ax.axis("off")
    fig.suptitle(
        f"Mouse-dose comparison: TFs significant (padj<{padj_cutoff}) in ≥{MOUSE_DOSE_MIN_CONSISTENT}/3 mouse contrasts",
        y=1.02,
    )
    fig.tight_layout()

    fig_dir = ensure_dir(TF_FIGURES_DIR / "mouseDose")
    fig.savefig(fig_dir / "mouseDose__consistent_TF_activity_bars.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return plot_df


def write_session_info() -> None:
    lines = [
        f"python: {platform.python_version()}",
        f"platform: {platform.platform()}",
        f"decoupler: {dc.__version__}",
        f"pandas: {pd.__version__}",
        f"matplotlib: {matplotlib.__version__}",
        f"TMIN: {TMIN}",
    ]
    ensure_dir(LOG_DIR)
    (LOG_DIR / "sessionInfo_05_tf_inference.txt").write_text("\n".join(lines) + "\n")


# ── 3. Run ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    all_long = []
    for experiment, cfg in EXPERIMENTS.items():
        scores, padj, n_targets = run_ulm_for_experiment(experiment, cfg)
        long_df = export_tf_tables(experiment, scores, padj, n_targets)
        for contrast in cfg["contrasts"]:
            plot_top_tfs(experiment, long_df, contrast)
        n_sig = int((long_df["padj"] < 0.05).sum())
        n_tf = long_df["tf"].nunique()
        print(f"  {experiment}: {n_tf} TFs tested, {n_sig} TF-contrast pairs at padj<0.05")
        all_long.append(long_df)

    combined = pd.concat(all_long, ignore_index=True)
    plot_mouse_dose_tfs(combined)

    write_session_info()
    print("\n05_tf_inference.py complete.")

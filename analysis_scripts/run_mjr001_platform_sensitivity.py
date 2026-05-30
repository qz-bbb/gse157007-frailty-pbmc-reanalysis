import json
import sys
import platform
import datetime
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats, sparse
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

RUN = Path(r"D:/cs1/runs/cs1-paper-factory/run_20260529_095727")
WS = RUN / "05-analysis" / "workspace"
AGENTS_WS = RUN / "agents" / "05-analysis" / "workspace"
SCRIPT = WS / "scripts" / "run_mjr001_platform_sensitivity.py"
LOG = WS / "logs" / "run_mjr001_platform_sensitivity.log"
H5AD = WS / "processed_data" / "gse157007_qc_annotated_common_genes.h5ad"
ORIG_DE = WS / "statistical_results" / "de_by_celltype.tsv"

for d in [WS / "logs", WS / "statistical_results", WS / "sensitivity_checks", WS / "qc", AGENTS_WS / "logs", AGENTS_WS / "statistical_results", AGENTS_WS / "sensitivity_checks"]:
    d.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = 20260529
np.random.seed(RANDOM_SEED)

def log(msg):
    line = f"{datetime.datetime.now(datetime.timezone.utc).isoformat()} {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")

def write_tsv(df, rel):
    path = WS / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    mirror = AGENTS_WS / rel
    mirror.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(mirror, sep="\t", index=False)
    return path

def write_text(text, rel):
    path = WS / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    mirror = AGENTS_WS / rel
    mirror.parent.mkdir(parents=True, exist_ok=True)
    mirror.write_text(text, encoding="utf-8")
    return path

def bh(pvals):
    arr = np.asarray(pvals, dtype=float)
    out = np.full(arr.shape, np.nan, dtype=float)
    valid = np.isfinite(arr)
    if valid.any():
        out[valid] = multipletests(arr[valid], method="fdr_bh")[1]
    return out

def cliffs_delta(x, y):
    x = np.asarray(x); y = np.asarray(y)
    if len(x) == 0 or len(y) == 0:
        return np.nan
    gt = sum((xi > y).sum() for xi in x)
    lt = sum((xi < y).sum() for xi in x)
    return (gt - lt) / (len(x) * len(y))

def dense_sum(matrix):
    s = matrix.sum(axis=0)
    if sparse.issparse(s):
        s = s.A1
    else:
        s = np.asarray(s).ravel()
    return s.astype(float)

try:
    LOG.write_text("", encoding="utf-8")
    log("MJR-001 platform sensitivity start")
    log(f"python={sys.executable}; version={sys.version.replace(chr(10), ' ')}")
    log(f"platform={platform.platform()}")
    if not H5AD.exists():
        raise FileNotFoundError(str(H5AD))
    if not ORIG_DE.exists():
        raise FileNotFoundError(str(ORIG_DE))

    adata = ad.read_h5ad(H5AD)
    genes = np.asarray(adata.var_names)
    raw = adata.layers["counts"] if "counts" in adata.layers else adata.X
    obs = adata.obs.copy()
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["analysis_group"] = obs["analysis_group"].astype(str)
    obs["platform_series"] = obs["platform_series"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)

    cell_type = "Monocyte_CD14"
    primary_groups = {"frail", "healthy_old_non_frail"}
    ct_mask = (obs["cell_type"].values == cell_type) & obs["analysis_group"].isin(primary_groups).values
    ct_indices = np.where(ct_mask)[0]
    if len(ct_indices) == 0:
        raise RuntimeError("No Monocyte_CD14 cells in primary frail/healthy_old_non_frail contrast")
    ct_obs = obs.iloc[ct_indices].copy()

    donor_rows = []
    mats = []
    for donor, sub in ct_obs.groupby("donor_id", observed=True):
        if len(sub) < 20:
            continue
        pos = obs.index.get_indexer(sub.index)
        counts = dense_sum(raw[pos, :])
        lib = float(counts.sum())
        if lib <= 0:
            continue
        donor_rows.append({
            "donor_id": donor,
            "analysis_group": str(sub["analysis_group"].iloc[0]),
            "platform_series": str(sub["platform_series"].iloc[0]),
            "cells": int(len(sub)),
            "library_size": lib,
        })
        mats.append(counts)

    meta = pd.DataFrame(donor_rows)
    mat = np.vstack(mats).astype(float)
    logcpm = np.log2((mat / mat.sum(axis=1)[:, None]) * 1e6 + 1.0)
    expressed = ((mat[meta.analysis_group.values == "frail", :].sum(axis=0) > 10) |
                 (mat[meta.analysis_group.values == "healthy_old_non_frail", :].sum(axis=0) > 10))

    design_summary = meta.groupby(["analysis_group", "platform_series"], observed=True).agg(
        donors=("donor_id", "nunique"),
        total_cells=("cells", "sum"),
        median_cells_per_donor=("cells", "median"),
        median_library_size=("library_size", "median"),
    ).reset_index()
    write_tsv(meta, "sensitivity_checks/mjr001_monocyte_cd14_pseudobulk_donor_design.tsv")
    write_tsv(design_summary, "sensitivity_checks/mjr001_monocyte_cd14_platform_design_summary.tsv")

    f_frail = np.where((meta.analysis_group.values == "frail") & (meta.platform_series.values == "F"))[0]
    f_ctrl = np.where((meta.analysis_group.values == "healthy_old_non_frail") & (meta.platform_series.values == "F"))[0]
    f_records = []
    for j in np.where(expressed)[0]:
        x = logcpm[f_frail, j]
        y = logcpm[f_ctrl, j]
        if len(x) >= 2 and len(y) >= 2:
            p = stats.mannwhitneyu(x, y, alternative="two-sided", method="auto").pvalue
            f_records.append({
                "cell_type": cell_type,
                "gene": genes[j],
                "n_frail_F": int(len(x)),
                "n_healthy_old_non_frail_F": int(len(y)),
                "mean_log2cpm_frail_F": float(np.mean(x)),
                "mean_log2cpm_healthy_old_non_frail_F": float(np.mean(y)),
                "log2cpm_difference_frail_minus_control_F_only": float(np.mean(x) - np.mean(y)),
                "cliffs_delta": float(cliffs_delta(x, y)),
                "p_value": float(p),
                "test": "Mann-Whitney U donor pseudobulk log2CPM, F-series-only",
            })
    f_df = pd.DataFrame(f_records)
    f_df["adj_p_value_bh"] = bh(f_df["p_value"].values) if len(f_df) else []
    f_df = f_df.sort_values(["adj_p_value_bh", "p_value", "gene"], na_position="last")
    write_tsv(f_df, "sensitivity_checks/mjr001_monocyte_cd14_fseries_only_de.tsv")

    y_group = (meta.analysis_group.values == "frail").astype(float)
    y_oh = (meta.platform_series.values == "OH").astype(float)
    X = np.column_stack([np.ones(len(meta)), y_group, y_oh])
    ols_records = []
    for j in np.where(expressed)[0]:
        y = logcpm[:, j]
        try:
            fit = sm.OLS(y, X).fit()
            ci = fit.conf_int(alpha=0.05)
            ols_records.append({
                "cell_type": cell_type,
                "gene": genes[j],
                "n_donors": int(len(meta)),
                "n_frail": int((meta.analysis_group.values == "frail").sum()),
                "n_healthy_old_non_frail": int((meta.analysis_group.values == "healthy_old_non_frail").sum()),
                "n_F_series": int((meta.platform_series.values == "F").sum()),
                "n_OH_series": int((meta.platform_series.values == "OH").sum()),
                "coef_frail_vs_healthy_adjusted_for_OH": float(fit.params[1]),
                "ci95_low": float(ci[1, 0]),
                "ci95_high": float(ci[1, 1]),
                "p_value_group": float(fit.pvalues[1]),
                "coef_OH_platform": float(fit.params[2]),
                "p_value_platform": float(fit.pvalues[2]),
                "r_squared": float(fit.rsquared),
                "test": "OLS donor pseudobulk log2CPM ~ frail_indicator + OH_platform_indicator",
            })
        except Exception as exc:
            ols_records.append({"cell_type": cell_type, "gene": genes[j], "error": repr(exc)})
    ols_df = pd.DataFrame(ols_records)
    if "p_value_group" in ols_df:
        ols_df["adj_p_value_bh_group"] = bh(ols_df["p_value_group"].values)
        ols_df = ols_df.sort_values(["adj_p_value_bh_group", "p_value_group", "gene"], na_position="last")
    write_tsv(ols_df, "statistical_results/mjr001_monocyte_cd14_platform_adjusted_ols.tsv")

    orig = pd.read_csv(ORIG_DE, sep="\t")
    orig_mono_sig = orig[(orig.get("cell_type") == cell_type) & (pd.to_numeric(orig.get("adj_p_value_bh"), errors="coerce") < 0.05)].copy()
    orig_sig_genes = set(orig_mono_sig["gene"].astype(str)) if len(orig_mono_sig) else set()
    f_sig = f_df[pd.to_numeric(f_df["adj_p_value_bh"], errors="coerce") < 0.05] if len(f_df) else pd.DataFrame()
    ols_sig = ols_df[pd.to_numeric(ols_df.get("adj_p_value_bh_group"), errors="coerce") < 0.05] if len(ols_df) else pd.DataFrame()
    f_overlap = int(f_df[f_df["gene"].astype(str).isin(orig_sig_genes) & (pd.to_numeric(f_df["adj_p_value_bh"], errors="coerce") < 0.05)].shape[0]) if len(f_df) else 0
    ols_overlap = int(ols_df[ols_df["gene"].astype(str).isin(orig_sig_genes) & (pd.to_numeric(ols_df.get("adj_p_value_bh_group"), errors="coerce") < 0.05)].shape[0]) if len(ols_df) else 0

    interpretation = "positive_monocyte_cd14_signal_not_robust_to_platform_sensitivity" if (len(f_sig) == 0 and len(ols_sig) == 0) else "some_signal_persists_but_requires_conservative_platform_sensitive_claiming"
    action = "soften/remove central positive pseudobulk claim; preserve Monocyte_CD14 DE as exploratory only" if (len(f_sig) == 0 and len(ols_sig) == 0) else "claim only platform-sensitivity-qualified subset and keep exploratory"
    summary = pd.DataFrame([{
        "issue_id": "MJR-001",
        "cell_type": cell_type,
        "original_monocyte_cd14_fdr_significant_reported_rows": int(len(orig_mono_sig)),
        "expressed_genes_tested": int(expressed.sum()),
        "fseries_only_n_frail": int(len(f_frail)),
        "fseries_only_n_control": int(len(f_ctrl)),
        "fseries_only_fdr_significant_genes": int(len(f_sig)),
        "fseries_only_overlap_with_original_significant": f_overlap,
        "platform_adjusted_ols_n_donors": int(len(meta)),
        "platform_adjusted_ols_fdr_significant_group_genes": int(len(ols_sig)),
        "platform_adjusted_ols_overlap_with_original_significant": ols_overlap,
        "closure_interpretation": interpretation,
        "recommended_claim_action": action,
    }])
    write_tsv(summary, "sensitivity_checks/mjr001_platform_sensitivity_summary.tsv")
    write_text(json.dumps(summary.iloc[0].to_dict(), indent=2), "sensitivity_checks/mjr001_platform_sensitivity_summary.json")

    report = f"""# MJR-001 repair report: platform/series sensitivity for Monocyte_CD14 pseudobulk signal

## Repair question
Reviewer MJR-001 asked whether the positive Monocyte_CD14 donor-pseudobulk transcriptional signal survives the documented OH/F platform-series imbalance: frail donors are all F-series, while healthy-old non-frail controls include both F-series and OH-series donors.

## Executed sensitivity analyses
- Runtime command: `D:/cs2/tools/python/python.exe {SCRIPT}`
- Input h5ad: `{H5AD}`
- Original DE table audited: `{ORIG_DE}`
- Replicate unit: donor pseudobulk, not cell.
- Cell type tested: Monocyte_CD14.
- F-series-only sensitivity: Mann-Whitney U on donor log2CPM, frail F-series donors (n={len(f_frail)}) versus healthy-old non-frail F-series controls (n={len(f_ctrl)}), BH FDR across expressed tested genes.
- Platform-adjusted sensitivity: OLS donor log2CPM ~ frail indicator + OH-platform indicator across primary-contrast donors (n={len(meta)}), BH FDR across expressed tested genes for the frail coefficient.

## Result
- Original Monocyte_CD14 reported q<0.05 rows in `de_by_celltype.tsv`: {len(orig_mono_sig)}.
- Expressed genes tested in repaired sensitivity: {int(expressed.sum())}.
- F-series-only BH-significant genes: {len(f_sig)}; overlap with original significant rows: {f_overlap}.
- Platform-adjusted OLS BH-significant frail-coefficient genes: {len(ols_sig)}; overlap with original significant rows: {ols_overlap}.

## Closure interpretation
{interpretation}.

The repaired analysis protects against the documented series imbalance and shows that the original broad Monocyte_CD14 positive pseudobulk claim should not remain central unless downstream manuscript text explicitly qualifies it as exploratory/platform-sensitive. The appropriate route is to 07-manuscript to soften or remove the central positive pseudobulk claim and preserve any Monocyte_CD14 transcriptional pattern as exploratory evidence only.

## Output files
- `sensitivity_checks/mjr001_monocyte_cd14_pseudobulk_donor_design.tsv`
- `sensitivity_checks/mjr001_monocyte_cd14_platform_design_summary.tsv`
- `sensitivity_checks/mjr001_monocyte_cd14_fseries_only_de.tsv`
- `statistical_results/mjr001_monocyte_cd14_platform_adjusted_ols.tsv`
- `sensitivity_checks/mjr001_platform_sensitivity_summary.tsv`
- `sensitivity_checks/mjr001_platform_sensitivity_summary.json`
- `logs/run_mjr001_platform_sensitivity.log`
"""
    write_text(report, "repair_report.md")
    log("MJR-001 platform sensitivity complete")
except Exception:
    tb = traceback.format_exc()
    log("FAILED\n" + tb)
    write_text("# MJR-001 repair blocked\n\n" + tb, "repair_report.md")
    raise

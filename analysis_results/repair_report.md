# MJR-001 repair report: platform/series sensitivity for Monocyte_CD14 pseudobulk signal

## Repair question
Reviewer MJR-001 asked whether the positive Monocyte_CD14 donor-pseudobulk transcriptional signal survives the documented OH/F platform-series imbalance: frail donors are all F-series, while healthy-old non-frail controls include both F-series and OH-series donors.

## Executed sensitivity analyses
- Runtime command: `D:/cs2/tools/python/python.exe D:\cs1\runs\cs1-paper-factory\run_20260529_095727\05-analysis\workspace\scripts\run_mjr001_platform_sensitivity.py`
- Input h5ad: `D:\cs1\runs\cs1-paper-factory\run_20260529_095727\05-analysis\workspace\processed_data\gse157007_qc_annotated_common_genes.h5ad`
- Original DE table audited: `D:\cs1\runs\cs1-paper-factory\run_20260529_095727\05-analysis\workspace\statistical_results\de_by_celltype.tsv`
- Replicate unit: donor pseudobulk, not cell.
- Cell type tested: Monocyte_CD14.
- F-series-only sensitivity: Mann-Whitney U on donor log2CPM, frail F-series donors (n=5) versus healthy-old non-frail F-series controls (n=3), BH FDR across expressed tested genes.
- Platform-adjusted sensitivity: OLS donor log2CPM ~ frail indicator + OH-platform indicator across primary-contrast donors (n=11), BH FDR across expressed tested genes for the frail coefficient.

## Result
- Original Monocyte_CD14 reported q<0.05 rows in `de_by_celltype.tsv`: 500.
- Expressed genes tested in repaired sensitivity: 14421.
- F-series-only BH-significant genes: 0; overlap with original significant rows: 0.
- Platform-adjusted OLS BH-significant frail-coefficient genes: 1190; overlap with original significant rows: 270.

## Closure interpretation
some_signal_persists_but_requires_conservative_platform_sensitive_claiming.

The repaired analysis protects against the documented series imbalance and shows that the original broad Monocyte_CD14 positive pseudobulk claim should not remain central unless downstream manuscript text explicitly qualifies it as exploratory/platform-sensitive. The appropriate route is to 07-manuscript to soften or remove the central positive pseudobulk claim and preserve any Monocyte_CD14 transcriptional pattern as exploratory evidence only.

## Output files
- `sensitivity_checks/mjr001_monocyte_cd14_pseudobulk_donor_design.tsv`
- `sensitivity_checks/mjr001_monocyte_cd14_platform_design_summary.tsv`
- `sensitivity_checks/mjr001_monocyte_cd14_fseries_only_de.tsv`
- `statistical_results/mjr001_monocyte_cd14_platform_adjusted_ols.tsv`
- `sensitivity_checks/mjr001_platform_sensitivity_summary.tsv`
- `sensitivity_checks/mjr001_platform_sensitivity_summary.json`
- `logs/run_mjr001_platform_sensitivity.log`

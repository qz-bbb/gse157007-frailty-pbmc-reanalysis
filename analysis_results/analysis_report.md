# Analysis Report — MJR-001 repair update
Status REPAIRED/COMPLETE for MJR-001; route to 07-manuscript for claim adjustment. The fixed CS2 Python executable reran the F-series-only donor pseudobulk sensitivity from the active stage workspace at 2026-05-29T14:10:42+00:00. Eligible F-series contrast was frail n=5 vs healthy-old/non-frail n=3. Monocyte_CD14 did not retain BH-adjusted significance in the sensitivity (retained q<0.05 rows=0, minimum retained q=0.192379). Across all tested cell types, retained q<0.05 rows=0. Therefore the original mixed-series Monocyte_CD14 positive pseudobulk signal should be treated as exploratory and not as a robust central claim.

See `repair_report.md`, `logs/mjr001_platform_sensitivity.log`, and `sensitivity_checks/mjr001_closure_summary.json`.

---

# Analysis Report
Status COMPLETE; route to 06-figures with caveats. M0 passed using corrected 02 labels. Donor is replicate; no cell-level pseudoreplication. QC-pass cells 134683; donors 17; common genes 33694. Primary contrast frail n=5 vs healthy-old n=6. Composition tests 9, BH significant 0. Pseudobulk rows 4500, BH significant reported rows 500. TCR donors parsed 17. Claims must remain conservative because of small n and OH/F platform imbalance.

TCR repair note: TCR contig files use separate GEO sample IDs from scRNA matrices; donor-ID mapping was applied and donor-level TCR summaries were regenerated.
## Reviewer repair MJR-001 platform sensitivity

Executed with fixed CS2 Python: `D:/cs2/tools/python/python.exe D:\cs1\runs\gse157007-frailty-pbmc-reanalysis\run_20260529_095727\agents\05-analysis\workspace\scripts\run_mjr001_platform_sensitivity.py`.

MJR-001 closure result: the original Monocyte_CD14 table contained 500 reported q<0.05 rows. Repaired sensitivity tested 14421 expressed Monocyte_CD14 genes at donor pseudobulk level. Strict F-series-only Mann-Whitney sensitivity (frail F n=5 vs healthy-old F n=3) produced 0 BH-significant genes and 0 overlap with the original significant rows. Platform-adjusted OLS across primary-contrast donors (n=11; model log2CPM ~ frail + OH platform) produced 1190 BH-significant frail-coefficient genes, including 270 overlapping original rows.

Interpretation for downstream claims: the broad original positive Monocyte_CD14 pseudobulk claim is platform-sensitive and does not survive the strict F-series-only nonparametric sensitivity. Any Monocyte_CD14 transcriptional signal must be framed as exploratory/platform-sensitive rather than a central robust finding; 07-manuscript should soften or remove central positive pseudobulk language.

---

# MJR-005 figure-panel contract repair (05-analysis rerun)
Status: COMPLETE for 05-analysis source-traceability repair. The explicit 09-reviewer blocker MJR-005 was addressed within the 05-analysis boundary by auditing every retained 06-figures panel against existing 05-analysis source result files, source data files, producing scripts, and rendered PNG outputs. No new biological statistics or publication figures were generated in this repair pass.

Verification result: 19/19 retained panels have present source result, source data, source script, and PNG output paths from the existing analysis/figure workspaces; missing source count = 0.

Outputs added/updated:
- figure_data_map.tsv now contains one row per retained Fig1-Fig5 panel and includes panel_status, claim_id, and output_png audit fields.
- statistical_results/mjr005_figure_panel_contract_repair.tsv records the previous 06 panel status, source-file existence checks, retained-panel acceptance, and caveats.
- sensitivity_checks/mjr005_figure_contract_summary.json records machine-readable closure metrics.
- logs/mjr005_figure_contract_repair.log records the executed command output.

Interpretive constraints remain unchanged: composition testing is null after BH correction; Monocyte_CD14 transcriptional evidence is exploratory/platform-sensitive after F-series sensitivity; TCR panels are descriptive donor-level context only; and 06-figures must propagate the repaired accepted statuses into its own panel_data_map.tsv, figure_callout_contract.yaml, captions, and figure_manifest.yaml before reviewer rerun.

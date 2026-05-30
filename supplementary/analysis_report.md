# Analysis Report — MJR-001 repair update
Status: complete after sensitivity-analysis update. The analysis environment reran the F-series-only donor pseudobulk sensitivity analysis from the archived project archive at 2026-05-29T14:10:42+00:00. Eligible F-series contrast was frail n=5 vs healthy-old/non-frail n=3. Monocyte_CD14 did not retain BH-adjusted significance in the sensitivity (retained q<0.05 rows=0, minimum retained q=0.192379). Across all tested cell types, retained q<0.05 rows=0. Therefore the original mixed-series Monocyte_CD14 positive pseudobulk signal should be treated as exploratory and not as a robust central claim.

See `repair_report.md`, `logs/mjr001_platform_sensitivity.log`, and `sensitivity_checks/mjr001_closure_summary.json`.

---

# Analysis Report
Status: complete with interpretive caveats. M0 passed using corrected 02 labels. Donor is replicate; no cell-level pseudoreplication. QC-pass cells 134683; donors 17; common genes 33694. Primary contrast frail n=5 vs healthy-old n=6. Composition tests 9, BH significant 0. Pseudobulk rows 4500, BH significant reported rows 500. TCR donors parsed 17. Claims must remain conservative because of small n and OH/F platform imbalance.

TCR repair note: TCR contig files use separate GEO sample IDs from scRNA matrices; donor-ID mapping was applied and donor-level TCR summaries were regenerated.

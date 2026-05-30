# Figure integrity report

Status: CONDITIONAL

Generated PNG figures: 5
Rendered panels: 19 (11 accepted; 8 conditional; 0 omitted)

Traceability:
- 05-analysis handoff, analysis manifest/report, figure_data_map.tsv, mechanism_map.tsv, and evidence_chain.tsv were checked for availability.
- 03-architect figure_plan.tsv and 05-analysis figure_data_map.tsv were both read and copied into the 06 workspace.
- Every rendered panel maps to an existing 05-analysis source file/result/script path in the manifest; missing verified paths: 0.
- PNGs were regenerated from figure_scripts/render_figures.py during this stage execution.
- Figure-only reshaping/ordering is limited to plotting preparation inside the render script.

Unsupported/limited claims:
- 05-analysis status is conditional after MJR-001 repair. The Monocyte_CD14 mixed-series pseudobulk signal must not be presented as a robust central positive finding.
- F-series-only sensitivity found no retained BH-significant Monocyte_CD14 pseudobulk rows; Fig3 presents this as sensitivity/repair context.
- Current evidence supports five main figures, not the aspirational seven or more main figures.

Omitted panels:
- No decorative placeholders were produced.
- No rendered panel has a missing source path. Panels not identically present in both 03 and 05 contracts are marked conditional and listed in figure_discrepancy_report.md.

Upstream repair route if stricter figure-plan compliance is required:
- Return to 03-architect/05-analysis to reconcile figure_plan.tsv with figure_data_map.tsv and supply additional mapped analysis outputs if seven or more main figures are mandatory.

Updated: 2026-05-29T15:15:10Z

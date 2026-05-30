# Figure integrity report

Status: PASS WITH LIMITATIONS

Generated PNG figures: 5
Rendered panels: 19 (11 accepted; 8 evidence-limited; 0 omitted rendered placeholders)

Traceability:
- analysis workflow transfer, analysis manifest/report, figure_data_map.tsv, mechanism_map.tsv, and evidence_chain.tsv were checked for availability.
- 03-architect figure_plan.tsv and analysis workflow figure_data_map.tsv were both read.
- Every rendered panel maps to an existing analysis workflow source file and an existing source script/result path; no missing source paths were detected during manifest verification.
- The reproducible render script is `figure_scripts/render_figures.py`; figure-only reshaping/ordering was performed inside that script for plotting.

Unsupported/limited claims:
- analysis workflow status is limited after MJR-001 repair. The Monocyte_CD14 mixed-series pseudobulk signal must not be presented as a robust central positive finding.
- F-series-only sensitivity found no retained BH-significant Monocyte_CD14 pseudobulk rows; Fig3 presents this as sensitivity/repair context.
- Current evidence supports five main figures, not the aspirational seven or more main figures.

Omitted panels:
- No decorative placeholders were produced.
- Panels absent from either the 03 or 05 contracts but plotted from concrete 05 files are marked accepted with explicit limitations.

source repair route if stricter figure-plan compliance is required:
- Return to 03-architect/analysis workflow to reconcile figure_plan.tsv with figure_data_map.tsv and supply additional mapped analysis results if seven or more main figures are mandatory.

Updated: 2026-05-29T14:13:49Z

# Figure discrepancy report

Status: conditional. 03-architect figure_plan.tsv and 05-analysis figure_data_map.tsv were both read. Differences are documented rather than silently redesigned.

Panels in 03 plan not explicitly mapped in 05 figure_data_map.tsv:
- Fig1D
- Fig2B
- Fig2C
- Fig3D
- Fig4B
- Fig4C
- Fig5C
- Fig5D

Panels in 05 figure_data_map.tsv not present in 03 figure_plan.tsv:
- ReviewerRepairMJR-001A
- ReviewerRepairMJR-001B

Repair route: reconcile 03-architect figure_plan.tsv and 05-analysis figure_data_map.tsv if strict panel-contract identity is required. Current PNGs use only existing 05-analysis source files/results.

# GSE157007 frailty PBMC single-cell reanalysis artifacts

This repository contains analysis scripts, derived result tables, figure-generation scripts, figure panels, and reproducibility metadata for the manuscript:

**Donor-aware reanalysis of public single-cell PBMC data identifies constrained frailty-associated immune signals after platform sensitivity analysis**

Author: Zhun Qiu, Northeast Forestry University, Harbin, China.  
ORCID: https://orcid.org/0009-0005-6891-1841  
Contact: qz@nefu.edu.cn

## Data source

The source single-cell sequencing data are publicly available from GEO under accession **GSE157007**. Raw 10x matrices and the processed `.h5ad` object are not committed here because of repository size constraints. This repository provides scripts, manifests, derived tables, quality-control summaries, sensitivity analyses, figure-source artifacts, and final figure panels needed to audit the reported analytical logic.

## Repository contents

- `analysis_scripts/` — single-cell preprocessing, donor-level summaries, differential analysis, TCR parsing/repair, and sensitivity checks.
- `analysis_results/` — derived tables, QC summaries, logs, manifests, and sensitivity-check outputs.
- `figure_scripts/` — scripts used to render figure panels.
- `figures/` — PNG figure panels used in the manuscript submission package.
- `manuscript_statements/` — data/code availability and administrative statements.
- `docs/` — figure manifests, panel data maps, integrity reports, and package manifest.

## Reproducibility notes

The analysis manifest records the exact local command and package versions used in the CS1 run. Key runtime components included Python 3.10, scanpy 1.11.5, anndata 0.11.4, numpy 2.2.6, pandas 2.3.3, scipy 1.15.3, statsmodels 0.14.6, matplotlib 3.10.9, and seaborn 0.13.2.

Because this repository intentionally excludes large raw/intermediate matrix files, a full re-run requires downloading GSE157007 source data from GEO and reconstructing the local input layout described in `analysis_results/analysis_manifest.yaml` and `analysis_results/processed_data/input_manifest.tsv`.

## Repository URL

https://github.com/qz-bbb/gse157007-frailty-pbmc-reanalysis

## License

No explicit reuse license has been selected yet. Contact the author before reusing unpublished manuscript-specific derived materials.

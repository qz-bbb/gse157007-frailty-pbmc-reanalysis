# Donor-aware reanalysis of public single-cell PBMC data identifies constrained frailty-associated immune signals after platform sensitivity analysis

## Abstract

**Background:** Frailty and chronological aging overlap clinically but are not identical biological states. Recent single-cell immune-aging studies across human longevity cohorts, animal immune tissues, primate immune organs, and frailty-focused T-cell work show that immune aging is context-, tissue-, and cell-state-dependent. GSE157007 provides public human PBMC single-cell RNA and TCR sequencing data across cord blood, young adult, healthy-old/non-frail, and frail donors, enabling a transparent donor-aware reanalysis of immune aging and frailty.

**Methods:** We reaudited local raw GSE157007 provenance, applied corrected donor-group labels, harmonized raw scRNA-seq matrices across 17 donors, annotated major PBMC populations using canonical marker support, and treated donor as the replicate unit. The primary contrast was frail donors (n = 5) versus healthy-old/non-frail donors (n = 6). We tested donor-level cell-type composition, performed exploratory donor pseudobulk transcriptional comparisons, parsed donor-level scTCR contig files, and conducted platform/series sensitivity analyses including an F-series-only contrast.

**Results:** The analysis retained 134,683 QC-pass cells from 17 donors and 33,694 common genes. Corrected groups comprised 5 frail, 6 healthy-old/non-frail, 3 young-adult, and 3 cord-blood donors. Nine major PBMC populations had requested marker support. Donor-level composition testing across nine cell types found no Benjamini–Hochberg (BH)-significant frail versus healthy-old/non-frail differences. The original mixed-series pseudobulk screen produced Benjamini–Hochberg (BH)-significant Monocyte_CD14 rows, but targeted NEAT1/MALAT1 monocyte checks were not Benjamini–Hochberg (BH)-significant and an F-series-only donor pseudobulk sensitivity found no retained Benjamini–Hochberg (BH)-significant rows in Monocyte_CD14 or any tested cell type. TCR files were parsed for all 17 donors and are reported as descriptive repertoire context.

**Conclusions:** This reproducible reanalysis supports a transparent, label-corrected view of GSE157007 and shows that frailty-associated PBMC signals are strongly constrained by small donor numbers and platform/series imbalance. The data do not support a robust platform-independent monocyte transcriptional claim; instead, they support cautious reporting of composition-null results, exploratory transcriptional patterns, and descriptive TCR context for future validation.

**Keywords:** frailty; immune aging; PBMC; single-cell RNA sequencing; TCR repertoire; donor-aware analysis; GSE157007

## Introduction

Frailty in later life is related to aging but is not reducible to chronological age. Immune remodeling is one plausible component of this divergence because peripheral blood immune populations integrate inflammatory exposure, antigen experience, and age-associated changes. Public single-cell datasets make it possible to revisit such questions with transparent computational workflows, but they also require careful attention to donor labels, sample provenance, biological replication, and platform confounding.

GSE157007 is a public human dataset generated for multidimensional single-cell analysis of peripheral blood immune aging and frailty [1]. The original study included cord blood, young-adult PBMC, healthy-old/non-frail PBMC, and frail PBMC donors, with scRNA-seq and scTCR-seq components and reported frailty-specific immune features [1]. Recent immune-aging and frailty-immunology studies provide broader context for single-cell immune aging and frailty-associated immune dysfunction [2-5]. However, the appropriate contribution of the present study is not rediscovery of the dataset or replacement of the original analysis. The contribution is a analysis-auditable, donor-aware reanalysis using corrected labels, raw-data provenance checks, major-cell-type composition tests, exploratory pseudobulk transcriptional analyses, TCR repertoire summaries, and explicit sensitivity analysis for platform/series imbalance.

The central question was whether frail older donors in GSE157007 show immune-cell composition, transcriptional-state, or TCR repertoire differences compared with healthy-old/non-frail donors after enforcing corrected labels and donor-level inference. We prespecified a conservative claim ceiling: candidate frailty-associated immune remodeling could be reported only where supported by donor-aware analysis; causal mechanisms, clinical biomarkers, and broad population claims were excluded.

## Methods

### Data source and provenance audit

This study used GSE157007, a public human PBMC/cord-blood single-cell dataset associated with the primary source paper [1]. We verified the public GEO accession GSE157007, including the raw archive, GEO SOFT and series-matrix metadata, and a source-access receipt retained in supporting provenance records. The public-data provenance audit recorded the source archive size and checksum, with reproducibility details retained in the accompanying analysis artifacts. Input provenance was checked against GEO files and accession metadata.

### Donor labels and analysis groups

Corrected donor labels were taken from the curated donor-group assignment table and used consistently for all donor-level analyses. The analysis groups were frail (F002, F004, F006, F007, F008; n = 5), healthy-old/non-frail (F020, F021, F023, OH14, OH15, OH17; n = 6), young adult (F012, F013, F014; n = 3), and cord blood (F016, F017, F024; n = 3). The original Supplementary Tables 1 and 2 were not present as separate machine-readable standalone public supplements, so label provenance and this caveat were retained throughout.

### scRNA-seq loading, harmonization, and QC

Raw 10X matrices were loaded by donor, harmonized to common genes, and QC filtered by the reproducible analysis workflow. The quality-controlled common-gene dataset contains the QC-annotated common-gene object used for downstream summaries. Donor-level QC summaries were written to `qc/qc_by_sample.tsv` and `qc/cell_qc_summary.tsv`. No new statistical analysis was performed during manuscript drafting.

### Cell-type annotation

Major PBMC populations were annotated using marker support summarized in `statistical_results/celltype_marker_support.tsv` and `statistical_results/annotation_confidence.tsv`. The annotated classes were B, CD4_T, CD8_T, NK, Monocyte_CD14, Monocyte_FCGR3A, DC, Platelet, and Plasma. Requested canonical markers were present for each of these major populations.

### Donor-level composition testing

The primary biological replicate was the donor. For frail versus healthy-old/non-frail composition analysis, analysis workflow used donor fractions and Mann–Whitney U tests with BH adjustment across the nine tested cell types. The result was `statistical_results/composition_differential.tsv`. Cell-level pseudoreplication was not used for inferential composition claims.

### Donor pseudobulk transcriptional analysis and sensitivity testing

Exploratory cell-type-specific pseudobulk comparisons were generated by donor for the frail versus healthy-old/non-frail contrast and written to `statistical_results/de_by_celltype.tsv`. The analysis workflow report identified the initial mixed-series Monocyte_CD14 positive rows as exploratory because frail donors were all F-series while healthy-old/non-frail controls included both F-series and OH-series donors. A repair analysis then performed an F-series-only donor pseudobulk sensitivity with frail n = 5 and healthy-old/non-frail F-series n = 3. A platform covariate model was not identifiable because frailty status was nested within F-series. The F-series-only sensitivity used donor pseudobulk log2(CPM+1), two-sided Mann–Whitney U tests per gene, and BH adjustment within cell type, as documented in `repair_report.md`.

### TCR repertoire parsing

scTCR contig annotation files were parsed for all 17 donors after applying donor-ID mapping between TCR GEO sample identifiers and scRNA donor identifiers. Donor-level summaries included total contigs, productive contigs, unique clonotypes, expanded clonotypes, Shannon diversity, and Simpson diversity. Because scRNA/TCR barcode linkage was limited, TCR results are presented as descriptive donor-level repertoire context rather than inferential mechanistic evidence.

### Figures

Only accepted PNG files from figure workflow were used. Figure-linked claims follow the figure workflow callout contract.

- Fig1: Fig1_dataset_qc.png (submitted as a separate PNG figure file)
- Fig2: Fig2_annotation_composition.png (submitted as a separate PNG figure file)
- Fig3: Fig3_transcription_sensitivity.png (submitted as a separate PNG figure file)
- Fig4: Fig4_tcr_diversity.png (submitted as a separate PNG figure file)
- Fig5: Fig5_context_sensitivity.png (submitted as a separate PNG figure file)

## Results

### Provenance and corrected labels support a reproducible but small donor-aware reanalysis

The analysis retained 134,683 QC-pass cells from 17 donors with 33,694 common genes. Group-level QC-pass cell counts were: 

```
       analysis_group  donors  cells  median_cells
           cord_blood       3  21826        7808.0
                frail       5  22136        3702.0
healthy_old_non_frail       6  58231       10088.5
          young_adult       3  32490        7006.0
```

The corrected donor label lock supported the primary frail versus healthy-old/non-frail contrast, with frail n = 5 and healthy-old/non-frail n = 6. The full dataset also included young-adult and cord-blood donors for context. These provenance and label results are shown in Figure 1.

**Figure 1 callout:** Dataset provenance, corrected donor groups, and QC summaries are shown in Fig. 1.

### Major PBMC populations had marker support, but composition differences were not Benjamini–Hochberg (BH)-significant

Nine major PBMC populations had requested marker support in the analysis result: B, CD4_T, CD8_T, NK, Monocyte_CD14, Monocyte_FCGR3A, DC, Platelet, and Plasma. Donor-level composition tests across these nine cell types did not identify any Benjamini–Hochberg (BH)-significant frail versus healthy-old/non-frail differences. The observed mean fraction differences were therefore interpreted as descriptive/candidate effects only: B: mean difference 0.036, BH q=0.643; CD4_T: mean difference -0.043, BH q=1.000; CD8_T: mean difference -0.054, BH q=0.643; DC: mean difference 0.017, BH q=0.643; Monocyte_CD14: mean difference 0.055, BH q=0.643; Monocyte_FCGR3A: mean difference -0.005, BH q=0.643; NK: mean difference 0.007, BH q=1.000; Plasma: mean difference 0.003, BH q=0.643; Platelet: mean difference -0.016, BH q=0.852.

**Figure 2 callout:** Cell-type marker support and donor-level composition contrasts are shown in Fig. 2.

### Mixed-series Monocyte_CD14 pseudobulk signal was exploratory and not robust to F-series-only sensitivity

The mixed-series donor pseudobulk screen generated 4,500 rows across nine cell types. In that initial screen, 500 Monocyte_CD14 rows had BH-adjusted q < 0.05, with a minimum q-value of 0.046381. However, this signal was not treated as a robust central finding because the primary contrast was vulnerable to platform/series imbalance: all frail donors were F-series, whereas healthy-old/non-frail donors included three F-series and three OH-series donors.

Targeted monocyte NEAT1/MALAT1 checks did not support a Benjamini–Hochberg (BH)-significant frailty-associated increase in this reanalysis: Monocyte_CD14 NEAT1: effect=-0.392, BH q=0.355; Monocyte_CD14 MALAT1: effect=-0.250, BH q=0.329; Monocyte_FCGR3A NEAT1: effect=-0.172, BH q=0.635; Monocyte_FCGR3A MALAT1: effect=0.003, BH q=0.762. The MJR-001 F-series-only sensitivity found no retained Benjamini–Hochberg (BH)-significant rows in Monocyte_CD14 or any other tested cell type: B: retained q < 0.05 rows=0, min q=0.235; CD4_T: retained q < 0.05 rows=0, min q=0.297; CD8_T: retained q < 0.05 rows=0, min q=0.323; DC: retained q < 0.05 rows=0, min q=0.268; Monocyte_CD14: retained q < 0.05 rows=0, min q=0.192; Monocyte_FCGR3A: retained q < 0.05 rows=0, min q=0.541; NK: retained q < 0.05 rows=0, min q=0.229; Plasma: retained q < 0.05 rows=0, min q=0.643; Platelet: retained q < 0.05 rows=0, min q=0.864. Therefore, the manuscript does not claim a robust platform-independent Monocyte_CD14 transcriptional program.

**Figure 3 callout:** Exploratory mixed-series transcriptional results, targeted NEAT1/MALAT1 checks, and the F-series-only sensitivity are shown in Fig. 3.

### TCR repertoire summaries provide descriptive donor-level context

TCR contig files were parsed for all 17 donors. Group-level descriptive summaries were:

```
                group  donors  productive_contigs  unique_clonotypes  median_shannon  median_expanded_clonotypes
           cord_blood       3               19564              10624        8.114661                      2773.0
                frail       5               14943               6286        6.387083                       607.0
healthy_old_non_frail       6               51236              18407        6.945241                      2115.5
          young_adult       3               21657              11725        8.001757                      2480.0
```

These summaries show donor-level repertoire measurements but were not used to infer a statistically validated frailty-specific TCR mechanism. The TCR analysis is most appropriately viewed as orthogonal descriptive context supporting reproducible reuse of the multimodal dataset.

**Figure 4 callout:** Donor-level TCR diversity, clonotype expansion summaries, and parse status are shown in Fig. 4.

### Age-context and platform summaries constrain interpretation

Age/context and platform summaries showed the central limitation of the frailty contrast. Donor counts and median total cells differed by group, and platform/series distribution was imbalanced:

```
       analysis_group platform_series  donors  median_total_cells
           cord_blood               F       3              7808.0
                frail               F       5              3702.0
healthy_old_non_frail               F       3              6934.0
healthy_old_non_frail              OH       3             12260.0
          young_adult               F       3              7006.0
```

The F-series-only composition sensitivity also found no Benjamini–Hochberg (BH)-significant composition differences. These results support the conclusion that apparent frailty-associated signals in this dataset must be interpreted cautiously and validated independently.

**Figure 5 callout:** Age/context donor counts, platform/series balance, and F-series composition sensitivity are shown in Fig. 5.

## Discussion

This donor-aware reanalysis of GSE157007 provides a reproducible, corrected-label view of public PBMC single-cell data in aging and frailty. The main contribution is not a new strong biomarker claim; rather, it is the transparent evaluation of what remains supportable after provenance auditing, donor-level inference, and platform/series sensitivity analysis.

The composition analysis did not identify Benjamini–Hochberg (BH)-significant differences across the nine tested major PBMC populations. This result does not exclude frailty-associated immune remodeling in larger or independently collected cohorts, but it prevents strong composition claims from this dataset alone. Similarly, the original mixed-series Monocyte_CD14 pseudobulk signal cannot be advanced as a robust frailty program because the F-series-only sensitivity found no retained Benjamini–Hochberg (BH)-significant rows. The targeted NEAT1/MALAT1 checks were also not Benjamini–Hochberg (BH)-significant in this reanalysis, which is important because the original dataset paper emphasized a frailty-specific long-noncoding-RNA-high monocyte signal [1]. The appropriate interpretation is therefore refinement and constraint of prior observations, not contradiction by a better-powered independent cohort.

The TCR analysis demonstrates that donor-level repertoire summaries can be regenerated from the public files, including productive contig and clonotype diversity metrics. However, without robust barcode-level linkage and with small group sizes, these summaries should be used as descriptive context rather than mechanistic proof of frailty-associated T-cell repertoire remodeling.

This study also illustrates a broader lesson for public single-cell reuse. Small human donor cohorts can contain strong technical or sampling structure. When all cases in a contrast are nested within one series or platform, positive transcriptional findings require sensitivity analyses before biological interpretation. In this analysis, the sensitivity analysis materially changed the claim: the manuscript moves from a candidate positive Monocyte_CD14 program to a more conservative conclusion that the current evidence is exploratory and platform-constrained.

## Limitations

1. The primary frail versus healthy-old/non-frail comparison used only 5 frail and 6 healthy-old/non-frail donors.
2. Frail donors were all F-series, while healthy-old/non-frail donors included both F-series and OH-series donors; a platform covariate model was not identifiable for frailty because of nesting.
3. The F-series-only sensitivity reduced the healthy-old/non-frail control group to n = 3, limiting power.
4. Composition tests found no Benjamini–Hochberg (BH)-significant differences.
5. The mixed-series Monocyte_CD14 pseudobulk signal was not robust to F-series-only sensitivity and is exploratory only.
6. TCR results are descriptive because robust scRNA/TCR linkage was not available in the downstream claim set.
7. Original Supplementary Tables 1 and 2 were not present as separate machine-readable artifacts in this analysis; corrected labels were derived from source evidence and documented with this caveat.
8. No new wet-lab validation or independent cohort replication was performed.

## Conclusions

A corrected-label, donor-aware reanalysis of GSE157007 supports reproducible reuse of public single-cell PBMC aging/frailty data but does not support a robust platform-independent frailty-associated Monocyte_CD14 transcriptional claim. The strongest evidence is a transparent provenance/QC/annotation workflow, null donor-level composition testing after BH adjustment, exploratory transcriptional results constrained by F-series-only sensitivity, and descriptive TCR repertoire context. Future studies should validate frailty-associated immune remodeling in larger, prospectively balanced cohorts.

## Data availability



The source single-cell sequencing data are publicly available from GEO under accession GSE157007. Derived analysis tables, figure-source artifacts, reproducibility metadata, and analysis scripts from this reanalysis are available in the public GitHub repository: https://github.com/qz-bbb/gse157007-frailty-pbmc-reanalysis.

## Ethics approval and consent to participate


This study is a secondary analysis of publicly available data from GEO accession GSE157007. No new human participants were recruited and no new human biospecimens were collected for this reanalysis. The source study by Luo et al. (Nature Aging, 2022; DOI: 10.1038/s43587-022-00198-9) reported that the study was approved by the institutional review board at the First Affiliated Hospital, Jinan University (approval KY-2020-027), that participants were enrolled at Guangzhou First People's Hospital and the First Affiliated Hospital, Jinan University, and that written consents were obtained. This reanalysis used only public, de-identified repository data and therefore required no additional ethics approval or consent.

## Consent for publication


Not applicable for this secondary public-data reanalysis.

## Competing interests


The author declares no competing interests.

## Funding


This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## Authors' contributions


Zhun Qiu (邱准) is the sole author and was responsible for conceptualization, methodology, formal analysis, investigation, data curation, visualization, writing – original draft, writing – review and editing, and project administration.

## Acknowledgements


None.

## References

1. Luo OJ, Lei W, Zhu G, Ren Z, Xu Y, Xiao C, et al. Multidimensional single-cell analysis of human peripheral blood reveals characteristic features of the immune system landscape in aging and frailty. Nature Aging. 2022;2(4):348-364. doi:10.1038/s43587-022-00198-9.

2. Zhang Z, Su H, Fu S, Ji F, Liang L, Song W, et al. The Immune Cell Atlas of "Longevity Molecular Tag": Identification of Principal Immune Cell Subsets and Their Underlying Molecular Regulatory Mechanisms. Aging Cell. 2026;25(3):e70431. doi:10.1111/acel.70431.

3. Philip CS, Filippov I, Haljasorg U, Peterson P. Single-Cell Profiling of Splenic Immune Ageing and Chronic Stress Adaptations in Mice With Natural Microbiota. European Journal of Immunology. 2026;56(3):e70170. doi:10.1002/eji.70170.

4. Wang S, Zhu Z, Yang H, Yan Y, Ran L, Yang N, et al. Single-Cell Profiling Across Immune Tissues and Organs Reveals Immunosenescence Signatures in Male Rhesus Monkeys. Advanced Science. 2026;13(17):e14353. doi:10.1002/advs.202514353.

5. Li H, Li X, Chen H, Yang Z, Zheng M, Xue Y, et al. Impaired function of Vgamma9Vdelta2 T cells in frail elderly. Immunity & Ageing. 2026;23(1):7. doi:10.1186/s12979-026-00558-8.

## Figure legends


Figure 1. Dataset provenance, corrected donor groups, and quality-control summaries. (A) Verified GEO source-file provenance and corrected group assignments. (B) Donor counts by corrected analysis group and platform series. (C) QC-pass cells by donor. (D) Donor-level QC diagnostic showing median genes and mitochondrial percentage, with point size proportional to QC-pass cells.

Figure 2. Cell-type annotation support and donor-level composition contrast. (A) Requested canonical markers present for each annotated PBMC cell type. (B) Mean donor fraction differences for frail versus healthy-old/non-frail donors with BH-adjusted q-values from analysis workflow. No composition test was Benjamini–Hochberg (BH)-significant.

Figure 3. Exploratory transcriptional results and platform-aware sensitivity. (A) Mixed-series Monocyte_CD14 donor pseudobulk results are shown as exploratory. (B) Targeted NEAT1/MALAT1 monocyte checks with BH-adjusted q-values. (C) F-series-only sensitivity showed no retained BH-adjusted q < 0.05 rows in Monocyte_CD14 or other tested cell types; central positive pseudobulk claims must therefore be softened or removed.

Figure 4. Descriptive TCR diversity summaries. (A) Donor-level Shannon diversity by analysis group. (B) Unique versus expanded clonotypes by donor. (C) TCR parse status summary. These panels are descriptive donor-level summaries from analysis workflow.

Figure 5. Context and platform sensitivity summaries. (A) Donor counts and median total cells across analysis groups. (B) Platform/series distribution by analysis group. (C) F-series-only composition sensitivity effect estimates with BH-adjusted q-values.



## Code availability

Analysis and figure-generation scripts used for this reanalysis are available in the public GitHub repository: https://github.com/qz-bbb/gse157007-frailty-pbmc-reanalysis.

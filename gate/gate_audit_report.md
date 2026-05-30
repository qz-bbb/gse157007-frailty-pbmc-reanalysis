# 10-gate final audit

Generated UTC: 2026-05-30T10:31:10Z

Final verdict: pass
Delivery ready: true

Reader PDF quality: single-source ReportLab build; no legacy PDF append; page_count=13; forbidden_text_hits=0; title_occurrence_count=1.

Checks:
- PACKAGE_MANIFEST: PASS — entries=37 mismatches=0 []
- INTERNAL_LEAKAGE: PASS — leaks=0 sample=[]
- READER_PDF: PASS — exists=True size=2016453 manifest_size=2016453 pages=13 text_layer_verified=True quality=single_source_no_recursive_append title_occurrences=1 forbidden_hits=0
- REVIEWER_CLOSURE: PASS — status=pass verdict=pass_minor_accepted_for_gate fatal=0 major=0
- PIPELINE_VALIDATOR: PASS — VALIDATION PASSED: No issues found.
- PUBLIC_REPOSITORY: PASS — https://github.com/qz-bbb/gse157007-frailty-pbmc-reanalysis

#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from datetime import datetime, timezone
from pathlib import Path
PROJECT_ID="cs1-paper-factory"; RUN_ID="run_20260529_095727"
RUN_DIR=Path("/mnt/d/cs1/runs/cs1-paper-factory/run_20260529_095727")
WS=RUN_DIR/"agents/05-analysis/workspace"
PANEL_MAP=RUN_DIR/"agents/06-figures/workspace/panel_data_map.tsv"
FIGURE_CONTRACT=RUN_DIR/"agents/06-figures/workspace/figure_callout_contract.yaml"
FIGURE_MANIFEST=RUN_DIR/"agents/06-figures/workspace/figure_manifest.yaml"
def win_path(path: Path)->str:
    s=str(path); return "D:/"+s[len("/mnt/d/"):] if s.startswith("/mnt/d/") else s
def rel_to_ws(path_text: str)->Path:
    marker="/agents/05-analysis/workspace/"; normalized=path_text.replace("\\","/")
    if marker in normalized: return WS/normalized.split(marker,1)[1]
    if normalized.startswith("D:/"): return Path("/mnt/d")/normalized[3:]
    return WS/normalized
def out_path(path_text: str)->Path:
    normalized=path_text.replace("\\","/")
    if normalized.startswith("D:/"): return Path("/mnt/d")/normalized[3:]
    return RUN_DIR/normalized
def rel(p: Path)->str:
    try: return str(p.relative_to(WS)).replace("\\","/")
    except ValueError: return win_path(p)
def caveat_for(title: str, stats: str)->str:
    if "TCR" in title or "TCR" in stats: return "descriptive donor-level TCR context only"
    if "pseudobulk" in title.lower() or "F-series" in title or "Monocyte" in stats: return "exploratory/platform-sensitive; strict F-series Monocyte_CD14 sensitivity has no retained BH-significant rows"
    if "Composition" in title or "composition" in stats or "fractions" in stats: return "donor-level small-n composition inference; no BH-significant primary composition tests"
    if "QC" in title or "provenance" in stats or "labels" in stats: return "provenance/QC support; original machine-readable Supplementary Tables 1/2 absent"
    return "small donor n and OH/F platform imbalance limit interpretation"
def main()->int:
    with PANEL_MAP.open(newline="",encoding="utf-8") as handle: rows=list(csv.DictReader(handle, delimiter="	"))
    required_inputs=[PANEL_MAP,FIGURE_CONTRACT,FIGURE_MANIFEST]
    missing_required=[win_path(p) for p in required_inputs if not p.exists()]
    if missing_required: raise SystemExit(f"Missing required 06 contract inputs: {missing_required}")
    rich_rows=[]; figure_rows=[]; missing_sources=[]
    for row in rows:
        sr=rel_to_ws(row["source_result"]); sd=rel_to_ws(row["source_data"]); sp=rel_to_ws(row["source_script_or_analysis_manifest_entry"]); op=out_path(row["output_png"])
        checks={"source_result":sr.exists(),"source_data":sd.exists(),"source_script":sp.exists(),"output_png":op.exists()}
        status="accepted" if all(checks.values()) else "blocked_missing_source"
        if status!="accepted":
            paths={"source_result":sr,"source_data":sd,"source_script":sp,"output_png":op}
            for label, ok in checks.items():
                if not ok: missing_sources.append(f"{row['figure_id']}{row['panel_id']} {label}: {win_path(paths[label])}")
        stats=row.get("statistics_used",""); title=row.get("panel_title",""); caveat=caveat_for(title, stats)
        base={"figure_id":row["figure_id"],"panel_id":row["panel_id"],"source_result":rel(sr),"source_script":rel(sp),"data_file":rel(sd),"statistic":stats,"caveat":caveat,"panel_status":status,"claim_id":row.get("claim_id_or_evidence_id",""),"output_png":win_path(op)}
        figure_rows.append(base)
        rich=dict(base); rich.update({"panel_title":title,"plot_type":row.get("plot_type",""),"exists_source_result":str(checks["source_result"]),"exists_data_file":str(checks["source_data"]),"exists_source_script":str(checks["source_script"]),"exists_output_png":str(checks["output_png"]),"previous_06_status":row.get("status",""),"repair_action":"accepted_traceable_retained_panel" if status=="accepted" else "route_to_06_remove_or_replan"})
        rich_rows.append(rich)
    figure_fields=["figure_id","panel_id","source_result","source_script","data_file","statistic","caveat","panel_status","claim_id","output_png"]
    with (WS/"figure_data_map.tsv").open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h, fieldnames=figure_fields, delimiter="	", lineterminator=chr(10)); w.writeheader(); w.writerows([{k:r[k] for k in figure_fields} for r in figure_rows])
    repair_fields=["figure_id","panel_id","panel_title","panel_status","previous_06_status","source_result","data_file","source_script","output_png","exists_source_result","exists_data_file","exists_source_script","exists_output_png","claim_id","plot_type","statistic","caveat","repair_action"]
    repair_path=WS/"statistical_results/mjr005_figure_panel_contract_repair.tsv"
    with repair_path.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h, fieldnames=repair_fields, delimiter="	", lineterminator=chr(10)); w.writeheader(); w.writerows([{k:r[k] for k in repair_fields} for r in rich_rows])
    summary={"agent":"05-analysis","project_id":PROJECT_ID,"run_id":RUN_ID,"issue_id":"MJR-005","timestamp_utc":datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00","Z"),"panels_checked":len(rich_rows),"panels_accepted_by_05_source_audit":sum(r["panel_status"]=="accepted" for r in rich_rows),"panels_blocked_missing_source":sum(r["panel_status"]!="accepted" for r in rich_rows),"missing_sources":missing_sources,"figure_data_map":win_path(WS/"figure_data_map.tsv"),"repair_table":win_path(repair_path),"inputs_checked":[win_path(p) for p in required_inputs],"interpretive_constraints":["No publication figures were generated by 05-analysis.","Panels marked accepted are accepted only for 05-analysis source traceability; 06-figures must rerun panel_data_map/figure_manifest/callout status propagation.","Composition claims remain null after BH correction.","Monocyte_CD14 transcriptional claims remain exploratory/platform-sensitive after strict F-series sensitivity.","TCR panels remain descriptive donor-level context."]}
    summary_path=WS/"sensitivity_checks/mjr005_figure_contract_summary.json"; summary_path.write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps(summary,indent=2)); return 0 if not missing_sources else 2
if __name__=="__main__": raise SystemExit(main())

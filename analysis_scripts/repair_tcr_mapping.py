import tarfile,gzip,json,datetime
from pathlib import Path
import pandas as pd, numpy as np
RUN=Path('D:/cs1/runs/cs1-paper-factory/run_20260529_095727'); WS=RUN/'agents/05-analysis/workspace'; TAR=Path('D:/cs2/data/raw/ds_geo_gse157007/GSE157007_RAW.tar')
smap=pd.read_csv(WS/'processed_data/sample_file_map.tsv',sep='\t')
with tarfile.open(TAR) as tf:
    names=tf.getnames()
    tcr_names=[n for n in names if 'filtered_contig_annotations' in n]
    fixed=[]; tsum=[]; exp=[]
    for _,r in smap.iterrows():
        donor=str(r.donor_id)
        tcr=[n for n in tcr_names if ('_'+donor+'_') in n or ('_'+donor+'_filtered') in n]
        rec=r.to_dict(); rec['tcr_files']=';'.join(tcr); fixed.append(rec)
        if not tcr:
            tsum.append(dict(donor_id=donor,sample_id=r.sample_id,group=r.analysis_group_confirmed,status='missing_tcr_file'))
            continue
        df=pd.read_csv(gzip.GzipFile(fileobj=tf.extractfile(tcr[0])))
        prod=df[df.productive.astype(str).str.lower().isin(['true','t','1'])] if 'productive' in df.columns else df
        col='raw_clonotype_id' if 'raw_clonotype_id' in prod.columns else ('clonotype_id' if 'clonotype_id' in prod.columns else None)
        vc=prod[col].dropna().astype(str).value_counts() if col else pd.Series(dtype=int)
        fq=vc/vc.sum() if vc.sum() else pd.Series(dtype=float)
        for cid,n in vc.head(50).items(): exp.append(dict(donor_id=donor,group=r.analysis_group_confirmed,clonotype_id=cid,count=int(n),frequency=float(n/vc.sum())))
        tsum.append(dict(donor_id=donor,sample_id=r.sample_id,group=r.analysis_group_confirmed,status='parsed',contigs_total=len(df),productive_contigs=len(prod),unique_clonotypes=int(len(vc)),expanded_clonotypes=int((vc>=2).sum()),proportion_productive_in_expanded_clonotypes=float(vc[vc>=2].sum()/vc.sum()) if vc.sum() else np.nan,shannon=float(-(fq*np.log(fq)).sum()) if len(fq) else np.nan,simpson=float(1-(fq**2).sum()) if len(fq) else np.nan,source_file=tcr[0]))
fixed=pd.DataFrame(fixed); fixed.to_csv(WS/'processed_data/sample_file_map.tsv',sep='\t',index=False)
tcr=pd.DataFrame(tsum)
tcr.to_csv(WS/'statistical_results/tcr_parse_summary.tsv',sep='\t',index=False)
tcr.to_csv(WS/'statistical_results/clonotype_diversity_by_donor.tsv',sep='\t',index=False)
pd.DataFrame(exp).to_csv(WS/'statistical_results/expanded_clonotypes.tsv',sep='\t',index=False)
pd.DataFrame([dict(donor_id=d,status='not_integrated',reason='donor-level TCR summaries repaired by donor-ID mapping; barcode linkage not used for core claims') for d in fixed.donor_id]).to_csv(WS/'statistical_results/tcr_scrna_linkage_report.tsv',sep='\t',index=False)
# update report and handoff metrics
n=int((tcr.status=='parsed').sum())
hr=WS/'handoff.json'; hand=json.loads(hr.read_text()); hand['key_metrics']['tcr_donors_parsed']=n; hand['timestamp']=datetime.datetime.now(datetime.timezone.utc).isoformat(); hand.setdefault('repair_notes',[]).append('TCR mapping repaired by donor_id because TCR GEO sample IDs differ from scRNA GEO sample IDs. Donor-level summaries regenerated.') ; hr.write_text(json.dumps(hand,indent=2))
rep=WS/'analysis_report.md'; txt=rep.read_text(); import re; txt=re.sub(r'TCR donors parsed \d+', f'TCR donors parsed {n}', txt); txt+='\nTCR repair note: TCR contig files use separate GEO sample IDs from scRNA matrices; donor-ID mapping was applied and donor-level TCR summaries were regenerated.\n'; rep.write_text(txt)
print('repaired TCR parsed donors',n)

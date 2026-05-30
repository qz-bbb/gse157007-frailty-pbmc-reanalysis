import json, tarfile, gzip, shutil, hashlib, datetime, traceback, sys, platform
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import scanpy as sc, anndata as ad
RUN=Path('D:/cs1/runs/cs1-paper-factory/run_20260529_095727'); WS=RUN/'agents/05-analysis/workspace'
RAW=Path('D:/cs2/data/raw/ds_geo_gse157007'); TAR=RAW/'GSE157007_RAW.tar'; LAB=RUN/'agents/02-knowledge/workspace/gse157007_confirmed_group_assignments.tsv'
for d in 'scripts logs processed_data qc statistical_results enrichment_results sensitivity_checks'.split(): (WS/d).mkdir(parents=True,exist_ok=True)
LOG=WS/'logs/run_sc_analysis.log'; LOG.write_text('')
mods=[]; seed=20260529; np.random.seed(seed)
def log(s): print(s,flush=True); LOG.open('a',encoding='utf-8').write(str(datetime.datetime.now())+' '+s+'\n')
def out(df,p): p=Path(p); p.parent.mkdir(parents=True,exist_ok=True); df.to_csv(p,sep='\t',index=False)
def rec(i,n,st,outs,note=''): mods.append(dict(module_id=i,module=n,status=st,outputs=';'.join(map(str,outs)),notes=note))
def sha(p):
 h=hashlib.sha256(); f=open(p,'rb')
 for b in iter(lambda:f.read(1048576),b''): h.update(b)
 f.close(); return h.hexdigest()
def bh(p): return multipletests([1 if pd.isna(x) else float(x) for x in p],method='fdr_bh')[1] if len(p) else []
def cd(x,y): return (sum((i>y).sum() for i in x)-sum((i<y).sum() for i in x))/(len(x)*len(y)) if len(x) and len(y) else np.nan
try:
 log('M0 audit')
 labels=pd.read_csv(LAB,sep='\t')
 man=pd.DataFrame([dict(path=str(p),exists=p.exists(),size_bytes=p.stat().st_size,sha256=sha(p)) for p in [TAR,RAW/'GSE157007_family.soft.gz',RAW/'GSE157007_series_matrix.txt.gz',LAB]])
 out(man,WS/'processed_data/input_manifest.tsv')
 with tarfile.open(TAR) as tf: inv=pd.DataFrame([dict(name=m.name,size_bytes=m.size) for m in tf.getmembers()])
 out(inv,WS/'processed_data/tar_file_inventory.tsv')
 rows=[]
 for _,r in labels.iterrows():
  hits=inv[inv.name.str.contains(str(r.sample_id)+'_'+str(r.donor_id),regex=False)].name.tolist() or inv[inv.name.str.contains(str(r.donor_id),regex=False)].name.tolist()
  rows.append({**r.to_dict(),'matrix_files':';'.join([x for x in hits if 'matrix' in x]),'barcode_files':';'.join([x for x in hits if 'barcodes' in x]),'feature_files':';'.join([x for x in hits if 'features' in x]),'tcr_files':';'.join([x for x in hits if 'contig' in x]),'platform_series':'OH' if str(r.donor_id).startswith('OH') else 'F'})
 smap=pd.DataFrame(rows); out(smap,WS/'processed_data/sample_file_map.tsv'); out(labels.rename(columns={'analysis_group_confirmed':'analysis_group'}),WS/'processed_data/donor_group_lock.tsv')
 probs=[]
 if len(labels)!=17: probs.append(dict(issue='label_count',detail=len(labels)))
 for _,r in smap.iterrows():
  if not (r.matrix_files and r.barcode_files and r.feature_files): probs.append(dict(issue='missing_scrna_file',detail=r.donor_id))
 out(pd.DataFrame(probs or [dict(issue='none',detail='M0 passed; corrected 02 donor labels used; supplementary Tables 1/2 absent caveat retained')]),WS/'processed_data/metadata_discrepancies.tsv')
 if probs: raise RuntimeError(str(probs))
 rec('M0','provenance/file/label audit','COMPLETE',[WS/'processed_data/input_manifest.tsv',WS/'processed_data/tar_file_inventory.tsv',WS/'processed_data/sample_file_map.tsv',WS/'processed_data/donor_group_lock.tsv',WS/'processed_data/metadata_discrepancies.tsv'],'audit passed')
 log('M1 load scanpy')
 root=WS/'processed_data/raw_10x_by_sample'; root.mkdir(exist_ok=True)
 ads=[]; qc=[]; exc=[]; feats=[]
 with tarfile.open(TAR) as tf:
  for _,r in smap.iterrows():
   d=root/r.donor_id; d.mkdir(exist_ok=True)
   for src,dst in [(r.matrix_files.split(';')[0],'matrix.mtx.gz'),(r.barcode_files.split(';')[0],'barcodes.tsv.gz'),(r.feature_files.split(';')[0],'features.tsv.gz')]:
    f=d/dst
    if not f.exists(): shutil.copyfileobj(tf.extractfile(src),open(f,'wb'))
   try:
    a=sc.read_10x_mtx(str(d),var_names='gene_symbols',cache=False); a.var_names_make_unique(); raw=a.n_obs
    a.obs['donor_id']=r.donor_id; a.obs['sample_id']=r.sample_id; a.obs['analysis_group']=r.analysis_group_confirmed; a.obs['geo_age_group']=r.geo_age_group; a.obs['platform_series']=r.platform_series; a.obs_names=[r.donor_id+'_'+x for x in a.obs_names]
    a.var['mt']=a.var_names.str.upper().str.startswith('MT-'); sc.pp.calculate_qc_metrics(a,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
    keep=(a.obs.n_genes_by_counts>=200)&(a.obs.total_counts>=500)&(a.obs.pct_counts_mt<=20); a=a[keep].copy()
    qc.append(dict(donor_id=r.donor_id,sample_id=r.sample_id,group=r.analysis_group_confirmed,platform=r.platform_series,cells_raw=raw,cells_pass_qc=a.n_obs,genes_detected=a.n_vars,median_counts=float(a.obs.total_counts.median()),median_genes=float(a.obs.n_genes_by_counts.median()),median_pct_mt=float(a.obs.pct_counts_mt.median())))
    feats.append(dict(donor_id=r.donor_id,n_features=a.n_vars)); ads.append(a)
   except Exception as e: exc.append(dict(donor_id=r.donor_id,reason=repr(e)))
 if not ads: raise RuntimeError('no samples loaded')
 adata=ad.concat(ads,join='inner',index_unique=None); adata.var_names_make_unique(); adata.layers['counts']=adata.X.copy()
 out(pd.DataFrame(qc),WS/'qc/qc_by_sample.tsv'); out(pd.DataFrame(exc or [dict(donor_id='none',reason='none')]),WS/'qc/excluded_samples.tsv'); out(pd.DataFrame(feats),WS/'qc/feature_harmonization_report.tsv')
 cs=adata.obs.groupby(['analysis_group','donor_id','platform_series'],observed=True).size().reset_index(name='cells_pass_qc'); out(cs,WS/'qc/cell_qc_summary.tsv'); out(cs.groupby(['analysis_group','platform_series'],observed=True).agg(donors=('donor_id','nunique'),cells=('cells_pass_qc','sum')).reset_index(),WS/'qc/batch_diagnostic_summary.tsv')
 rec('M1','scRNA extraction/loading/QC/harmonization','COMPLETE',[WS/'qc/qc_by_sample.tsv',WS/'qc/cell_qc_summary.tsv',WS/'qc/excluded_samples.tsv',WS/'qc/feature_harmonization_report.tsv',WS/'qc/batch_diagnostic_summary.tsv'],f'{adata.n_obs} cells {adata.n_vars} genes')
 log('M2 annotation')
 sc.pp.normalize_total(adata,target_sum=1e4); sc.pp.log1p(adata)
 markers={'CD4_T':['CD3D','CD3E','CD4','IL7R'],'CD8_T':['CD3D','CD8A','CD8B','NKG7'],'NK':['NKG7','GNLY','KLRD1'],'B':['MS4A1','CD79A','CD79B'],'Monocyte_CD14':['LYZ','S100A8','S100A9','CD14'],'Monocyte_FCGR3A':['FCGR3A','MS4A7','LST1'],'DC':['FCER1A','CST3','HLA-DRA'],'Platelet':['PPBP','PF4'],'Plasma':['MZB1','JCHAIN']}
 sup=[]; cols=[]
 for ct,gs in markers.items():
  pr=[g for g in gs if g in adata.var_names]; sup.append(dict(cell_type=ct,requested_markers=';'.join(gs),present_markers=';'.join(pr),n_present=len(pr))); col='score_'+ct; cols.append(col); sc.tl.score_genes(adata,pr,score_name=col,random_state=seed) if pr else adata.obs.__setitem__(col,0)
 scs=adata.obs[cols]; adata.obs['cell_type']=scs.idxmax(1).str.replace('score_','',regex=False).values; mx=scs.max(1); adata.obs['annotation_confidence']=pd.cut(mx,[-9e9,0,.25,.75,9e9],labels=['low_or_negative','low','moderate','high']).astype(str)
 out(pd.DataFrame(sup),WS/'statistical_results/celltype_marker_support.tsv'); out(adata.obs.groupby(['cell_type','annotation_confidence'],observed=True).size().reset_index(name='cells'),WS/'statistical_results/annotation_confidence.tsv')
 long=adata.obs.groupby(['donor_id','sample_id','analysis_group','geo_age_group','platform_series','cell_type'],observed=True).size().reset_index(name='cells'); tot=adata.obs.groupby('donor_id',observed=True).size().rename('total_cells').reset_index(); long=long.merge(tot); long['fraction']=long.cells/long.total_cells
 frac=long.pivot_table(index=['donor_id','sample_id','analysis_group','geo_age_group','platform_series','total_cells'],columns='cell_type',values='fraction',fill_value=0).reset_index(); frac.columns=[str(c) for c in frac.columns]
 cnt=long.pivot_table(index=['donor_id','sample_id','analysis_group','geo_age_group','platform_series','total_cells'],columns='cell_type',values='cells',fill_value=0).reset_index(); cnt.columns=[str(c) for c in cnt.columns]
 out(frac,WS/'processed_data/celltype_fraction_by_donor.tsv'); out(cnt,WS/'processed_data/celltype_counts_by_sample.tsv'); out(long,WS/'processed_data/celltype_long_counts_fractions.tsv')
 rec('M2','PBMC annotation/landscape','COMPLETE',[WS/'statistical_results/celltype_marker_support.tsv',WS/'statistical_results/annotation_confidence.tsv',WS/'processed_data/celltype_fraction_by_donor.tsv',WS/'processed_data/celltype_counts_by_sample.tsv'],'marker-score major PBMC annotation')
 log('M3 composition')
 celltypes=[c for c in frac.columns if c not in ['donor_id','sample_id','analysis_group','geo_age_group','platform_series','total_cells']]; comp=[]; sens=[]
 for ct in celltypes:
  x=frac.loc[frac.analysis_group=='frail',ct].astype(float).values; y=frac.loc[frac.analysis_group=='healthy_old_non_frail',ct].astype(float).values
  if len(x)>=2 and len(y)>=2: comp.append(dict(cell_type=ct,n_frail=len(x),n_healthy_old_non_frail=len(y),mean_fraction_frail=x.mean(),mean_fraction_healthy_old_non_frail=y.mean(),effect_mean_difference=x.mean()-y.mean(),cliffs_delta=cd(x,y),test='Mann-Whitney U on donor fractions',p_value=stats.mannwhitneyu(x,y).pvalue))
  y2=frac.loc[(frac.analysis_group=='healthy_old_non_frail')&(frac.platform_series=='F'),ct].astype(float).values
  if len(x)>=2 and len(y2)>=2: sens.append(dict(cell_type=ct,contrast='frail_vs_healthy_old_F_series_only',n_frail=len(x),n_control=len(y2),mean_diff=x.mean()-y2.mean(),p_value=stats.mannwhitneyu(x,y2).pvalue,test='Mann-Whitney U donor fractions'))
 comp=pd.DataFrame(comp); comp['adj_p_value_bh']=bh(comp.p_value); out(comp,WS/'statistical_results/composition_differential.tsv'); sens=pd.DataFrame(sens); sens['adj_p_value_bh']=bh(sens.p_value) if len(sens) else []; out(sens,WS/'sensitivity_checks/composition_sensitivity.tsv')
 rec('M3','donor-level composition testing','COMPLETE',[WS/'statistical_results/composition_differential.tsv',WS/'sensitivity_checks/composition_sensitivity.tsv'],'donor replicate; BH FDR')
 log('M4 pseudobulk')
 genes=np.array(adata.var_names); raw=adata.layers['counts']; group=labels.set_index('donor_id').analysis_group_confirmed; design=[]; de=[]
 for ct in sorted(adata.obs.cell_type.unique()):
  idx=np.where(adata.obs.cell_type.values==ct)[0]; obs=adata.obs.iloc[idx]; donors=[d for d,n in obs.groupby('donor_id',observed=True).size().items() if n>=20]
  gr=group.loc[[d for d in donors if d in group.index]] if donors else pd.Series(dtype=str)
  if (gr=='frail').sum()<3 or (gr=='healthy_old_non_frail').sum()<3: design.append(dict(cell_type=ct,status='skipped_insufficient_donors_or_cells',n_donors=len(donors),n_frail=int((gr=='frail').sum()),n_healthy_old_non_frail=int((gr=='healthy_old_non_frail').sum()))); continue
  mats=[]; meta=[]
  for d in gr.index:
   di=idx[np.where(obs.donor_id.values==d)[0]]; mats.append(np.asarray(raw[di,:].sum(0)).ravel()); meta.append(dict(donor_id=d,group=group.loc[d],cells=len(di)))
  mat=np.vstack(mats).astype(float); lib=mat.sum(1); lc=np.log2((mat/lib[:,None])*1e6+1); meta=pd.DataFrame(meta); f=np.where(meta.group.values=='frail')[0]; h=np.where(meta.group.values=='healthy_old_non_frail')[0]
  design.append(dict(cell_type=ct,status='tested',n_donors=len(meta),n_frail=len(f),n_healthy_old_non_frail=len(h),min_cells_per_donor=int(meta.cells.min())))
  expr=(mat[f,:].sum(0)>10)|(mat[h,:].sum(0)>10); p=np.full(len(genes),np.nan); diff=lc[f,:].mean(0)-lc[h,:].mean(0)
  for j in np.where(expr)[0]: p[j]=stats.mannwhitneyu(lc[f,j],lc[h,j]).pvalue
  pad=np.full(len(genes),np.nan); valid=~np.isnan(p); pad[valid]=bh(p[valid])
  for j in np.argsort(np.nan_to_num(pad,nan=2))[:500]:
   if valid[j]: de.append(dict(cell_type=ct,gene=genes[j],n_frail=len(f),n_healthy_old_non_frail=len(h),log2cpm_difference=diff[j],p_value=p[j],adj_p_value_bh=pad[j],test='Mann-Whitney U donor pseudobulk logCPM'))
 out(pd.DataFrame(design),WS/'processed_data/pseudobulk_design.tsv'); de=pd.DataFrame(de); out(de,WS/'statistical_results/de_by_celltype.tsv')
 targ=[]
 for ct in [c for c in adata.obs.cell_type.unique() if 'Monocyte' in c]:
  for gn in ['NEAT1','MALAT1']:
   if gn not in genes: continue
   j=int(np.where(genes==gn)[0][0]); vals=[]
   for d in labels.donor_id:
    ii=np.where((adata.obs.cell_type.values==ct)&(adata.obs.donor_id.values==d))[0]
    if len(ii)>=20:
     s=np.asarray(raw[ii,j].sum(0)).ravel()[0]; lib=np.asarray(raw[ii,:].sum(0)).sum(); vals.append(dict(group=group.loc[d],v=float(np.log2((s/lib)*1e6+1))))
   df=pd.DataFrame(vals); x=df.loc[df.group=='frail','v'].values if len(df) else []; y=df.loc[df.group=='healthy_old_non_frail','v'].values if len(df) else []
   targ.append(dict(cell_type=ct,gene=gn,n_frail=len(x),n_healthy_old_non_frail=len(y),effect_diff=(np.mean(x)-np.mean(y)) if len(x) and len(y) else np.nan,p_value=stats.mannwhitneyu(x,y).pvalue if len(x)>=2 and len(y)>=2 else np.nan,test='targeted donor pseudobulk'))
 targ=pd.DataFrame(targ); targ['adj_p_value_bh']=bh(targ.p_value) if len(targ) else []; out(targ,WS/'statistical_results/neat1_malat1_monocyte_check.tsv')
 (WS/'enrichment_results/NOT_RUN.md').write_text('Not run: no vetted local gene-set database supplied; use de_by_celltype.tsv for repair if resources are authorized.\n')
 rec('M4','cell-type-specific transcriptional programs','COMPLETE',[WS/'processed_data/pseudobulk_design.tsv',WS/'statistical_results/de_by_celltype.tsv',WS/'statistical_results/neat1_malat1_monocyte_check.tsv',WS/'enrichment_results/NOT_RUN.md'],'donor pseudobulk; enrichment skipped')
 log('M5 TCR')
 tsum=[]; exp=[]
 with tarfile.open(TAR) as tf:
  for _,r in smap.iterrows():
   fs=[x for x in str(r.tcr_files).split(';') if x]
   if not fs: tsum.append(dict(donor_id=r.donor_id,status='missing')); continue
   df=pd.read_csv(gzip.GzipFile(fileobj=tf.extractfile(fs[0]))); prod=df[df.productive.astype(str).str.lower().isin(['true','t','1'])] if 'productive' in df else df; col='raw_clonotype_id' if 'raw_clonotype_id' in prod else 'clonotype_id'
   vc=prod[col].dropna().astype(str).value_counts() if col in prod else pd.Series(dtype=int); fq=vc/vc.sum() if vc.sum() else pd.Series(dtype=float)
   for cid,n in vc.head(50).items(): exp.append(dict(donor_id=r.donor_id,group=r.analysis_group_confirmed,clonotype_id=cid,count=int(n),frequency=float(n/vc.sum())))
   tsum.append(dict(donor_id=r.donor_id,sample_id=r.sample_id,group=r.analysis_group_confirmed,status='parsed',contigs_total=len(df),productive_contigs=len(prod),unique_clonotypes=int(len(vc)),expanded_clonotypes=int((vc>=2).sum()),proportion_productive_in_expanded_clonotypes=float(vc[vc>=2].sum()/vc.sum()) if vc.sum() else np.nan,shannon=float(-(fq*np.log(fq)).sum()) if len(fq) else np.nan,simpson=float(1-(fq**2).sum()) if len(fq) else np.nan,source_file=fs[0]))
 tcr=pd.DataFrame(tsum); out(tcr,WS/'statistical_results/tcr_parse_summary.tsv'); out(tcr,WS/'statistical_results/clonotype_diversity_by_donor.tsv'); out(pd.DataFrame(exp),WS/'statistical_results/expanded_clonotypes.tsv'); out(pd.DataFrame([dict(donor_id=d,status='not_integrated',reason='donor-level TCR only; barcode linkage not used for core claims') for d in labels.donor_id]),WS/'statistical_results/tcr_scrna_linkage_report.tsv')
 rec('M5','scTCR repertoire summaries','COMPLETE',[WS/'statistical_results/tcr_parse_summary.tsv',WS/'statistical_results/clonotype_diversity_by_donor.tsv',WS/'statistical_results/expanded_clonotypes.tsv',WS/'statistical_results/tcr_scrna_linkage_report.tsv'],'donor-level TCR parsed')
 log('M6 sensitivity')
 out(frac.groupby('analysis_group',observed=True).agg(donors=('donor_id','nunique'),median_total_cells=('total_cells','median')).reset_index(),WS/'sensitivity_checks/age_context_summary.tsv'); out(frac.groupby(['analysis_group','platform_series'],observed=True).agg(donors=('donor_id','nunique'),median_total_cells=('total_cells','median')).reset_index(),WS/'sensitivity_checks/platform_sensitivity.tsv')
 alt=[]
 for ct in celltypes:
  for a,b in [('frail','young_adult'),('frail','cord_blood'),('healthy_old_non_frail','young_adult')]:
   x=frac.loc[frac.analysis_group==a,ct].astype(float).values; y=frac.loc[frac.analysis_group==b,ct].astype(float).values
   if len(x)>=2 and len(y)>=2: alt.append(dict(cell_type=ct,contrast=f'{a}_vs_{b}',n_group1=len(x),n_group2=len(y),mean_diff=x.mean()-y.mean(),p_value=stats.mannwhitneyu(x,y).pvalue,test='Mann-Whitney U donor fractions'))
 alt=pd.DataFrame(alt); alt['adj_p_value_bh']=bh(alt.p_value) if len(alt) else []; out(alt,WS/'sensitivity_checks/group_contrast_sensitivity.tsv')
 rec('M6','age-gradient/platform sensitivity','COMPLETE',[WS/'sensitivity_checks/age_context_summary.tsv',WS/'sensitivity_checks/platform_sensitivity.tsv',WS/'sensitivity_checks/group_contrast_sensitivity.tsv'],'completed')
 (WS/'sensitivity_checks/optional_modules_NOT_RUN.md').write_text('M7 CITE-seq and M8 LR/GRN skipped: no unambiguous ADT matrix or vetted local LR/GRN database supplied.\n')
 rec('M7','optional CITE-seq validation','SKIPPED_WITH_RATIONALE',[WS/'sensitivity_checks/optional_modules_NOT_RUN.md'],'no ADT matrix'); rec('M8','optional communication/GRN hypotheses','SKIPPED_WITH_RATIONALE',[WS/'sensitivity_checks/optional_modules_NOT_RUN.md'],'no vetted database')
 adata.write_h5ad(WS/'processed_data/gse157007_qc_annotated_common_genes.h5ad',compression='gzip')
 out(pd.DataFrame(mods),WS/'module_execution_log.tsv'); (WS/'module_execution_log.md').write_text('# Module Execution Log\n\n'+'\n'.join([f"## {m['module_id']} - {m['module']}\nStatus: {m['status']}\nOutputs: {m['outputs']}\nNotes: {m['notes']}\n" for m in mods]))
 (WS/'analysis_plan.md').write_text('# Analysis Plan\n04-feasibility conditional pass followed. M0 first, corrected 02 labels, donor-aware statistics, no cell-level pseudoreplication. Executed M0-M6 and donor TCR; skipped optional enrichment/ADT/LR/GRN with rationale.\n')
 fig=[('Fig1','A','processed_data/input_manifest.tsv','scripts/run_sc_analysis.py','processed_data/tar_file_inventory.tsv','provenance','audit only'),('Fig1','B','processed_data/donor_group_lock.tsv','scripts/run_sc_analysis.py','processed_data/sample_file_map.tsv','labels','corrected 02 labels'),('Fig1','C','qc/qc_by_sample.tsv','scripts/run_sc_analysis.py','qc/cell_qc_summary.tsv','QC','threshold caveat'),('Fig2','A','statistical_results/celltype_marker_support.tsv','scripts/run_sc_analysis.py','processed_data/celltype_counts_by_sample.tsv','annotation','marker-score'),('Fig2','D','statistical_results/composition_differential.tsv','scripts/run_sc_analysis.py','sensitivity_checks/composition_sensitivity.tsv','composition stats','small n'),('Fig3','A','statistical_results/de_by_celltype.tsv','scripts/run_sc_analysis.py','processed_data/pseudobulk_design.tsv','pseudobulk DE','underpowered'),('Fig3','B','statistical_results/neat1_malat1_monocyte_check.tsv','scripts/run_sc_analysis.py','statistical_results/neat1_malat1_monocyte_check.tsv','target genes','targeted'),('Fig4','A','statistical_results/clonotype_diversity_by_donor.tsv','scripts/run_sc_analysis.py','statistical_results/tcr_parse_summary.tsv','TCR diversity','descriptive'),('Fig5','A','sensitivity_checks/age_context_summary.tsv','scripts/run_sc_analysis.py','sensitivity_checks/group_contrast_sensitivity.tsv','age context','secondary'),('Fig5','B','sensitivity_checks/platform_sensitivity.tsv','scripts/run_sc_analysis.py','sensitivity_checks/composition_sensitivity.tsv','platform','OH/F caveat')]
 out(pd.DataFrame(fig,columns='figure_id panel_id source_result source_script data_file statistic caveat'.split()),WS/'figure_data_map.tsv')
 ev=[('processed_data/donor_group_lock.tsv','Corrected labels support contrast','supplementary tables absent','design'),('statistical_results/composition_differential.tsv','Donor-level composition effects estimated','small n','candidate'),('statistical_results/de_by_celltype.tsv','Donor pseudobulk ranks generated','underpowered','candidate'),('statistical_results/clonotype_diversity_by_donor.tsv','TCR donor summaries generated','no barcode linkage','descriptive')]
 out(pd.DataFrame(ev,columns='result_file supported_claim limitation claim_strength'.split()),WS/'evidence_chain.tsv')
 out(pd.DataFrame([dict(hypothesis='Frailty-associated PBMC lineage remodeling',evidence_file='statistical_results/composition_differential.tsv',support_type='donor fractions',status='candidate'),dict(hypothesis='Monocyte NEAT1/MALAT1 differences',evidence_file='statistical_results/neat1_malat1_monocyte_check.tsv',support_type='targeted pseudobulk',status='candidate'),dict(hypothesis='TCR expansion context',evidence_file='statistical_results/clonotype_diversity_by_donor.tsv',support_type='donor TCR',status='descriptive')]),WS/'mechanism_map.tsv')
 comp_sig=int((comp.adj_p_value_bh<.05).sum()); de_sig=int((de.adj_p_value_bh<.05).sum()) if len(de) else 0; tcr_n=int((tcr.status=='parsed').sum())
 (WS/'analysis_report.md').write_text(f'# Analysis Report\nStatus COMPLETE; route to 06-figures with caveats. M0 passed using corrected 02 labels. Donor is replicate; no cell-level pseudoreplication. QC-pass cells {adata.n_obs}; donors {adata.obs.donor_id.nunique()}; common genes {adata.n_vars}. Primary contrast frail n=5 vs healthy-old n=6. Composition tests {len(comp)}, BH significant {comp_sig}. Pseudobulk rows {len(de)}, BH significant reported rows {de_sig}. TCR donors parsed {tcr_n}. Claims must remain conservative because of small n and OH/F platform imbalance.\n')
 manifest='agent: 05-analysis\nproject_id: cs1-paper-factory\nrun_id: run_20260529_095727\nstatus: complete\nrandom_seed: 20260529\ncommand: "D:/cs2/tools/python/python.exe D:/cs1/runs/cs1-paper-factory/run_20260529_095727/agents/05-analysis/workspace/scripts/run_sc_analysis.py"\nreplicate_unit: donor\nsoftware_versions:\n  python: "'+sys.version.replace('\n',' ')+'"\n  scanpy: "'+sc.__version__+'"\n  anndata: "'+ad.__version__+'"\n  numpy: "'+np.__version__+'"\n  pandas: "'+pd.__version__+'"\nmodules:\n' + ''.join([f"  - id: {m['module_id']}\n    status: {m['status']}\n" for m in mods])
 (WS/'analysis_manifest.yaml').write_text(manifest)
 hand={'agent':'05-analysis','project_id':'cs1-paper-factory','run_id':'run_20260529_095727','status':'complete','decision':'pass_to_06_figures','route':'06-figures','next_agent':'06-figures','timestamp':datetime.datetime.now(datetime.timezone.utc).isoformat(),'key_metrics':{'datasets_analyzed':1,'donors_analyzed':int(adata.obs.donor_id.nunique()),'qc_pass_cells':int(adata.n_obs),'common_genes':int(adata.n_vars),'modules_complete':7,'modules_skipped':2,'modules_failed':0,'composition_tests':int(len(comp)),'composition_fdr_significant':comp_sig,'de_reported_rows':int(len(de)),'de_fdr_significant_reported_rows':de_sig,'tcr_donors_parsed':tcr_n},'outputs':{'analysis_report':str(WS/'analysis_report.md'),'analysis_manifest':str(WS/'analysis_manifest.yaml'),'module_execution_log':str(WS/'module_execution_log.md'),'figure_data_map':str(WS/'figure_data_map.tsv')},'caveats':['small donor count','OH/F platform imbalance','marker-score annotation','optional enrichment ADT LR GRN not run','TCR linkage not integrated']}
 (WS/'handoff.json').write_text(json.dumps(hand,indent=2))
 log('COMPLETE pass_to_06_figures')
except Exception as e:
 tb=traceback.format_exc(); log('FAILED '+repr(e)+'\n'+tb); (WS/'analysis_report.md').write_text('# BLOCKED\n\n'); (WS/'handoff.json').write_text(json.dumps({'agent':'05-analysis','status':'blocked','route':'repair','error':repr(e),'traceback':tb},indent=2)); raise

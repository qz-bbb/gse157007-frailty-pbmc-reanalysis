
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
BASE = Path('D:/cs1/runs/cs1-paper-factory/run_20260529_095727') if os.name == 'nt' else Path('/mnt/d/cs1/runs/cs1-paper-factory/run_20260529_095727')
AN = BASE/'05-analysis/workspace'
OUT = BASE/'agents/06-figures/workspace/figures'
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'figure.dpi':150,'savefig.dpi':300})
COLORS={'frail':'#D55E00','healthy_old_non_frail':'#0072B2','young_adult':'#009E73','cord_blood':'#CC79A7','F':'#56B4E9','OH':'#E69F00'}
def read(rel): return pd.read_csv(AN/rel, sep='\t')
def label(ax, s): ax.text(-0.12,1.08,s,transform=ax.transAxes,fontsize=14,fontweight='bold',va='top')
def save(fig, name):
    p=OUT/name; fig.savefig(p, bbox_inches='tight', facecolor='white'); plt.close(fig); return str(p)
def fig1():
    im=read('processed_data/input_manifest.tsv'); sample=read('processed_data/sample_file_map.tsv'); qc=read('qc/qc_by_sample.tsv')
    fig,axs=plt.subplots(2,2,figsize=(10,7.5),constrained_layout=True)
    ax=axs[0,0]; label(ax,'A'); im2=im.copy(); im2['source']=im2['path'].apply(lambda x: Path(str(x).replace(chr(92),'/')).name[:34]); ax.barh(im2['source'], im2['size_bytes']/1e6, color='#0072B2'); ax.set_xlabel('File size (MB)'); ax.set_title('Input provenance files verified')
    for i,v in enumerate(im2['exists']): ax.text(0.02,i,'exists' if v else 'missing',va='center',fontsize=7)
    ax=axs[0,1]; label(ax,'B'); ct=sample.groupby(['analysis_group_confirmed','platform_series']).size().unstack(fill_value=0); ct.plot(kind='bar', stacked=True, ax=ax, color=[COLORS.get(c,'gray') for c in ct.columns]); ax.set_ylabel('Donors'); ax.set_xlabel('Analysis group'); ax.set_title('Corrected donor group lock'); ax.tick_params(axis='x',rotation=30)
    ax=axs[1,0]; label(ax,'C'); qc2=qc.sort_values('cells_pass_qc', ascending=False); ax.bar(qc2['donor_id'], qc2['cells_pass_qc'], color=[COLORS.get(g,'gray') for g in qc2['group']]); ax.set_ylabel('QC-pass cells'); ax.set_xlabel('Donor'); ax.set_title('QC-pass cells by donor'); ax.tick_params(axis='x',rotation=90)
    ax=axs[1,1]; label(ax,'D')
    for g,d in qc.groupby('group'):
        ax.scatter(d['median_genes'], d['median_pct_mt'], s=d['cells_pass_qc']/100, label=g, color=COLORS.get(g,'gray'), alpha=.75, edgecolor='black', linewidth=.4)
    ax.set_xlabel('Median genes per cell'); ax.set_ylabel('Median mitochondrial %'); ax.set_title('QC diagnostic by donor'); ax.legend(frameon=False, loc='best')
    return save(fig,'Fig1_dataset_qc.png')
def fig2():
    m=read('statistical_results/celltype_marker_support.tsv').sort_values('n_present'); comp=read('statistical_results/composition_differential.tsv').sort_values('effect_mean_difference')
    fig,axs=plt.subplots(1,2,figsize=(10,4.6),constrained_layout=True)
    ax=axs[0]; label(ax,'A'); ax.barh(m['cell_type'],m['n_present'],color='#009E73'); ax.set_xlabel('Requested markers present'); ax.set_title('Cell-type marker support')
    for y,row in enumerate(m.itertuples()): ax.text(row.n_present+.05,y,str(row.present_markers),va='center',fontsize=6)
    ax=axs[1]; label(ax,'B'); colors=['#D55E00' if v>0 else '#0072B2' for v in comp['effect_mean_difference']]; ax.barh(comp['cell_type'],comp['effect_mean_difference']*100,color=colors); ax.axvline(0,color='black',lw=.8); ax.set_xlabel('Mean fraction difference, frail - healthy-old (percentage points)'); ax.set_title('Donor-level composition contrast')
    for y,row in enumerate(comp.itertuples()): ax.text(row.effect_mean_difference*100,y,f' q={row.adj_p_value_bh:.2g}',va='center',fontsize=7)
    return save(fig,'Fig2_annotation_composition.png')
def fig3():
    de=read('statistical_results/de_by_celltype.tsv'); mon=de[de['cell_type'].eq('Monocyte_CD14')].copy(); neat=read('statistical_results/neat1_malat1_monocyte_check.tsv'); sens=read('sensitivity_checks/mjr001_f_series_summary.tsv').sort_values('min_adj_p_value_retained')
    fig,axs=plt.subplots(1,3,figsize=(13,4.3),constrained_layout=True)
    ax=axs[0]; label(ax,'A'); ax.scatter(mon['log2cpm_difference'], -np.log10(mon['p_value'].clip(lower=1e-300)), c=mon['adj_p_value_bh']<0.05, cmap='cividis', s=16, alpha=.8); ax.axvline(0,color='gray',lw=.8); ax.set_xlabel('log2CPM difference'); ax.set_ylabel('-log10 p'); ax.set_title('Monocyte_CD14 mixed-series pseudobulk\n(exploratory)')
    for _,r in mon.nsmallest(5,'adj_p_value_bh').iterrows(): ax.text(r['log2cpm_difference'], -np.log10(max(r['p_value'],1e-300)), r['gene'], fontsize=6)
    ax=axs[1]; label(ax,'B'); neat['label']=neat['cell_type']+' '+neat['gene']; colors=['#D55E00' if v>0 else '#0072B2' for v in neat['effect_diff']]; ax.barh(neat['label'], neat['effect_diff'], color=colors); ax.axvline(0,color='black',lw=.8); ax.set_xlabel('Effect difference'); ax.set_title('Targeted NEAT1/MALAT1 check')
    for y,r in enumerate(neat.itertuples()): ax.text(r.effect_diff,y,f' q={r.adj_p_value_bh:.2g}',va='center',fontsize=7)
    ax=axs[2]; label(ax,'C'); ax.barh(sens['cell_type'], sens['min_adj_p_value_retained'], color='#CC79A7'); ax.axvline(.05,color='black',ls='--',lw=.8); ax.set_xlabel('Minimum retained BH q-value'); ax.set_title('F-series-only sensitivity: no retained q<0.05')
    return save(fig,'Fig3_transcription_sensitivity.png')
def fig4():
    t=read('statistical_results/clonotype_diversity_by_donor.tsv'); fig,axs=plt.subplots(1,3,figsize=(12,4.2),constrained_layout=True)
    ax=axs[0]; label(ax,'A'); groups=list(t['group'].drop_duplicates()); data=[t.loc[t.group==g,'shannon'] for g in groups]; bp=ax.boxplot(data, labels=groups, patch_artist=True)
    for patch,g in zip(bp['boxes'],groups): patch.set_facecolor(COLORS.get(g,'lightgray'))
    ax.set_ylabel('Shannon diversity'); ax.set_title('Donor TCR diversity'); ax.tick_params(axis='x',rotation=30)
    ax=axs[1]; label(ax,'B'); ax.scatter(t['unique_clonotypes'], t['expanded_clonotypes'], c=[COLORS.get(g,'gray') for g in t['group']], s=50, edgecolor='black', linewidth=.4); ax.set_xlabel('Unique clonotypes'); ax.set_ylabel('Expanded clonotypes'); ax.set_title('Clonotype expansion by donor')
    ax=axs[2]; label(ax,'C'); parsed=t.groupby('status').size(); ax.bar(parsed.index, parsed.values, color='#56B4E9'); ax.set_ylabel('Donors'); ax.set_title('TCR parse status')
    return save(fig,'Fig4_tcr_diversity.png')
def fig5():
    age=read('sensitivity_checks/age_context_summary.tsv'); plat=read('sensitivity_checks/platform_sensitivity.tsv'); cs=read('sensitivity_checks/composition_sensitivity.tsv').sort_values('mean_diff')
    fig,axs=plt.subplots(1,3,figsize=(13,4.3),constrained_layout=True)
    ax=axs[0]; label(ax,'A'); ax.bar(age['analysis_group'],age['donors'],color=[COLORS.get(g,'gray') for g in age['analysis_group']]); ax.set_ylabel('Donors'); ax.set_title('Age/context donor counts'); ax.tick_params(axis='x',rotation=30); ax2=ax.twinx(); ax2.plot(age['analysis_group'],age['median_total_cells'],color='black',marker='o'); ax2.set_ylabel('Median total cells')
    ax=axs[1]; label(ax,'B'); pt=plat.pivot(index='analysis_group',columns='platform_series',values='donors').fillna(0); pt.plot(kind='bar',stacked=True,ax=ax,color=[COLORS.get(c,'gray') for c in pt.columns]); ax.set_ylabel('Donors'); ax.set_title('Platform/series balance'); ax.tick_params(axis='x',rotation=30)
    ax=axs[2]; label(ax,'C'); ax.barh(cs['cell_type'],cs['mean_diff']*100,color=['#D55E00' if x>0 else '#0072B2' for x in cs['mean_diff']]); ax.axvline(0,color='black',lw=.8); ax.set_xlabel('F-series mean fraction difference (percentage points)'); ax.set_title('Composition sensitivity (F series only)')
    for y,r in enumerate(cs.itertuples()): ax.text(r.mean_diff*100,y,f' q={r.adj_p_value_bh:.2g}',va='center',fontsize=7)
    return save(fig,'Fig5_context_sensitivity.png')
if __name__=='__main__':
    for f in [fig1,fig2,fig3,fig4,fig5]: print(f())

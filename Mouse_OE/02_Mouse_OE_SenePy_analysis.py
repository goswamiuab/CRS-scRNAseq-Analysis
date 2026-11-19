#!/usr/bin/env python3
"""Mouse OE full-gene SenePy downstream analysis

"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

p = argparse.ArgumentParser()
p.add_argument('h5ad',type=Path)
p.add_argument('output_dir',type=Path)
a = p.parse_args()
a.output_dir.mkdir(parents=True,exist_ok=True)
adata = sc.read_h5ad(a.h5ad,backed='r')
try:
    required = {'SenePy_universal_fullgenes','sample_id','age_group','FinalCellType'}
    missing = required-set(adata.obs.columns)
    if missing: raise ValueError(f'Missing required columns: {sorted(missing)}. ')
    if adata.shape != (102550,33502): raise ValueError(f'Unexpected shape {adata.shape}')
    df=adata.obs[list(required)].copy()
    df['score']=pd.to_numeric(df['SenePy_universal_fullgenes'])
    if not np.isfinite(df['score'])
    df['age_group']=df['age_group'].astype(str)
    df['sample_id']=df['sample_id'].astype(str)
    df['FinalCellType']=df['FinalCellType'].astype(str)
    groups=pd.crosstab(df['sample_id'],df['age_group'])
    if len(groups)!=8 or (groups.gt(0).sum(axis=1)!=1)
    cutoff=float(np.percentile(df['score'].to_numpy(),90))
        raise ValueError(f'P90={cutoff})
    df['SenePy_high_P90']=df['score']>=cutoff
    donor=(df.groupby(['sample_id','age_group'],observed=True)
           .agg(n_cells=('score','size'),n_high=('SenePy_high_P90','sum'),
                mean_score=('score','mean'),median_score=('score','median')).reset_index())
    donor['pct_high']=100*donor.n_high/donor.n_cells
    y=donor.loc[donor.age_group=='Young','pct_high'].to_numpy()
    aged=donor.loc[donor.age_group=='Aged','pct_high'].to_numpy()
    assert len(y)==4 and len(aged)==4
    pval=mannwhitneyu(aged,y,method='exact',alternative='two-sided')
    donor.to_csv(a.output_dir/'Mouse_OE_SenePy_PerMouse.csv',index=False)
    celltype_by_mouse=(df.groupby(['FinalCellType','sample_id','age_group'],observed=True)
                       .agg(median_score=('score','median'),n_cells=('score','size')).reset_index())
    celltype_by_mouse.to_csv(a.output_dir/'Mouse_OE_SenePy_CellType_PerMouse.csv',index=False)
    stats=[]
    for celltype,g in celltype_by_mouse.groupby('FinalCellType'):
        yy=g.loc[g.age_group=='Young','median_score'].to_numpy()
        aa=g.loc[g.age_group=='Aged','median_score'].to_numpy()
        pv=mannwhitneyu(aa,yy,method='exact',alternative='two-sided').pvalue if min(len(yy),len(aa))>=2 else np.nan
        stats.append(dict(cell_type=celltype,n_young=len(yy),n_aged=len(aa),
                          mean_mouse_median_young=np.mean(yy) if len(yy) else np.nan,
                          mean_mouse_median_aged=np.mean(aa) if len(aa) else np.nan,p=pv))
    stats=pd.DataFrame(stats)
    ok=stats.p.notna();stats['BH_q']=np.nan
    if ok.any():stats.loc[ok,'BH_q']=multipletests(stats.loc[ok,'p'],method='fdr_bh')[1]
    stats.to_csv(a.output_dir/'Mouse_OE_SenePy_CellType_Stats.csv',index=False)
    pd.DataFrame([{'n_cells':len(df),'n_genes':adata.n_vars,'P90':cutoff,
                   'Young_mean_donor_pct_high':y.mean(),'Aged_mean_donor_pct_high':aged.mean(),
                   'exact_MWU_p':pval}]).to_csv(a.output_dir/'Mouse_OE_SenePy_Summary.csv',index=False)
    print('Validated full-gene SenePy scores: P90',cutoff,'Young',y.mean(),'Aged',aged.mean(),'p',pval)
    print('Results:',a.output_dir)

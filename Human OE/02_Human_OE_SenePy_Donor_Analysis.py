#!/usr/bin/env python3
"""Human OE: Analysis


"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests

DEFAULT = '/data/user/goswamis/Human_OE_Reanalysis_2026/Human_OE_64249_SenePy.h5ad'

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--input', default=DEFAULT)
    p.add_argument('--output', default='Human_OE_SenePy_Analysis')
    p.add_argument('--manifest', help='CSV with sample,dataset columns ()')
    p.add_argument('--score-column', default='senepy_score')
    args=p.parse_args()
    out=Path(args.output); out.mkdir(parents=True,exist_ok=True)
    a=sc.read_h5ad(args.input, backed='r')
    need=['sample','group','cell_type_final',args.score_column]
    absent=[x for x in need if x not in a.obs]
    if absent: raise ValueError(f'Missing required obs columns: {absent})
    df=a.obs[need].copy().reset_index().rename(columns={'index':'cell_id','sample':'donor','group':'age_group','cell_type_final':'cell_type',args.score_column:'senepy_score'})
    df['senepy_score']=pd.to_numeric(df['senepy_score'],errors='raise')
    assert np.isfinite(df['senepy_score']).all()
    assert df.age_group.astype(str).isin(['Young','Aged']).all()
    p90=float(df.senepy_score.quantile(.90))
    df['senepy_high']=df.senepy_score>=p90
    df[['cell_id','donor','age_group','cell_type','senepy_score','senepy_high']].to_csv(out/'Human_OE_SenePy_cell_scores.csv',index=False)
    donor=df.groupby(['donor','age_group'],observed=True).agg(n_cells=('senepy_score','size'), n_high=('senepy_high','sum'), mean_score=('senepy_score','mean'),median_score=('senepy_score','median')).reset_index()
    donor['pct_high']=100*donor.n_high/donor.n_cells
    if donor.donor.nunique()!=10: print(expected 10 donors; found',donor.donor.nunique())
    donor.to_csv(out/'Human_OE_SenePy_donor_summary.csv',index=False)
    young=donor.loc[donor.age_group=='Young','pct_high'].to_numpy()
    aged=donor.loc[donor.age_group=='Aged','pct_high'].to_numpy()
    if min(len(young),len(aged))>=2:
        test=mannwhitneyu(aged,young,alternative='two-sided',method='exact')
        pd.DataFrame([dict(comparison='Aged_vs_Young',n_young=len(young),n_aged=len(aged),young_mean_pct=young.mean(),aged_mean_pct=aged.mean(),p_value=test.pvalue)]).to_csv(out/'Human_OE_SenePy_donor_age_test.csv',index=False)
    # Donor-by-cell-type summaries and donor-level tests, not pseudoreplicated cell tests.
    ct=df.groupby(['donor','age_group','cell_type'],observed=True).agg(n_cells=('senepy_score','size'),median_score=('senepy_score','median'),mean_score=('senepy_score','mean'),n_high=('senepy_high','sum')).reset_index()
    ct['pct_high']=100*ct.n_high/ct.n_cells
    ct.to_csv(out/'Human_OE_SenePy_celltype_by_donor.csv',index=False)
    rows=[]
    for name,g in ct.groupby('cell_type',observed=True):
        y=g.loc[g.age_group=='Young','median_score'].to_numpy()
        o=g.loc[g.age_group=='Aged','median_score'].to_numpy()
        pv=mannwhitneyu(o,y,alternative='two-sided',method='exact').pvalue if min(len(y),len(o))>=2 else np.nan
        rows.append(dict(cell_type=name,n_young=len(y),n_aged=len(o),p_value=pv))
    tests=pd.DataFrame(rows)
    tests['q_value']=np.nan
    mask=tests.p_value.notna()
    if mask.any(): tests.loc[mask,'q_value']=multipletests(tests.loc[mask,'p_value'],method='fdr_bh')[1]
    tests.to_csv(out/'Human_OE_SenePy_celltype_donor_statistics.csv',index=False)
    
    if args.manifest:
        m=pd.read_csv(args.manifest,dtype=str)
        if not {'sample','dataset'}.issubset(m.columns): raise ValueError('Manifest requires sample,dataset columns')
        if m['sample'].duplicated().any(): raise ValueError('Duplicate sample in manifest')
        donor=donor.merge(m[['sample','dataset']],left_on='donor',right_on='sample',validate='one_to_one',how='left')
        if donor.dataset.isna().any(): raise ValueError('Some donors missing GEO source in manifest')
        donor['aged']=(donor.age_group=='Aged').astype(int)
        for response in ('mean_score','pct_high'):
            fit=ols(f'{response} ~ aged + C(dataset)',data=donor).fit(cov_type='HC3')
            pd.DataFrame({'term':fit.params.index,'beta':fit.params.values,'SE':fit.bse.values,'p_value':fit.pvalues.values,'CI_low':fit.conf_int()[0].values,'CI_high':fit.conf_int()[1].values}).to_csv(out/f'Human_OE_adjusted_{response}.csv',index=False)
        donor.to_csv(out/'Human_OE_SenePy_donor_with_source.csv',index=False)
    else:
        pd.DataFrame({'sample':sorted(donor.donor.unique()),'dataset':['VERIFY_GEO_SOURCE']*len(donor)}).to_csv(out/'GEO_source_manifest_TEMPLATE.csv',index=False)
        print('Source-adjusted regression skipped: --manifest with verified GEO mapping not supplied.')
    (out/'RUN_INFO.txt').write_text(f'input={args.input}\ncells={len(df)}\ndonors={donor.donor.nunique()}\nglobal_senepy_p90={p90:.16g}\nscore_column={args.score_column}\n',encoding='utf-8')
    a.file.close()
    print('Analysis complete. Global P90:',p90,'; output:',out)
if __name__=='__main__': main()

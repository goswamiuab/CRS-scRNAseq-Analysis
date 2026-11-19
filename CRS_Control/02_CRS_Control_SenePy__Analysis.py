#!/usr/bin/env python3
"""Human CRS/control full-gene SenePy v1.0.1 and biological-donor analysis.
 python 02_CRS_SenePy_and_Donor_Analysis.py --input CRS_Control_72262_FULLGENES_SenePy_FINAL.h5ad 
 

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--out', required=True)
parser.add_argument('--rescore',action='store_true')
a=parser.parse_args()
out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
adata=sc.read_h5ad(a.input)
assert adata.shape == (72262, 27204), f'Unexpected shape {adata.shape}'
for col in ('orig.ident','group4','cell_type_final','scDblFinder.class'):
    assert col in adata.obs, f'Missing {col}'
assert adata.obs['orig.ident'].nunique() == 16
assert adata.obs['scDblFinder.class'].astype(str).eq('singlet').all()
score_col='SenePy_universal_fullgenes'
if a.rescore:
    import senepy as sp
    from scipy import sparse
    assert sparse.issparse(adata.X), '
    hubs=sp.load_hubs(species='Human',sig_type='cell_type')
    hubs.merge_hubs(hubs.metadata,new_name='Universal',calculate_thresh=True,p_thres=0.01)
    hub=hubs.hubs['Universal']
    genes=[str(g) for g,w in hub]
    coverage=len(set(genes)&set(adata.var_names))
    print(f'Signature coverage {coverage}/{len(genes)}')
    assert len(genes)==734 and coverage==733, 'Signature coverage'
    # Official method densifies internally. Run only within sufficiently provisioned HPC allocation.
    np.random.seed(0)
    scores=np.asarray(sp.score_hub(adata,hub,n_bins=25,ctrl_size=50,binarize=True,importance=True,translator=None,verbose=True),dtype=float)
    assert len(scores)==adata.n_obs and np.isfinite(scores).all()
    adata.obs[score_col]=scores
else:
    if score_col not in adata.obs: raise ValueError('Saved final SenePy')
    scores=pd.to_numeric(adata.obs[score_col],errors='raise').to_numpy(dtype=float)
    assert np.isfinite(scores).all()
p90=float(np.percentile(scores,90))
adata.obs['SenePy_high_P90']=scores>=p90
print('Global P90:',p90,'high cells:',int(adata.obs['SenePy_high_P90'].sum()))
if not np.isclose(p90,1.505203997154884,atol=1e-8):
    ('P90 differs from verified final CRS analysis')
obs=adata.obs
per_donor=(obs.groupby(['orig.ident','group4'],observed=True).agg(n_cells=(score_col,'size'),n_high=('SenePy_high_P90','sum'),mean_score=(score_col,'mean'),median_score=(score_col,'median')).reset_index())
per_donor['pct_high']=100*per_donor['n_high']/per_donor['n_cells']
per_donor.to_csv(out/'CRS_SenePy_per_donor.csv',index=False)
order=['Control_Young','Control_Aged','CRS_Young','CRS_Aged']
group=per_donor.groupby('group4',observed=True).agg(n_donors=('orig.ident','nunique'),mean_pct_high=('pct_high','mean'),median_pct_high=('pct_high','median'),mean_donor_score=('mean_score','mean')).reindex(order)
group.to_csv(out/'CRS_SenePy_group_summary.csv')
rows=[]
for y,g in [('Control_Young','Control_Aged'),('CRS_Young','CRS_Aged')]:
    v1=per_donor.loc[per_donor.group4==y,'pct_high'].to_numpy()
    v2=per_donor.loc[per_donor.group4==g,'pct_high'].to_numpy()
    rows.append(dict(comparison=f'{g} vs {y}',young_n=len(v1),aged_n=len(v2),young_mean_pct=v1.mean(),aged_mean_pct=v2.mean(),p_value=mannwhitneyu(v2,v1,alternative='two-sided',method='exact').pvalue))
tests=pd.DataFrame(rows)
tests['q_value']=multipletests(tests.p_value,method='fdr_bh')[1]
tests.to_csv(out/'CRS_SenePy_age_comparisons.csv',index=False)
# Cell-type score per independent donor; equal weight donors per group.
ct=(obs.groupby(['cell_type_final','orig.ident','group4'],observed=True).agg(n_cells=(score_col,'size'),median_score=(score_col,'median'),mean_score=(score_col,'mean'),pct_high=('SenePy_high_P90','mean')).reset_index())
ct['pct_high']*=100
ct.to_csv(out/'CRS_SenePy_celltype_per_donor.csv',index=False)
heat=ct.groupby(['cell_type_final','group4'],observed=True)['median_score'].mean().unstack('group4').reindex(columns=order)
heat.to_csv(out/'CRS_SenePy_celltype_heatmap_values.csv')
stat=[]
for y,g in [('Control_Young','Control_Aged'),('CRS_Young','CRS_Aged')]:
    for typ in sorted(ct.cell_type_final.unique()):
        sub=ct.loc[ct.cell_type_final==typ]
        a1=sub.loc[sub.group4==y,'median_score'].to_numpy()
        a2=sub.loc[sub.group4==g,'median_score'].to_numpy()
        p=mannwhitneyu(a2,a1,method='exact',alternative='two-sided').pvalue if min(len(a1),len(a2))>=2 else np.nan
        stat.append(dict(comparison=f'{g} vs {y}',cell_type=typ,young_n=len(a1),aged_n=len(a2),p_value=p,delta_mean_donor_median=a2.mean()-a1.mean() if len(a1) and len(a2) else np.nan))
st=pd.DataFrame(stat);st['q_value']=np.nan
for name in st.comparison.unique():
    mask=(st.comparison==name)&st.p_value.notna()
    st.loc[mask,'q_value']=multipletests(st.loc[mask,'p_value'],method='fdr_bh')[1]
st.to_csv(out/'CRS_SenePy_celltype_age_statistics.csv',index=False)
if a.rescore:
    adata.uns['SenePy_analysis']={'signature':'Human universal','package':'SenePy v1.0.1','signature_genes':734,'matched_signature_genes':733,'n_bins':25,'ctrl_size':50,'binarize':True,'importance':True,'global_p90':p90}
    adata.write_h5ad(out/'CRS_Control_72262_FULLGENES_SenePy_FINAL.h5ad',compression='gzip')
print(group.round(4).to_string())
print(tests.round(5).to_string(index=False))
print('Saved:',out)

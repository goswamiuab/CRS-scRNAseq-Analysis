
"""Human OE reanalysis.

.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import harmonypy as hm
from scipy import sparse


def parse_args():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--counts', type=Path, required=True, help='Merged donor-labeled raw counts h5ad')
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--annotation-map', type=Path, help='CSV: leiden_0.8,cell_type_final')
    p.add_argument('--scrublet-expected-doublet-rate', type=float, default=0.06,
                   )
    p.add_argument('--scrublet-seed', type=int, default=42,)
    p.add_argument('--hvg-flavor', default='seurat')
    p.add_argument('--harmony-theta', type=float, default=2.0)
    p.add_argument('--skip-scrublet', action='store_true')
    return p.parse_args()


def main():
    a=parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
    adata=sc.read_h5ad(a.counts)
    for c in ['sample','group','dataset']:
        if c not in adata.obs: raise ValueError(f'Missing required obs[{c!r}]')
    if adata.obs_names.has_duplicates or adata.var_names.has_duplicates:
        raise ValueError('Cell IDs and gene symbols must be unique before running.')
    if set(adata.obs['group'].astype(str).unique()) - {'Young','Aged'}:
        raise ValueError('Expected group Young/Aged; check manifest.')
    if not sparse.issparse(adata.X):
        raise ValueError(')
    data=adata.X.data
    if data.size and (np.min(data)<0 or not np.allclose(data,np.rint(data))):
        raise ValueError('X.')
    adata.obs['sample']=adata.obs['sample'].astype(str)
    adata.obs['group']=adata.obs['group'].astype(str)
    adata.obs['dataset']=adata.obs['dataset'].astype(str)
    if (adata.obs.groupby('sample')['group'].nunique()>1).any():
        raise ValueError('A donor assigned to multiple age groups.')
    if (adata.obs.groupby('sample')['dataset'].nunique()>1).any():
        raise ValueError('A donor assigned to multiple source datasets.')
    sc.pp.calculate_qc_metrics(adata, qc_vars=[], inplace=True, percent_top=None)
    mt=adata.var_names.str.upper().str.startswith('MT-')
    adata.obs['pct_counts_mt']=np.asarray(adata.X[:,mt].sum(axis=1)).ravel()/np.maximum(np.asarray(adata.X.sum(axis=1)).ravel(),1)*100
    # 'approximately' 200-8000 genes in Methods
    qc=(adata.obs.n_genes_by_counts.between(200,8000)
        & (adata.obs.total_counts>400) & (adata.obs.pct_counts_mt<10))
    adata.obs.assign(qc_pass=qc).to_csv(a.outdir/'preQC_cell_metrics.csv')
    adata=adata[qc].copy()
    sc.pp.filter_genes(adata,min_cells=3)
    print('Cells after QC and gene filtering:',adata.shape)
    if not a.skip_scrublet:
        # Scrublet estimates each donor separately.
        doublet=pd.Series(False,index=adata.obs_names,dtype=bool)
        scores=pd.Series(np.nan,index=adata.obs_names,dtype=float)
        for donor in adata.obs['sample'].unique():
            ids=adata.obs_names[adata.obs['sample']==donor]
            x=adata[ids,:].X
            model=scr.Scrublet(x, expected_doublet_rate=a.scrublet_expected_doublet_rate,
                               random_state=a.scrublet_seed)
            s,pred=model.scrub_doublets()
            scores.loc[ids]=s; doublet.loc[ids]=pred
        adata.obs['scrublet_score']=scores.reindex(adata.obs_names).to_numpy()
        adata.obs['predicted_doublet']=doublet.reindex(adata.obs_names).to_numpy()
        adata.obs[['sample','group','dataset','scrublet_score','predicted_doublet']].to_csv(a.outdir/'scrublet_predictions.csv')
        adata=adata[~adata.obs.predicted_doublet].copy()

    # Capture full gene expression (full-gene SenePy) 
    adata.layers['counts']=adata.X.copy()
    sc.pp.normalize_total(adata,target_sum=1e4)  
    sc.pp.log1p(adata)
    adata.raw=adata.copy() 
    sc.pp.highly_variable_genes(adata,n_top_genes=2000,flavor=a.hvg_flavor)
    adata_hvg=adata[:,adata.var.highly_variable].copy()
    sc.pp.scale(adata_hvg,max_value=10) 
    sc.tl.pca(adata_hvg,n_comps=50,svd_solver='arpack')
    corrected=hm.run_harmony(adata_hvg.obsm['X_pca'],adata_hvg.obs,'sample',theta=a.harmony_theta)
    z=corrected.Z_corr
    if z.shape==(50,adata_hvg.n_obs): z=z.T
    if z.shape!=(adata_hvg.n_obs,50)
    adata_hvg.obsm['X_pca_harmony']=z
    sc.pp.neighbors(adata_hvg,n_neighbors=15,n_pcs=30,use_rep='X_pca_harmony',metric='euclidean')
    sc.tl.umap(adata_hvg,random_state=42)  
    for resolution in [0.3,0.5,0.8,1.0]:
        sc.tl.leiden(adata_hvg,resolution=resolution,key_added=f'leiden_{resolution}',flavor='igraph',n_iterations=2)
    sc.tl.rank_genes_groups(adata_hvg,groupby='leiden_0.8',method='wilcoxon',corr_method='benjamini-hochberg')
    pd.DataFrame(adata_hvg.uns['rank_genes_groups']['names']).to_csv(a.outdir/'Leiden08_marker_gene_rankings.csv',index=False)
    if a.annotation_map:
        mapping=pd.read_csv(a.annotation_map,dtype=str)
        if not {'leiden_0.8','cell_type_final'}.issubset(mapping):
        lookup=mapping.set_index('leiden_0.8')['cell_type_final']
        adata_hvg.obs['cell_type_final']=adata_hvg.obs['leiden_0.8'].astype(str).map(lookup)
        if adata_hvg.obs['cell_type_final'].isna().any():
        clusters=sorted(adata_hvg.obs['leiden_0.8'].astype(str).unique(),key=lambda x:int(x))
        pd.DataFrame({'leiden_0.8':clusters,'cell_type_final':['VERIFY_MANUALLY']*len(clusters)}).to_csv(a.outdir/'annotation_map_TEMPLATE.csv',index=False)
        print('Cell annotation  assigned')
    #  full-feature AnnData, preserving counts score verification.
    adata.obs=adata.obs.join(adata_hvg.obs[[f'leiden_{r}' for r in (0.3,0.5,0.8,1.0)]],how='left')
    if 'cell_type_final' in adata_hvg.obs: adata.obs['cell_type_final']=adata_hvg.obs['cell_type_final']
    adata.obsm['X_umap']=adata_hvg.obsm['X_umap'].copy()
    adata.obsm['X_pca_harmony']=adata_hvg.obsm['X_pca_harmony'].copy()
    adata.uns['analysis_provenance']={
        'source':'Reconstructed Methods implementation; not the original code',
        'scrublet_expected_doublet_rate':a.scrublet_expected_doublet_rate,
        'hvg_flavor':a.hvg_flavor, 'harmony_theta':a.harmony_theta,
    path=a.outdir/'Human_OE_PostQC_Fullgenes_Reanalysis.h5ad'
    adata.write_h5ad(path,compression='gzip')
    pd.crosstab(adata.obs['sample'],adata.obs['group']).to_csv(a.outdir/'final_sample_cell_counts.csv')
    print('Saved:',path,'shape:',adata.shape)
    if adata.n_obs != 64249:
        print('')
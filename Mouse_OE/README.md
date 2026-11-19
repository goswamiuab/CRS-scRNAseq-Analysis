# Mouse olfactory epithelium — manuscript code (Seurat, SenMayo, SenePy)

## Scope
These scripts provide **standard, concise, executable downstream analysis** rather than every exploratory notebook or plotting instruction.

- `01_Mouse_OE_Seurat_SenMayo.R` reads the saved QC-passed Seurat object, joins **verified** final cell-type annotations by exact cell ID (annotation CSV required), runs UCell SenMayo scoring with a supplied **exact published 125-gene mouse-compatible signature**, and summarizes scores/burden at the mouse level.
- `02_Mouse_OE_SenePy_Validated_Results.py` reads the final **validated full-gene SenePy-scored** H5AD, checks the known global P90 and per-mouse results, and calculates donor/cell-type summaries and exact Mann–Whitney tests with BH adjustment.

## Confirmed original Seurat parameters
The saved Seurat `@commands` and final H5Seurat history establish: LogNormalize (scale.factor=10,000); 2,000 `vst` variable features; centered/scaled with scale.max=10; 50 PCA components (seed 42); PCA dims 1:20 for neighbors (k=20; Euclidean Annoy) and UMAP (30 neighbors; cosine, min.dist=0.3; seed 42); SNN clusters at 0.5 (29), 1.0 (45), 1.5 (57); final clusters correspond to 1.5 and 23 annotated populations. **No Harmony integration was verified** in the final Mouse OE pipeline.

## Limitations / reproducibility disclosures
1. The exact original pre-QC import, thresholds, doublet handling, and annotation mapping are not recovered. The final retained cells have nCount_RNA >=1000, nFeature_RNA 201–6000 and percent.mt <10%, but those observed boundaries do **not** establish the original commands. The Seurat script is started from the processed object and preserves its original reductions.
2. The Seurat RDS `E:/Public_GSE_data/mouse_data_output/Seurat_Mouse_YoungAged_filtered.rds` lacks `FinalCellType`. Export an **exact-ID** annotation CSV from `Mouse_fullgenes_zellkonverter.h5ad` via:

   ```python
   import scanpy as sc
   a=sc.read_h5ad(r'E:\Public_GSE_data\Mouse_data\Mice_data_senePy_analysis4young_4Aged\New folder\Mouse_fullgenes_zellkonverter.h5ad',backed='r')
   x=a.obs[['FinalCellType']].copy();x.index.name='cell_id';x.reset_index().to_csv('FinalCellType_by_cell.csv',index=False)
   a.file.close()
   ```
3. Provide the original, properly ortholog-mapped **SenMayo 125-gene CSV**, header `gene`. The script will not invent a signature. The exact SenMayo signature used for the reported manuscript must be the same one supplied at execution.
4. The final mouse SenePy analysis uses 102,550 cells × 33,502 genes, human-independent Mouse universal signature (635 genes; 582 represented), SenePy v1.0.1, 25 expression bins, 50 controls per bin, binarize=True, importance=True. A separate memory-efficient implementation was benchmarked against official score_hub on 200 identical cells with zero numerical difference. **That actual algorithm implementation is not included in this concise downstream package.** For a fully rerunnable score-from-counts repository, attach the original validated sparse-scoring notebook/script as well.
5. The older `Mouse_fullgenes_CURRENT_SenePy.h5ad` column `senepy_score` is **not** automatically accepted as validated. The downstream Python requires a file with `SenePy_universal_fullgenes`, checking global P90=0.9295258497, Young mean mouse burden ~0.4538%, Aged ~16.8175% and exact p=0.0285714. Do not rename an old score column to bypass checks.
6. Quantile P75 SenMayo and P90 SenePy are *analysis-specific operational high-score labels*, not proof of cellular senescence. Individual mice are the inferential replicates. The group comparison is unadjusted exact Mann–Whitney; cell-type p-values are BH-adjusted.

## Run

R:
```
Rscript 01_Mouse_OE_Seurat_SenMayo.R \
 "E:/Public_GSE_data/mouse_data_output/Seurat_Mouse_YoungAged_filtered.rds" \
 SenMayo_125_genes.csv results/SenMayo FinalCellType_by_cell.csv
```

Python:
```
python 02_Mouse_OE_SenePy_Validated_Results.py \
  "E:/Reviesed_with_comment response_10_04_2026/Final _Revised_10_09_2026/Mouse_OE_FINAL_PostQC_SenePy_Fullgenes.h5ad" \
  results/SenePy
```

Install R packages: `Seurat`, `SeuratObject`, `UCell`, `Matrix`; Python packages: `scanpy`, `scipy`, `pandas`, `numpy`, `statsmodels`.

Run scripts on your machine before public deposition, cite a versioned SenMayo signature file, and compare output tables to the exact manuscript panels.

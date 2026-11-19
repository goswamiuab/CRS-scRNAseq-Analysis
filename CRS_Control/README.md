# Human CRS/control scRNA-seq: compact manuscript analysis code

This supporting-code package provides a **Seurat/R + SenMayo UCell** workflow and a **Python/SenePy v1.0.1 donor-level analysis** for 16 human donors in four groups (Control Young/Aged, CRS Young/Aged). It intentionally excludes exploratory plotting and unrelated intermediate notebooks.

## Verified study population
- Original post-QC Seurat object: 81,608 cells / 27,204 RNA genes.
- `scDblFinder` classifications identified 9,346 predicted doublets; 72,262 predicted singlets retained from 16 donors (21 cell types).
- Full-gene SenePy: human universal signature 733 of 734 genes matched (99.86%), global P90=1.505203997154884, 7,227 high cells.
- SenePy v1.0.1 used 25 bins, 50 controls, binary expression, importance weighting. The official `score_hub` returns a list of scores in cell order.
- Inferential unit is the **donor**, never the cell. Within-control and within-CRS age comparisons are predefined exact two-sided Mann-Whitney tests, BH-adjusted across two comparisons. Cell-type medians are independently summarized per donor, BH adjusted per comparison across cell types.

## Files
- `01_CRS_Seurat_Doublet_SenMayo.R`: reuse existing `scDblFinder.class`, or calculate it from the original full-gene counts per sequencing library; then filter singlets, score a **provided, verified 125-gene human SenMayo signature** with UCell and generate donor-level tables. A signature file is required because the exact original SenMayo gene list/version is not present in the saved source records. This script is a standard documented workflow, not a reproduction of historical Seurat QC/integration/clustering.
- `02_CRS_SenePy_and_Donor_Analysis.py`: runs reproducible downstream analysis from the confirmed scored H5AD. `--rescore` optionally calls the official full-gene SenePy v1.0.1 function from the unscored full-gene sparse H5AD (requires a sufficiently large HPC memory allocation).

## Example commands
```bash
Rscript 01_CRS_Seurat_Doublet_SenMayo.R /path/crs.final.rds /path/senmayo_human_125_genes.txt ./results_R
python 02_CRS_SenePy_and_Donor_Analysis.py --input /path/CRS_Control_72262_FULLGENES_SenePy_FINAL.h5ad --out ./results_Python
# High-memory HPC only, for a new score computation:
python 02_CRS_SenePy_and_Donor_Analysis.py --input /path/CRS_Control_72262_Singlets_27204_Genes_PRE_SenePy.h5ad --out ./results_Python --rescore
```

## Dependencies
R: Seurat, SeuratObject, Matrix, UCell; optionally SingleCellExperiment + scDblFinder if doublets need to be called anew.
Python: scanpy, anndata, numpy, pandas, scipy, statsmodels; senepy==1.0.1 for `--rescore`.

## Reproducibility limits
- The original Seurat QC/integration command history for CRS has not been independently recovered here; no historical clustering or Harmony results are claimed to be reproduced by the compact R script.
- **SenMayo recomputation must be compared against the original saved manuscript scores**, because the originally used signature list and Seurat assay/layer scoring configuration were not verified during this reconstruction. No SenMayo numerical results are asserted by this package before execution.
- The **SenePy score-column analysis** is grounded in the successfully executed full-gene scoring run on Cheaha; run it against the user's actual final H5AD before public deposition. The full scorer is memory intensive.
- All figure formatting and exploratory visualization are excluded deliberately. These scripts can be shared as targeted reproducibility supplements, not a complete FASTQ-to-figure workflow.

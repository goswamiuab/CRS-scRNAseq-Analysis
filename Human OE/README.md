# Human OE single-cell RNA-seq: manuscript-supporting analysis code

## Status and scope
This concise package supports **human olfactory epithelium** senescence analyses associated with public datasets **GSE139522 (Durante)** and **GSE184117 (Oliva)**. It is not proof of exact end-to-end reproduction of the original 64,249-cell atlas. Do not claim the scripts were run as written until checked against final results.

### Files
1. `01_Human_OE_SenMayo_UCell.R`: Recompute SenMayo (published 125-gene input) in a user-supplied **annotated Seurat RDS** with verified columns (`sample`,`group`,`cell_type_final`). Generates global P75 and donor summaries. **No matching Human OE Seurat RDS nor original gene-set file has been verified.** If the Human OE data are currently stored only in an H5AD, transfer data/metadata to a Seurat object using a documented workflow before running. Do not confuse this with the original CRS RDS.
2. `02_Human_OE_SenePy_Donor_Analysis.py`: Recompute global P90 and donor-level summaries from the *saved* `senepy_score` in the final Human OE 64,249-cell H5AD. Source-adjusted regression requires a verified sample-to-GEO manifest (`sample,dataset`). No origin is inferred from age.
3. `03_Human_OE_SenePy_FullGene_Rescore_OPTIONAL.py`: Optional **independent** full-count rescoring using the official SenePy v1.0.1 function. High memory required. Returns a per-cell list, not AnnData. **Not validated against the historical Human OE score values.** Do not report these as original scores until independently reproduced and compared.
4. `04_Human_OE_QC_Harmony_Leiden_TEMPLATE.py`: Reconstructed Python preprocessing template from stated Human OE Methods, **NOT verified original executable analysis**. It should NOT be described as the exact original QC/annotation code. Parameters stated in manuscript: >400 UMI, roughly 200–8,000 genes, <10% mito, Scrublet, 2,000 HVGs, 50 PCs, Harmony donor batch, 30 Harmony PCs for neighbors, 15 neighbors Euclidean, Leiden 0.8 final. Verify actual original notebook and donor-to-GEO manifest for exact reproduction.

## Typical run (Python downstream)
```bash
python 02_Human_OE_SenePy_Donor_Analysis.py \
  --input /data/user/goswamis/Human_OE_Reanalysis_2026/Human_OE_64249_SenePy.h5ad \
  --manifest verified_human_oe_geo_manifest.csv \
  --output Human_OE_Manuscript_SenePy
```
Without `--manifest`, donor- and cell-type-level tables still run, but adjusted regression is skipped and a CSV template is created.

## Typical run (SenMayo R)
```bash
Rscript 01_Human_OE_SenMayo_UCell.R /verified/path/Human_OE_annotated.rds /verified/path/SenMayo_125_genes.txt Human_OE_SenMayo_results
```
This recalculates UCell on the given expression object; exact equivalence depends on matching gene nomenclature, cell universe, assay/data and historical UCell settings.

## Statistical interpretation
Donor is the biological unit. Global P75 (SenMayo) and P90 (SenePy) are *dataset-specific operational thresholds*, not formal diagnoses of senescence. Source origin is partially confounded with age in the public Human OE cohort, and within-GSE139522 there is only one young donor. P values are not based on treating cells as independent replicates.

## Before Figshare deposition
Verify input paths, original executed notebooks for QC/Harmony/Scrublet, original Human OE SenMayo signature/method and SenePy count input, GEO source manifest, run output consistency, software versions and data sharing permissions. Python scripts have been syntax-checked only.

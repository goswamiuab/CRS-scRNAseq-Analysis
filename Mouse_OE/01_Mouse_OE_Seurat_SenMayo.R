
# Mouse OE: Analysis
Seurat, SeuratObject, UCell, Matrix, ggplot2

# Rscript 01_Mouse_OE_Seurat_SenMayo.R \
#   "/mouse_data_output/Seurat_Mouse_YoungAged.rds" \
#   "SenMayo_125_genes.csv" "mouse_oe_results" \
#   "FinalCellType_by_cell.csv"
#
# SenMayo gene list CSV: gene with the EXACT 125 published
# mouse-ortholog gene symbols used for this study.

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3 || length(args) > 4) {
  stop("Usage: Rscript 01_Mouse_OE_Seurat_SenMayo.R <input.rds> <senmayo_genes.csv> <output_dir> [annotation.csv]")
}
stopifnot(requireNamespace("Seurat", quietly=TRUE), requireNamespace("SeuratObject", quietly=TRUE),
          requireNamespace("UCell", quietly=TRUE), requireNamespace("Matrix", quietly=TRUE))
input <- args[1]; signature_file <- args[2]; out <- args[3]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
mouse <- readRDS(input)
stopifnot(inherits(mouse,"Seurat"), ncol(mouse)==102550, nrow(mouse[["RNA"]])==33502)
stopifnot(all(c("sample_id","age_group") %in% colnames(mouse[[]])))
stopifnot(length(unique(as.character(mouse$sample_id)))==8)

# final annotation via explicit, cell-ID-validated
if (!"FinalCellType" %in% colnames(mouse[[]])) {
  if (length(args)<4) stop("FinalCellType absent from RDS; supply annotation CSV (cell_id, FinalCellType).")
  ann <- read.csv(args[4], stringsAsFactors=FALSE, check.names=FALSE)
  stopifnot(all(c("cell_id","FinalCellType") %in% names(ann)), !anyDuplicated(ann$cell_id))
  if (!all(colnames(mouse) %in% ann$cell_id))
  label <- ann$FinalCellType[match(colnames(mouse),ann$cell_id)]
  stopifnot(!anyNA(label))
  mouse$FinalCellType <- label
}


# LogNormalize (scale factor 10000), vst 2000 HVG, scaled max 10,
# PCA npcs 50 (seed 42), FindNeighbors PCA dims 1:20 k=20 (Euclidean/Annoy),
# UMAP PCA dims 1:20 n.neighbors=30 cosine min.dist=.3 (seed 42).
# H5Seurat RNA_snn clustering at resolution 1.5 (57 clusters).

qc <- mouse[[]][,intersect(c("sample_id","age_group","nCount_RNA","nFeature_RNA","percent.mt","FinalCellType"), colnames(mouse[[]])),drop=FALSE]
qc$cell_id <- rownames(qc)
write.csv(qc,file.path(out,"Mouse_OE_QC_and_annotations.csv"),row.names=FALSE)
write.csv(as.data.frame(table(sample_id=mouse$sample_id,group=mouse$age_group)),
          file.path(out,"Mouse_OE_Sample_Counts.csv"),row.names=FALSE)


# mouse <- Seurat::NormalizeData(mouse, normalization.method="LogNormalize", scale.factor=10000)
# mouse <- Seurat::FindVariableFeatures(mouse,selection.method="vst",nfeatures=2000)
# mouse <- Seurat::ScaleData(mouse,do.center=TRUE,do.scale=TRUE,scale.max=10)
# mouse <- Seurat::RunPCA(mouse,npcs=50,seed.use=42)
# mouse <- Seurat::FindNeighbors(mouse,reduction="pca",dims=1:20,k.param=20,
#                               nn.method="annoy",n.trees=50,annoy.metric="euclidean")
# mouse <- Seurat::RunUMAP(mouse,reduction="pca",dims=1:20,n.neighbors=30,
#                          metric="cosine",min.dist=.3,spread=1,seed.use=42)
# mouse <- Seurat::FindClusters(mouse,graph.name="RNA_snn",resolution=c(.5,1,1.5),
#                               algorithm=1,random.seed=0)

signature <- read.csv(signature_file, stringsAsFactors=FALSE, check.names=FALSE)
if (!"gene" %in% names(signature)) stop("Signature CSV requires column: gene")
genes <- unique(trimws(as.character(signature$gene)))
genes <- genes[nzchar(genes)]
if (length(genes)!=125) stop(sprintf("Expected exact 125-gene SenMayo set; got %d unique genes",length(genes)))
seen <- intersect(genes, rownames(mouse[["RNA"]]))
write.csv(data.frame(gene=genes,present=genes %in% seen),file.path(out,"SenMayo_gene_coverage.csv"),row.names=FALSE)
message("SenMayo matched genes: ", length(seen),"/125")
if (length(seen)<30)

# Apply UCell in Seurat
Seurat::DefaultAssay(mouse) <- "RNA"
mouse <- UCell::AddModuleScore_UCell(mouse, features=list(SenMayo=seen), assay="RNA")
if (!"SenMayo_UCell" %in% colnames(mouse[[]]))

m <- mouse[[]]
m$cell_id <- rownames(m)
m$age_group <- as.character(m$age_group)
m$sample_id <- as.character(m$sample_id)
m$FinalCellType <- as.character(m$FinalCellType)
cutoff <- as.numeric(stats::quantile(m$SenMayo_UCell,.75,na.rm=TRUE,names=FALSE))
m$SenMayo_high_P75 <- m$SenMayo_UCell >= cutoff
write.csv(m[,c("cell_id","sample_id","age_group","FinalCellType","SenMayo_UCell","SenMayo_high_P75")],
          file.path(out,"Mouse_OE_SenMayo_PerCell.csv"),row.names=FALSE)

by_mouse <- do.call(rbind,lapply(split(m,m$sample_id),function(z) data.frame(
  sample_id=as.character(z$sample_id[1]),age_group=as.character(z$age_group[1]),
  n_cells=nrow(z),n_high=sum(z$SenMayo_high_P75),
  pct_high=100*mean(z$SenMayo_high_P75),mean_score=mean(z$SenMayo_UCell),
  median_score=median(z$SenMayo_UCell))))
write.csv(by_mouse,file.path(out,"Mouse_OE_SenMayo_PerMouse.csv"),row.names=FALSE)

# Exact mouse-level Wilcoxon test (4 young, 4 aged).
y <- by_mouse$pct_high[by_mouse$age_group=="Young"]
a <- by_mouse$pct_high[by_mouse$age_group=="Aged"]
stopifnot(length(y)==4,length(a)==4)

test <- suppressWarnings(stats::wilcox.test(a,y,alternative="two.sided",exact=TRUE))
write.csv(data.frame(metric="SenMayo_high_pct",threshold_P75=cutoff,
                     Young_mean=mean(y),Aged_mean=mean(a),p_value=test$p.value),
          file.path(out,"Mouse_OE_SenMayo_Statistics.csv"),row.names=FALSE)
writeLines(c(paste("Seurat",as.character(packageVersion("Seurat"))),
             paste("UCell",as.character(packageVersion("UCell"))),
             paste("R",getRversion()),
             paste("SenMayo P75",signif(cutoff,10))),file.path(out,"versions.txt"))
message("Saved SenMayo results: ", normalizePath(out))

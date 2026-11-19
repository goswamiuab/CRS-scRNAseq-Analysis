
# CRS/control Seurat processing

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)


crs <- readRDS(rds_path)
stopifnot(inherits(crs, 'Seurat'), 'RNA' %in% SeuratObject::Assays(crs))

provenance <- capture.output({
  cat('Original dimensions (features x cells):', paste(dim(crs), collapse=' x '), '\n')
  cat('Assays:', paste(SeuratObject::Assays(crs), collapse=', '), '\n')
  cat('Default assay:', SeuratObject::DefaultAssay(crs), '\n')
  cat('Reductions:', paste(SeuratObject::Reductions(crs), collapse=', '), '\n')
  cat('SCT 0.7 cluster column present:', 'SCT_snn_res.0.7' %in% names(crs@meta.data), '\n')
  for (cn in c('NormalizeData.RNA','FindVariableFeatures.RNA','RunPCA.SCT','RunUMAP.SCT.pca')) {
    if (cn %in% names(crs@commands)) {
      cat('\n', cn, '\n', sep='')
      p <- crs@commands[[cn]]@params
      for (key in intersect(c('normalization.method','scale.factor','selection.method','nfeatures','npcs','dims','n.neighbors','metric','min.dist','seed.use'),names(p))) {
        cat(key, ': ', paste(p[[key]], collapse=', '), '\n', sep='')
      }
    }
  }
})

meta <- crs@meta.data
needed <- c('orig.ident','sample','cell_type')
if (!all(needed %in% names(meta))) stop('Missing metadata: ', paste(setdiff(needed,names(meta)),collapse=', '))
# scDblFinder 
if (!'scDblFinder.class' %in% names(meta)) {
  suppressPackageStartupMessages({library(SingleCellExperiment); library(scDblFinder)})
  counts <- SeuratObject::LayerData(crs, assay = 'RNA', layer = 'counts')
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  sce$sample <- as.character(crs$sample)
  set.seed(42)
  sce <- scDblFinder::scDblFinder(sce, samples = 'sample')
  crs$scDblFinder.class <- as.character(sce$scDblFinder.class)
}
print(table(crs$scDblFinder.class))
singlets <- subset(crs, subset = scDblFinder.class == 'singlet')
id <- as.character(singlets$orig.ident)
group <- ifelse(startsWith(id,'C4_'),'Control_Young',ifelse(startsWith(id,'C6_'),'Control_Aged',ifelse(startsWith(id,'T4_'),'CRS_Young',ifelse(startsWith(id,'T6_'),'CRS_Aged',NA_character_))))
stopifnot(!anyNA(group))
singlets$group4 <- group
singlets$cell_type_final <- as.character(singlets$cell_type)
# Signature verified set of human HGNC symbols.
signature <- unique(trimws(readLines(signature_path, warn = FALSE)))
signature <- signature[nzchar(signature)]
cat('Signature supplied:', length(signature), 'genes\n')
if (length(signature) != 125L) 
# UCell scores from RNA expression 
if (!requireNamespace('UCell', quietly=TRUE)) stop('Install UCell to generate SenMayo scores')
DefaultAssay(singlets) <- 'RNA'
# UCell function may read RNA data layer
if (length(grep('^data', SeuratObject::Layers(singlets[['RNA']]))) == 0L) {
  singlets <- NormalizeData(singlets, normalization.method='LogNormalize', scale.factor=10000, verbose=FALSE)
}
singlets <- UCell::AddModuleScore_UCell(singlets, features = list(SenMayo=signature), name=NULL)
score_name <- 'SenMayo_UCell'
if (!score_name %in% colnames(singlets@meta.data)) stop('UCell output missing ', score_name, '; inspect columns')
p75 <- unname(quantile(singlets@meta.data[[score_name]], .75, na.rm=FALSE))
singlets$SenMayo_high_P75 <- singlets@meta.data[[score_name]] >= p75
obs <- singlets@meta.data
# Explicit tidy summary
splits <- split(seq_len(nrow(obs)), interaction(obs$orig.ident,obs$group4, drop=TRUE))
tab <- do.call(rbind,lapply(splits,function(ix){ data.frame(donor=as.character(obs$orig.ident[ix[1]]),group=as.character(obs$group4[ix[1]]),n_cells=length(ix),n_high=sum(obs$SenMayo_high_P75[ix]),pct_high=100*mean(obs$SenMayo_high_P75[ix]),mean_score=mean(obs[[score_name]][ix]),median_score=median(obs[[score_name]][ix])) }))
rownames(tab)<-NULL
write.csv(tab,file.path(outdir,'CRS_SenMayo_per_donor.csv'),row.names=FALSE)
write.csv(data.frame(cell_id=colnames(singlets),score=obs[[score_name]],high=obs$SenMayo_high_P75),file.path(outdir,'CRS_SenMayo_per_cell.csv'),row.names=FALSE)
writeLines(sprintf('SenMayo P75 = %.12g; %d singlets, %d donors', p75, ncol(singlets),length(unique(obs$orig.ident))),file.path(outdir,'SenMayo_run_summary.txt'))
cat('Finished SenMayo; output:',outdir,'\n')

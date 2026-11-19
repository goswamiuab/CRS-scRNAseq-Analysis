
# Human OE SenMayo UCell scoring 

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3)
input <- args[1]; genesfile <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
obj <- readRDS(input)
stopifnot(inherits(obj,"Seurat"))
meta <- obj@meta.data
needed <- c("sample","group","cell_type_final")
if (!all(needed %in% colnames(meta))) 
stopifnot(!anyNA(meta[,needed]))
genes <- unique(trimws(readLines(genesfile, warn=FALSE)))
genes <- genes[nzchar(genes)]
if (length(genes)!=125L) stop("Expected exactly 125 unique SenMayo genes: got ",length(genes))
# Assay names
if (!("RNA" %in% SeuratObject::Assays(obj)))
DefaultAssay(obj) <- "RNA"
present <- intersect(genes, rownames(obj[["RNA"]]))
write.csv(data.frame(gene=genes, present=genes %in% present), file.path(outdir,"SenMayo_gene_coverage.csv"), row.names=FALSE)
if (length(present)<50) stop("Low SenMayo gene coverage; verify gene naming and input assay")
obj <- UCell::AddModuleScore_UCell(obj, features=list(SenMayo=present), name=NULL)
scorecol <- "SenMayo_UCell"
if (!scorecol %in% colnames(obj@meta.data))
score <- as.numeric(obj@meta.data[[scorecol]])
p75 <- unname(quantile(score,0.75,na.rm=FALSE))
df <- data.frame(cell_id=colnames(obj), donor=as.character(obj$sample), age_group=as.character(obj$group), cell_type=as.character(obj$cell_type_final), SenMayo_UCell=score, SenMayo_high=score>=p75, stringsAsFactors=FALSE)
stopifnot(all(df$age_group %in% c("Young","Aged")))
write.csv(df, file.path(outdir,"Human_OE_SenMayo_per_cell.csv"), row.names=FALSE)
by_donor <- aggregate(cbind(n_cells=rep(1L,nrow(df)), n_high=as.integer(df$SenMayo_high), score_sum=df$SenMayo_UCell), by=list(donor=df$donor,group=df$age_group), FUN=sum)
by_donor$pct_high <- 100*by_donor$n_high/by_donor$n_cells
by_donor$mean_score <- by_donor$score_sum/by_donor$n_cells
by_donor$score_sum <- NULL
write.csv(by_donor,file.path(outdir,"Human_OE_SenMayo_by_donor.csv"),row.names=FALSE)
writeLines(c(paste0("cells=",nrow(df)),paste0("donors=",length(unique(df$donor))),paste0("input_signature_genes=",length(genes)),paste0("matched_genes=",length(present)),paste0("global_p75=",format(p75,digits=16)),paste0("Seurat=",as.character(packageVersion("Seurat"))),paste0("UCell=",as.character(packageVersion("UCell")))), file.path(outdir,"SenMayo_run_provenance.txt"))
cat("SenMayo completed; global P75:",p75,"\n")

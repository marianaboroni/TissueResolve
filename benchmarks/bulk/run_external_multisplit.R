#!/usr/bin/env Rscript
# Run MuSiC + BisqueRNA across ALL (split, scenario) datasets exported by
# export_external_multisplit.py.  Writes per-(split,scenario,method) predictions
# and a status table.  Honest: records executed/failed/skipped, never fabricates.
suppressWarnings({
  .libPaths(c(normalizePath(file.path("benchmarks","envs","r_lib"), mustWork=FALSE), .libPaths()))
})
BASE <- file.path("benchmarks","outputs","holdout_bulk","external_inputs_multisplit")
SCENARIOS <- c("balanced","imbalanced","rare","similar_subtypes","missing_population")
REF_OF <- function(scen) if (scen=="missing_population") "reference_missing" else "reference"
have <- function(p) requireNamespace(p, quietly=TRUE)
read_mtx <- function(f) as.matrix(read.table(f, header=TRUE, sep="\t", row.names=1, check.names=FALSE))

run_music <- function(sc, meta, bulk, cts) {
  suppressMessages(library(Biobase)); suppressMessages(library(MuSiC))
  if ("bulk.eset" %in% names(formals(MuSiC::music_prop))) {
    pheno <- data.frame(cellType=meta$cellType, SubjectName=meta$SubjectName, row.names=meta$cell_id)
    sc.eset <- ExpressionSet(assayData=sc, phenoData=AnnotatedDataFrame(pheno[colnames(sc),,drop=FALSE]))
    est <- music_prop(bulk.eset=ExpressionSet(assayData=bulk), sc.eset=sc.eset,
                      clusters="cellType", samples="SubjectName", select.ct=cts, verbose=FALSE)
  } else {
    suppressMessages(library(SingleCellExperiment))
    sce <- SingleCellExperiment(assays=list(counts=sc),
            colData=DataFrame(cellType=meta$cellType, SubjectName=meta$SubjectName, row.names=meta$cell_id))
    est <- music_prop(bulk.mtx=bulk, sc.sce=sce, clusters="cellType",
                      samples="SubjectName", select.ct=cts, verbose=FALSE)
  }
  p <- est$Est.prop.weighted; p / rowSums(p)
}
run_bisque <- function(sc, meta, bulk, cts) {
  suppressMessages(library(Biobase)); suppressMessages(library(BisqueRNA))
  pheno <- data.frame(cellType=meta$cellType, SubjectName=meta$SubjectName, row.names=meta$cell_id)
  sc.eset <- ExpressionSet(assayData=sc, phenoData=AnnotatedDataFrame(pheno[colnames(sc),,drop=FALSE]))
  out <- ReferenceBasedDecomposition(ExpressionSet(assayData=bulk), sc.eset, markers=NULL,
            cell.types="cellType", subject.names="SubjectName", use.overlap=FALSE)
  p <- t(out$bulk.props); p / rowSums(p)
}

splits <- list.dirs(BASE, recursive=FALSE, full.names=FALSE)
splits <- splits[grepl("^split[0-9]+$", splits)]
rows <- list()
for (sp in splits) {
  for (scen in SCENARIOS) {
    refdir <- file.path(BASE, sp, REF_OF(scen))
    sc_f <- file.path(refdir, "reference_counts_genes_by_cells.tsv")
    bulk_f <- file.path(BASE, sp, scen, "bulk_counts_genes_by_samples.tsv")
    if (!file.exists(sc_f) || !file.exists(bulk_f)) next
    sc <- read_mtx(sc_f); bulk <- read_mtx(bulk_f)
    meta <- read.table(file.path(refdir,"reference_cell_metadata.tsv"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
    shared <- intersect(rownames(sc), rownames(bulk))
    sc <- sc[shared,,drop=FALSE]; bulk <- bulk[shared,,drop=FALSE]; cts <- unique(meta$cellType)
    for (m in c("MuSiC","BisqueRNA")) {
      key <- paste(sp, scen, m, sep="__")
      if (!have(m) || !have("Biobase")) {
        rows[[key]] <- data.frame(split=sp, scenario=scen, method=m, status="skipped",
          version="", runtime_seconds=0, stringsAsFactors=FALSE); next
      }
      t0 <- Sys.time()
      res <- tryCatch(list(p=if(m=="MuSiC") run_music(sc,meta,bulk,cts) else run_bisque(sc,meta,bulk,cts),
                           ver=as.character(packageVersion(m))), error=function(e) e)
      rt <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      if (inherits(res,"error")) {
        rows[[key]] <- data.frame(split=sp, scenario=scen, method=m, status="failed",
          version="", runtime_seconds=round(rt,2), stringsAsFactors=FALSE)
        cat(sprintf("%-34s FAILED %s\n", key, conditionMessage(res))); next
      }
      write.table(round(res$p,6), file.path(BASE, sp, scen, paste0(m,"_pred.tsv")),
                  sep="\t", quote=FALSE, col.names=NA)
      rows[[key]] <- data.frame(split=sp, scenario=scen, method=m, status="executed",
        version=res$ver, runtime_seconds=round(rt,2), stringsAsFactors=FALSE)
      cat(sprintf("%-34s executed %.1fs\n", key, rt))
    }
  }
}
write.table(do.call(rbind, rows), file.path(BASE, "external_multisplit_status.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("\nWrote", file.path(BASE,"external_multisplit_status.tsv"), "\n")

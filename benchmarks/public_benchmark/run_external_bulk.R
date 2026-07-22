#!/usr/bin/env Rscript
# Public benchmark — run external bulk tools (MuSiC, BisqueRNA, BayesPrism) on the SAME exported
# single-cell reference + mixtures as TissueResolve. Writes pred_<Method>.tsv (samples x cellTypes)
# per scenario + a status table. Never fabricates: failures are recorded, not filled in.
suppressWarnings({
  for (l in c("benchmarks/envs/Rlib","benchmarks/envs/r_lib"))
    .libPaths(c(normalizePath(l, mustWork=FALSE), .libPaths()))
})
args <- commandArgs(trailingOnly=TRUE)
BASE <- if (length(args)>=1) args[1] else "benchmarks/results/public_benchmark/bulk/external_inputs"
DATASETS <- if (length(args)>=2) strsplit(args[2],",")[[1]] else c("breast","lung")
SCEN <- c("base","low_depth","reduced_overlap")
# BayesPrism is very slow (min-scale per run); opt-in only. Off by default -> recorded deferred_runtime.
BAYESPRISM_SCEN <- if (nzchar(Sys.getenv("RUN_BAYESPRISM"))) c("base") else character(0)
.ms <- Sys.getenv("MAX_SHARED_GENES"); MAX_SHARED <- as.integer(if (nzchar(.ms)) .ms else "8000")
have <- function(p) requireNamespace(p, quietly=TRUE)
`%||%` <- function(a,b) if (is.null(a)||length(a)==0) b else a
read_mtx <- function(f) as.matrix(read.table(f, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
status <- list()

run_music <- function(sc, meta, bulk, cts) {
  suppressMessages(library(Biobase)); suppressMessages(library(MuSiC))
  if ("bulk.eset" %in% names(formals(MuSiC::music_prop))) {
    ph <- data.frame(cellType=meta$cellType, SubjectName=meta$SubjectName, row.names=meta$cell_id)
    sce <- ExpressionSet(assayData=sc, phenoData=AnnotatedDataFrame(ph[colnames(sc),,drop=FALSE]))
    est <- music_prop(bulk.eset=ExpressionSet(assayData=bulk), sc.eset=sce, clusters="cellType",
                      samples="SubjectName", select.ct=cts, verbose=FALSE)
  } else {
    suppressMessages(library(SingleCellExperiment))
    sce <- SingleCellExperiment(assays=list(counts=sc),
      colData=DataFrame(cellType=meta$cellType, SubjectName=meta$SubjectName, row.names=meta$cell_id))
    est <- music_prop(bulk.mtx=bulk, sc.sce=sce, clusters="cellType", samples="SubjectName",
                      select.ct=cts, verbose=FALSE)
  }
  p <- est$Est.prop.weighted; p / rowSums(p)
}

run_bisque <- function(sc, meta, bulk, cts) {
  suppressMessages(library(Biobase)); suppressMessages(library(BisqueRNA))
  ph <- data.frame(cellType=meta$cellType, SubjectName=meta$SubjectName, row.names=meta$cell_id)
  sce <- ExpressionSet(assayData=sc, phenoData=AnnotatedDataFrame(ph[colnames(sc),,drop=FALSE]))
  out <- ReferenceBasedDecomposition(ExpressionSet(assayData=bulk), sce, markers=NULL,
           cell.types="cellType", subject.names="SubjectName", use.overlap=FALSE)
  p <- t(out$bulk.props); p / rowSums(p)
}

run_bayesprism <- function(sc, meta, bulk, cts) {
  suppressMessages(library(BayesPrism))
  scc <- t(sc); bkk <- t(bulk)                        # cells x genes ; samples x genes
  pr <- new.prism(reference=scc, mixture=bkk, input.type="count.matrix",
                  cell.type.labels=meta$cellType, cell.state.labels=meta$cellType, key=NULL)
  res <- run.prism(pr, n.cores=2)
  theta <- get.fraction(res, which.theta="final", state.or.type="type")
  theta / rowSums(theta)
}

for (ds in DATASETS) {
  refdir <- file.path(BASE, ds, "reference")
  if (!file.exists(file.path(refdir, "reference_counts_genes_by_cells.tsv"))) {
    cat("skip dataset (no export):", ds, "\n"); next }
  sc0 <- read_mtx(file.path(refdir, "reference_counts_genes_by_cells.tsv"))
  meta <- read.table(file.path(refdir, "reference_cell_metadata.tsv"), header=TRUE, sep="\t",
                     stringsAsFactors=FALSE)
  cts <- unique(meta$cellType)
  for (scen in SCEN) {
    bf <- file.path(BASE, ds, scen, "bulk_counts_genes_by_samples.tsv")
    if (!file.exists(bf)) next
    bulk <- read_mtx(bf)
    shared <- intersect(rownames(sc0), rownames(bulk))
    sc <- sc0[shared,,drop=FALSE]; bk <- bulk[shared,,drop=FALSE]
    methods <- c("MuSiC","BisqueRNA")
    if (scen %in% BAYESPRISM_SCEN) methods <- c(methods, "BayesPrism")
    for (m in methods) {
      key <- paste(ds, scen, m, sep="__")
      pkg <- if (m=="BisqueRNA") "BisqueRNA" else m
      if (m %in% c("MuSiC","BisqueRNA") && length(shared) > MAX_SHARED) {
        status[[key]] <- list(dataset=ds,scenario=scen,method=m,status="deferred",version="",
          runtime_seconds=0,n_shared=length(shared),
          error=sprintf("deferred (runtime): %d shared genes > MAX_SHARED=%d", length(shared), MAX_SHARED))
        cat(key,"deferred (runtime, large gene set)\n"); next }
      if (!have(pkg)) { status[[key]] <- list(dataset=ds,scenario=scen,method=m,status="failed_install",
        version="",runtime_seconds=0,n_shared=length(shared),error="not loadable")
        cat(key,"failed_install\n"); next }
      t0 <- Sys.time()
      res <- tryCatch({
        p <- switch(m, MuSiC=run_music(sc,meta,bk,cts), BisqueRNA=run_bisque(sc,meta,bk,cts),
                    BayesPrism=run_bayesprism(sc,meta,bk,cts))
        list(p=p, ver=as.character(packageVersion(pkg)))
      }, error=function(e) e)
      rt <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      if (inherits(res,"error")) { status[[key]] <- list(dataset=ds,scenario=scen,method=m,
        status="failed_run",version=as.character(packageVersion(pkg)),runtime_seconds=rt,
        n_shared=length(shared),error=conditionMessage(res)); cat(key,"FAILED:",conditionMessage(res),"\n"); next }
      predf <- file.path(BASE, ds, scen, paste0("pred_", m, ".tsv"))
      write.table(round(res$p,6), predf, sep="\t", quote=FALSE, col.names=NA)
      status[[key]] <- list(dataset=ds,scenario=scen,method=m,status="executed",version=res$ver,
        runtime_seconds=rt,n_shared=length(shared),error="")
      cat(sprintf("%s executed v%s %.1fs shared=%d\n", key, res$ver, rt, length(shared)))
    }
  }
}
rows <- do.call(rbind, lapply(status, function(s) data.frame(dataset=s$dataset,scenario=s$scenario,
  method=s$method,status=s$status,version=s$version,runtime_seconds=round(s$runtime_seconds,2),
  n_shared=s$n_shared,error=substr(s$error %||% "",1,200),stringsAsFactors=FALSE)))
write.table(rows, file.path(BASE,"external_method_status.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("\nWrote", file.path(BASE,"external_method_status.tsv"), "\n")

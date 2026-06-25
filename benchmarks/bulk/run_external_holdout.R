#!/usr/bin/env Rscript
# Run external bulk methods (MuSiC, BisqueRNA) on the held-out-donor scenarios
# exported by export_external_inputs.py.  Same reference + bulk as TissueResolve.
# Records status/version/runtime per (scenario, method); never fabricates results.
suppressWarnings({
  lib <- file.path("benchmarks", "envs", "r_lib")
  .libPaths(c(normalizePath(lib, mustWork = FALSE), .libPaths()))
})
BASE <- file.path("benchmarks", "outputs", "holdout_bulk", "external_inputs")
SCENARIOS <- c("balanced", "imbalanced", "rare", "similar_subtypes", "missing_population")
REF_OF <- list(balanced="reference", imbalanced="reference", rare="reference",
               similar_subtypes="reference", missing_population="reference_missing")
status <- list()

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
have <- function(p) requireNamespace(p, quietly = TRUE)

read_mtx <- function(f) as.matrix(read.table(f, header=TRUE, sep="\t", row.names=1, check.names=FALSE))

run_music <- function(sc, meta, bulk, cts) {
  suppressMessages(library(Biobase)); suppressMessages(library(MuSiC))
  has_eset <- "bulk.eset" %in% names(formals(MuSiC::music_prop))
  if (has_eset) {
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

for (scen in SCENARIOS) {
  refdir <- file.path(BASE, REF_OF[[scen]])
  sc_f   <- file.path(refdir, "reference_counts_genes_by_cells.tsv")
  meta_f <- file.path(refdir, "reference_cell_metadata.tsv")
  bulk_f <- file.path(BASE, scen, "bulk_counts_genes_by_samples.tsv")
  sc <- read_mtx(sc_f); bulk <- read_mtx(bulk_f)
  meta <- read.table(meta_f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  shared <- intersect(rownames(sc), rownames(bulk))
  sc <- sc[shared,,drop=FALSE]; bulk <- bulk[shared,,drop=FALSE]
  cts <- unique(meta$cellType)

  for (m in c("MuSiC", "BisqueRNA")) {
    key <- paste(scen, m, sep="__")
    pkg_ok <- have(m) && have("Biobase")
    if (!pkg_ok) {
      status[[key]] <- list(scenario=scen, method=m, status="skipped",
        version="", runtime_seconds=0, n_shared=length(shared),
        error=paste(m, "not installed/loadable"))
      cat(sprintf("%-30s skipped (not installed)\n", key)); next
    }
    t0 <- Sys.time()
    res <- tryCatch({
      p <- if (m == "MuSiC") run_music(sc, meta, bulk, cts) else run_bisque(sc, meta, bulk, cts)
      list(p=p, ver=as.character(packageVersion(m)))
    }, error=function(e) e)
    rt <- as.numeric(difftime(Sys.time(), t0, units="secs"))
    if (inherits(res, "error")) {
      status[[key]] <- list(scenario=scen, method=m, status="failed",
        version=if(have(m)) as.character(packageVersion(m)) else "", runtime_seconds=rt,
        n_shared=length(shared), error=conditionMessage(res))
      cat(sprintf("%-30s FAILED: %s\n", key, conditionMessage(res))); next
    }
    predf <- file.path(BASE, scen, paste0(m, "_pred.tsv"))
    write.table(round(res$p, 6), predf, sep="\t", quote=FALSE, col.names=NA)
    status[[key]] <- list(scenario=scen, method=m, status="executed", version=res$ver,
      runtime_seconds=rt, n_shared=length(shared), error="", output=predf)
    cat(sprintf("%-30s executed v%s %.1fs shared=%d\n", key, res$ver, rt, length(shared)))
  }
}

# write a flat status table
rows <- do.call(rbind, lapply(status, function(s)
  data.frame(scenario=s$scenario, method=s$method, status=s$status, version=s$version,
             runtime_seconds=round(s$runtime_seconds,2), n_shared=s$n_shared,
             error=substr(s$error %||% "", 1, 200), stringsAsFactors=FALSE)))
write.table(rows, file.path(BASE, "external_method_status.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("\nWrote", file.path(BASE, "external_method_status.tsv"), "\n")

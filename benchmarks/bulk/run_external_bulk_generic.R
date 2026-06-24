#!/usr/bin/env Rscript
# Generic external bulk runner (MuSiC, BisqueRNA) over a base dir of scenario subdirs
# exported by export_external_bulk_generic.py. One shared "reference/" dir; each
# scenario subdir has bulk_counts_genes_by_samples.tsv. Writes <method>_pred.tsv per
# scenario + external_method_status.tsv. Records status/version/runtime; never fabricates.
#
# Usage: Rscript benchmarks/bulk/run_external_bulk_generic.R <BASE_DIR>
suppressWarnings({
  lib <- file.path("benchmarks", "envs", "r_lib")
  .libPaths(c(normalizePath(lib, mustWork = FALSE), .libPaths()))
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: run_external_bulk_generic.R <BASE_DIR>")
BASE <- args[[1]]
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

refdir <- file.path(BASE, "reference")
sc <- read_mtx(file.path(refdir, "reference_counts_genes_by_cells.tsv"))
meta <- read.table(file.path(refdir, "reference_cell_metadata.tsv"),
                   header=TRUE, sep="\t", stringsAsFactors=FALSE)
scen_dirs <- list.dirs(BASE, full.names=FALSE, recursive=FALSE)
scen_dirs <- scen_dirs[file.exists(file.path(BASE, scen_dirs, "bulk_counts_genes_by_samples.tsv"))]
status <- list()
for (scen in scen_dirs) {
  bulk <- read_mtx(file.path(BASE, scen, "bulk_counts_genes_by_samples.tsv"))
  shared <- intersect(rownames(sc), rownames(bulk))
  sc_s <- sc[shared,,drop=FALSE]; bulk_s <- bulk[shared,,drop=FALSE]
  cts <- unique(meta$cellType)
  for (m in c("MuSiC", "BisqueRNA")) {
    key <- paste(scen, m); t0 <- Sys.time()
    if (!have(m)) { status[[key]] <- list(scenario=scen, method=m, status="skipped",
                                          version="NA", runtime_seconds=NA, n_shared=length(shared),
                                          error="package not installed"); next }
    res <- tryCatch({
      p <- if (m == "MuSiC") run_music(sc_s, meta, bulk_s, cts) else run_bisque(sc_s, meta, bulk_s, cts)
      list(p=p, ver=as.character(packageVersion(m)))
    }, error=function(e) list(err=conditionMessage(e)))
    rt <- as.numeric(difftime(Sys.time(), t0, units="secs"))
    if (!is.null(res$err)) {
      status[[key]] <- list(scenario=scen, method=m, status="failed", version="NA",
                            runtime_seconds=round(rt,2), n_shared=length(shared), error=res$err)
      cat(scen, m, "FAILED:", res$err, "\n"); next
    }
    write.table(res$p, file.path(BASE, scen, paste0(m, "_pred.tsv")),
                sep="\t", quote=FALSE, col.names=NA)
    status[[key]] <- list(scenario=scen, method=m, status="executed", version=res$ver,
                          runtime_seconds=round(rt,2), n_shared=length(shared), error="")
    cat(scen, m, "ok", round(rt,1), "s\n")
  }
}
rows <- do.call(rbind, lapply(status, function(s) data.frame(
  scenario=s$scenario, method=s$method, status=s$status, version=s$version,
  runtime_seconds=s$runtime_seconds, n_shared=s$n_shared, error=s$error)))
write.table(rows, file.path(BASE, "external_method_status.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("\nWrote", file.path(BASE, "external_method_status.tsv"), "\n")

#!/usr/bin/env Rscript
# Real BisqueRNA bulk deconvolution runner.
# Uses harmonised prepared inputs (same reference + pseudobulk as TissueResolve).
# Fails gracefully (records skipped/failed) without aborting the benchmark.
.libPaths(c(file.path("benchmarks","envs","Rlib"), .libPaths()))
source(file.path("benchmarks","bulk","methods","_external_common.R"))
method <- "BisqueRNA"; t0 <- Sys.time()
logf <- file.path("benchmarks","logs","run_BisqueRNA.log")
dir.create(dirname(logf), recursive=TRUE, showWarnings=FALSE)

fail <- function(status, msg) {
  write_meta(method, list(method=method, version="", executed=FALSE, imported=FALSE,
    exported_only=FALSE, status=status,
    runtime_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")),
    input_normalization_used="counts", reference_level_used="fine",
    command_run="run_bisque.R", warnings="", error_message=msg, output_path=""))
  cat(method, ":", status, "-", msg, "\n"); quit(status=0)
}

if (!have_pkg("BisqueRNA") || !have_pkg("Biobase"))
  fail("skipped", "BisqueRNA/Biobase not installed (see benchmarks/envs/install_external_tools.sh)")

prep <- file.path("benchmarks","outputs","prepared_inputs")
sc_counts_f <- file.path(prep,"reference","reference_counts_genes_by_cells.tsv")
sc_meta_f   <- file.path(prep,"reference","reference_cell_metadata.tsv")
bulk_f      <- file.path(prep,"bulk","bulk_counts_genes_by_samples.tsv")
for (f in c(sc_counts_f, sc_meta_f, bulk_f))
  if (!file.exists(f)) fail("failed", paste("missing prepared input:", f))

res <- tryCatch({
  suppressMessages(library(Biobase)); suppressMessages(library(BisqueRNA))
  sc <- as.matrix(read.table(sc_counts_f, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  meta <- read.table(sc_meta_f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  bulk <- as.matrix(read.table(bulk_f, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  shared <- intersect(rownames(sc), rownames(bulk))
  if (length(shared) < 50) stop(paste("too few shared genes:", length(shared)))
  sc <- sc[shared, , drop=FALSE]; bulk <- bulk[shared, , drop=FALSE]

  pheno <- data.frame(cellType=meta$cellType, SubjectName=meta$SubjectName,
                      row.names=meta$cell_id)
  sc.eset <- ExpressionSet(assayData=sc,
              phenoData=AnnotatedDataFrame(pheno[colnames(sc), , drop=FALSE]))
  bulk.eset <- ExpressionSet(assayData=bulk)

  out <- ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL,
            cell.types="cellType", subject.names="SubjectName", use.overlap=FALSE)
  props <- t(out$bulk.props)             # samples x cell types
  props <- props / rowSums(props)         # renormalise rows to sum to 1
  list(props=props, version=as.character(packageVersion("BisqueRNA")),
       n_shared=length(shared))
}, error=function(e) e)

if (inherits(res, "error")) fail("failed", conditionMessage(res))

predf <- out_pred(method)
write.table(round(res$props, 6), predf, sep="\t", quote=FALSE, col.names=NA)
write_meta(method, list(method=method, version=res$version, executed=TRUE,
  imported=FALSE, exported_only=FALSE, status="executed",
  runtime_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")),
  input_normalization_used="counts (rows renormalised to sum 1)",
  reference_level_used="fine",
  command_run="BisqueRNA::ReferenceBasedDecomposition(use.overlap=FALSE)",
  warnings=paste("shared genes:", res$n_shared), error_message="",
  output_path=predf))
cat(method, ": executed -> ", predf, "\n")

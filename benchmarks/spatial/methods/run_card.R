#!/usr/bin/env Rscript
# Real CARD spatial deconvolution runner.
# Uses harmonised prepared inputs (same single-cell reference as TissueResolve
# + the exported Visium counts/coordinates). Real Visium has no ground truth,
# so downstream comparison is concordance/structure, not accuracy.
# Fails gracefully (records skipped/failed) without aborting the benchmark.
.libPaths(c(file.path("benchmarks","envs","Rlib"), .libPaths()))
source(file.path("benchmarks","spatial","methods","_external_common.R"))
method <- "CARD"; t0 <- Sys.time()
logf <- file.path("benchmarks","logs","run_CARD.log")
dir.create(dirname(logf), recursive=TRUE, showWarnings=FALSE)

fail <- function(status, msg) {
  write_meta(method, list(method=method, version="", executed=FALSE, imported=FALSE,
    exported_only=FALSE, status=status,
    runtime_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")),
    input_normalization_used="counts", reference_level_used="fine",
    command_run="run_card.R", warnings="", error_message=msg, output_path=""))
  cat(method, ":", status, "-", msg, "\n"); quit(status=0)
}

if (!have_pkg("CARD")) fail("skipped",
  "CARD not installed (see benchmarks/envs/install_external_tools.sh)")

prep <- file.path("benchmarks","outputs","prepared_inputs")
sc_counts_f <- file.path(prep,"reference","reference_counts_genes_by_cells.tsv")
sc_meta_f   <- file.path(prep,"reference","reference_cell_metadata.tsv")
sp_counts_f <- file.path(prep,"spatial","spatial_counts_genes_by_spots.tsv")
sp_loc_f    <- file.path(prep,"spatial","spatial_coordinates.tsv")
for (f in c(sc_counts_f, sc_meta_f, sp_counts_f, sp_loc_f))
  if (!file.exists(f)) fail("failed", paste("missing prepared input:", f))

res <- tryCatch({
  suppressMessages(library(CARD))
  sc <- as.matrix(read.table(sc_counts_f, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  meta <- read.table(sc_meta_f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  rownames(meta) <- meta$cell_id; meta$sampleInfo <- meta$SubjectName
  sp <- as.matrix(read.table(sp_counts_f, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  loc <- read.table(sp_loc_f, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
  loc <- loc[colnames(sp), c("x","y"), drop=FALSE]
  cts <- unique(meta$cellType)

  obj <- createCARDObject(sc_count=as(sc, "sparseMatrix"), sc_meta=meta[colnames(sc), ],
            spatial_count=as(sp, "sparseMatrix"), spatial_location=loc,
            ct.varname="cellType", ct.select=cts, sample.varname="sampleInfo",
            minCountGene=0, minCountSpot=0)
  obj <- CARD_deconvolution(CARD_object=obj)
  props <- as.matrix(obj@Proportion_CARD)   # spots x cell types
  props <- props / rowSums(props)
  list(props=props, version=as.character(packageVersion("CARD")),
       n_spots=nrow(props))
}, error=function(e) e)

if (inherits(res, "error")) fail("failed", conditionMessage(res))

predf <- out_pred(method)
write.table(round(res$props, 6), predf, sep="\t", quote=FALSE, col.names=NA)
write_meta(method, list(method=method, version=res$version, executed=TRUE,
  imported=FALSE, exported_only=FALSE, status="executed",
  runtime_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")),
  input_normalization_used="counts (rows renormalised to sum 1)",
  reference_level_used="fine",
  command_run="CARD::createCARDObject + CARD_deconvolution",
  warnings=paste("spots:", res$n_spots, "; real Visium has no ground truth"),
  error_message="", output_path=predf))
cat(method, ": executed -> ", predf, "\n")

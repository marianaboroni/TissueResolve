#!/usr/bin/env Rscript
# RCTD (spacexr) on the synthetic-spatial scenario exported by export_spatial_external.py.
# Outputs spot×celltype proportions (full mode weights, row-normalised). Honest status.
suppressWarnings(.libPaths(c(normalizePath(file.path("benchmarks","envs","r_lib"),
                                           mustWork=FALSE), .libPaths())))
BASE <- file.path("benchmarks","outputs","spatial_real_visium","external_inputs")
t0 <- Sys.time()
status <- function(st, ver, msg, predf="") {
  write.table(data.frame(method="RCTD", status=st, version=ver,
    runtime_seconds=round(as.numeric(difftime(Sys.time(),t0,units="secs")),1),
    error=substr(msg,1,200), output=predf),
    file.path(BASE,"RCTD_status.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  cat("RCTD:", st, ver, msg, "\n")
}
if (!requireNamespace("spacexr", quietly=TRUE)) { status("skipped","","spacexr not installed"); quit(status=0) }

res <- tryCatch({
  suppressMessages(library(spacexr)); suppressMessages(library(Matrix))
  sc <- as.matrix(read.table(file.path(BASE,"reference","reference_counts_genes_by_cells.tsv"),
                             header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  meta <- read.table(file.path(BASE,"reference","reference_cell_metadata.tsv"),
                     header=TRUE, sep="\t", stringsAsFactors=FALSE)
  spots <- as.matrix(read.table(file.path(BASE,"spot_counts_genes_by_spots.tsv"),
                                header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  coords <- read.table(file.path(BASE,"spot_coords.tsv"),
                       header=TRUE, sep="\t", stringsAsFactors=FALSE)
  ct <- factor(meta$cellType); names(ct) <- meta$cell_id
  # drop cell types with < 25 cells (RCTD requirement)
  keep_ct <- names(which(table(ct) >= 25))
  cells <- meta$cell_id[meta$cellType %in% keep_ct]
  sc <- sc[, cells, drop=FALSE]; ct <- droplevels(ct[cells])
  reference <- Reference(sc, ct)
  cm <- data.frame(row.names=coords$spot, x=coords$row, y=coords$col)
  puck <- SpatialRNA(cm, spots, colSums(spots))
  rctd <- create.RCTD(puck, reference, max_cores=1, CELL_MIN_INSTANCE=25)
  rctd <- run.RCTD(rctd, doublet_mode="full")
  w <- as.matrix(rctd@results$weights)            # spots × celltypes
  w <- sweep(w, 1, rowSums(w), "/")
  list(props=w, ver=as.character(packageVersion("spacexr")), nct=ncol(w))
}, error=function(e) e)

if (inherits(res,"error")) { status("failed","", conditionMessage(res)); quit(status=0) }
predf <- file.path(BASE, "RCTD_pred.tsv")
write.table(round(res$props,6), predf, sep="\t", quote=FALSE, col.names=NA)
status("executed", res$ver, paste(res$nct,"cell types"), predf)

#!/usr/bin/env Rscript
# RCTD across all seed<s>/external_inputs/ produced by run_spatial_multiseed.py.
suppressWarnings(.libPaths(c(normalizePath(file.path("benchmarks","envs","r_lib"),
                                           mustWork=FALSE), .libPaths())))
BASE <- file.path("benchmarks","outputs","spatial_multiseed")
have <- requireNamespace("spacexr", quietly=TRUE)
read_mtx <- function(f) as.matrix(read.table(f, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
seeds <- list.dirs(BASE, recursive=FALSE, full.names=FALSE)
seeds <- seeds[grepl("^seed[0-9]+$", seeds)]
rows <- list()
for (sd in seeds) {
  ei <- file.path(BASE, sd, "external_inputs")
  predf <- file.path(BASE, sd, "RCTD_pred.tsv")
  t0 <- Sys.time()
  if (!have) { rows[[sd]] <- data.frame(seed=sd, method="RCTD", status="skipped",
               version="", runtime_seconds=0); next }
  res <- tryCatch({
    suppressMessages(library(spacexr))
    sc <- read_mtx(file.path(ei,"reference","reference_counts_genes_by_cells.tsv"))
    meta <- read.table(file.path(ei,"reference","reference_cell_metadata.tsv"),
                       header=TRUE, sep="\t", stringsAsFactors=FALSE)
    spots <- read_mtx(file.path(ei,"spatial_spot_counts.tsv"))
    coords <- read.table(file.path(ei,"spatial_coords.tsv"), header=TRUE, sep="\t")
    ct <- factor(meta$cellType); names(ct) <- meta$cell_id
    keep_ct <- names(which(table(ct) >= 25))
    cells <- meta$cell_id[meta$cellType %in% keep_ct]
    sc <- sc[, cells, drop=FALSE]; ct <- droplevels(ct[cells])
    reference <- Reference(sc, ct)
    cm <- data.frame(row.names=coords$spot, x=coords$row, y=coords$col)
    rctd <- create.RCTD(SpatialRNA(cm, spots, colSums(spots)), reference,
                        max_cores=1, CELL_MIN_INSTANCE=25)
    rctd <- run.RCTD(rctd, doublet_mode="full")
    w <- as.matrix(rctd@results$weights); w <- sweep(w, 1, rowSums(w), "/")
    list(props=w, ver=as.character(packageVersion("spacexr")))
  }, error=function(e) e)
  rt <- as.numeric(difftime(Sys.time(), t0, units="secs"))
  if (inherits(res,"error")) {
    rows[[sd]] <- data.frame(seed=sd, method="RCTD", status="failed", version="",
                             runtime_seconds=round(rt,1))
    cat(sd, "RCTD FAILED:", conditionMessage(res), "\n"); next
  }
  write.table(round(res$props,6), predf, sep="\t", quote=FALSE, col.names=NA)
  rows[[sd]] <- data.frame(seed=sd, method="RCTD", status="executed",
                           version=res$ver, runtime_seconds=round(rt,1))
  cat(sd, "RCTD executed", round(rt,0), "s\n")
}
write.table(do.call(rbind, rows), file.path(BASE,"RCTD_multiseed_status.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote RCTD_multiseed_status.tsv\n")

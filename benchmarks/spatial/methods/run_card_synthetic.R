#!/usr/bin/env Rscript
# CARD spatial deconvolution on the synthetic-spatial scenario (ground truth).
# Reads the exported external_inputs/ layout (same reference + spots + truth as
# TissueResolve) and writes CARD_pred.tsv + CARD_status.tsv so that
# score_external_synthetic.py scores it with the same metrics. Fails gracefully.
.libPaths(c(file.path("benchmarks","envs","Rlib"), .libPaths()))
method <- "CARD"; t0 <- Sys.time()
base <- file.path("benchmarks","outputs","spatial_synthetic","external_inputs")

status <- function(st, ver, msg, predf="") {
  df <- data.frame(method=method, status=st, version=ver,
                   runtime_seconds=round(as.numeric(difftime(Sys.time(),t0,units="secs")),1),
                   error=substr(as.character(msg),1,200), output=predf,
                   stringsAsFactors=FALSE)
  write.table(df, file.path(base, "CARD_status.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  cat(method, ":", st, "-", as.character(msg), "\n")
}

if (!requireNamespace("CARD", quietly=TRUE)) { status("skipped","","CARD not installed"); quit(status=0) }

sc_counts_f <- file.path(base,"reference","reference_counts_genes_by_cells.tsv")
sc_meta_f   <- file.path(base,"reference","reference_cell_metadata.tsv")
sp_counts_f <- file.path(base,"spatial","spot_counts_genes_by_spots.tsv")
sp_loc_f    <- file.path(base,"spatial","spot_coords.tsv")
for (f in c(sc_counts_f, sc_meta_f, sp_counts_f, sp_loc_f))
  if (!file.exists(f)) { status("failed","",paste("missing input:",f)); quit(status=0) }

res <- tryCatch({
  suppressMessages(library(CARD))
  sc <- as.matrix(read.table(sc_counts_f, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  meta <- read.table(sc_meta_f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  rownames(meta) <- meta$cell_id
  meta$sampleInfo <- "synthetic"
  sp <- as.matrix(read.table(sp_counts_f, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  loc <- read.table(sp_loc_f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  rownames(loc) <- as.character(loc$spot)
  # CARD expects a spatial_location data.frame with columns x,y indexed by spot.
  locxy <- data.frame(x=loc$row, y=loc$col, row.names=rownames(loc))
  locxy <- locxy[colnames(sp), , drop=FALSE]
  cts <- unique(meta$cellType)

  obj <- createCARDObject(sc_count=as(sc, "sparseMatrix"), sc_meta=meta[colnames(sc), ],
            spatial_count=as(sp, "sparseMatrix"), spatial_location=locxy,
            ct.varname="cellType", ct.select=cts, sample.varname="sampleInfo",
            minCountGene=0, minCountSpot=0)
  obj <- CARD_deconvolution(CARD_object=obj)
  props <- as.matrix(obj@Proportion_CARD)        # spots x cell types
  props <- props / rowSums(props)
  list(props=props, version=as.character(packageVersion("CARD")), n=nrow(props))
}, error=function(e) e)

if (inherits(res, "error")) { status("failed","",conditionMessage(res)); quit(status=0) }

predf <- file.path(base, "CARD_pred.tsv")
write.table(round(res$props, 6), predf, sep="\t", quote=FALSE, col.names=NA)
status("executed", res$version, paste(res$n, "spots"), predf)

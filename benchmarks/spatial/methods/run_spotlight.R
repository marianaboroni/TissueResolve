#!/usr/bin/env Rscript
# Real SPOTlight spatial deconvolution runner.
# Uses harmonised prepared inputs (same single-cell reference as TissueResolve
# + exported Visium counts). Marker genes per cell type are derived from the
# reference (top mean-expression genes). Real Visium has no ground truth, so
# downstream comparison is concordance/structure, not accuracy.
# Fails gracefully (records skipped/failed) without aborting the benchmark.
.libPaths(c(file.path("benchmarks","envs","Rlib"), .libPaths()))
source(file.path("benchmarks","spatial","methods","_external_common.R"))
method <- "SPOTlight"; t0 <- Sys.time()
logf <- file.path("benchmarks","logs","run_SPOTlight.log")
dir.create(dirname(logf), recursive=TRUE, showWarnings=FALSE)

fail <- function(status, msg) {
  write_meta(method, list(method=method, version="", executed=FALSE, imported=FALSE,
    exported_only=FALSE, status=status,
    runtime_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")),
    input_normalization_used="counts", reference_level_used="fine",
    command_run="run_spotlight.R", warnings="", error_message=msg, output_path=""))
  cat(method, ":", status, "-", msg, "\n"); quit(status=0)
}

if (!have_pkg("SPOTlight")) fail("skipped",
  "SPOTlight not installed (see benchmarks/envs/install_external_tools.sh)")

prep <- file.path("benchmarks","outputs","prepared_inputs")
sc_counts_f <- file.path(prep,"reference","reference_counts_genes_by_cells.tsv")
sc_meta_f   <- file.path(prep,"reference","reference_cell_metadata.tsv")
sp_counts_f <- file.path(prep,"spatial","spatial_counts_genes_by_spots.tsv")
for (f in c(sc_counts_f, sc_meta_f, sp_counts_f))
  if (!file.exists(f)) fail("failed", paste("missing prepared input:", f))

res <- tryCatch({
  suppressMessages(library(SPOTlight))
  sc <- as.matrix(read.table(sc_counts_f, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  meta <- read.table(sc_meta_f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  sp <- as.matrix(read.table(sp_counts_f, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  shared <- intersect(rownames(sc), rownames(sp))
  if (length(shared) < 50) stop(paste("too few shared genes:", length(shared)))
  sc <- sc[shared, , drop=FALSE]; sp <- sp[shared, , drop=FALSE]
  labels <- meta$cellType[match(colnames(sc), meta$cell_id)]

  # marker genes per cell type from the reference (top mean-expression),
  # as the mgs data.frame SPOTlight expects (gene, cluster, weight)
  cts <- unique(labels)
  mgs_list <- lapply(cts, function(ct) {
    m <- rowMeans(sc[, labels == ct, drop=FALSE])
    g <- names(sort(m, decreasing=TRUE))[1:min(50, length(m))]
    data.frame(gene=g, cluster=ct, weight=as.numeric(m[g]), stringsAsFactors=FALSE)
  })
  mgs <- do.call(rbind, mgs_list)

  out <- SPOTlight(x=sc, y=sp, groups=as.character(labels), mgs=mgs,
                   gene_id="gene", group_id="cluster", weight_id="weight")
  mat <- if (is.list(out) && !is.null(out$mat)) out$mat else out[[1]]
  props <- as.matrix(mat)                  # spots x cell types
  props <- props / rowSums(props)
  list(props=props, version=as.character(packageVersion("SPOTlight")),
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
  command_run="SPOTlight::SPOTlight(x, y, groups, mgs)",
  warnings=paste("spots:", res$n_spots, "; markers from reference; no ground truth"),
  error_message="", output_path=predf))
cat(method, ": executed -> ", predf, "\n")

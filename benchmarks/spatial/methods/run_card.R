#!/usr/bin/env Rscript
source(file.path("benchmarks","spatial","methods","_external_common.R"))
method <- "CARD"; t0 <- Sys.time()
if (!have_pkg("CARD")) {
  write_meta(method, list(method=method, version="", executed=FALSE, imported=FALSE,
    exported_only=FALSE, status="skipped", runtime_seconds=0,
    input_normalization_used="counts", reference_level_used="fine",
    command_run="run_card.R", warnings="",
    error_message="CARD not installed (see benchmarks/envs/install_external_tools.sh)",
    output_path=""))
  cat("CARD: skipped (not installed)\n"); quit(status=0)
}
write_meta(method, list(method=method, version=as.character(packageVersion("CARD")),
  executed=FALSE, imported=FALSE, exported_only=FALSE, status="not_implemented_here",
  runtime_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")),
  input_normalization_used="counts", reference_level_used="fine",
  command_run="run_card.R", warnings="execution stub: installed but runner not wired",
  error_message="", output_path=""))
cat("CARD: installed; runner stub\n")

#!/usr/bin/env Rscript
source(file.path("benchmarks","bulk","methods","_external_common.R"))
method <- "BayesPrism"; t0 <- Sys.time()
if (!have_pkg("BayesPrism")) {
  write_meta(method, list(method=method, version="", executed=FALSE, imported=FALSE,
    exported_only=FALSE, status="skipped", runtime_seconds=0,
    input_normalization_used="counts", reference_level_used="fine",
    command_run="run_bayesprism.R", warnings="",
    error_message="BayesPrism not installed (see benchmarks/envs/install_external_tools.sh)",
    output_path=""))
  cat("BayesPrism: skipped (not installed)\n"); quit(status=0)
}
# Real run would go here, reading prepared_inputs and writing out_pred(method).
write_meta(method, list(method=method, version=as.character(packageVersion("BayesPrism")),
  executed=FALSE, imported=FALSE, exported_only=FALSE, status="not_implemented_here",
  runtime_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")),
  input_normalization_used="counts", reference_level_used="fine",
  command_run="run_bayesprism.R", warnings="execution stub: installed but runner not wired",
  error_message="", output_path=""))
cat("BayesPrism: installed; runner stub (wire BayesPrism::music_prop here)\n")

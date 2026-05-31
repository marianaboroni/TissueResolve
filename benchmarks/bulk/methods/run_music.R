#!/usr/bin/env Rscript
source(file.path("benchmarks","bulk","methods","_external_common.R"))
method <- "MuSiC"; t0 <- Sys.time()
if (!have_pkg("MuSiC")) {
  write_meta(method, list(method=method, version="", executed=FALSE, imported=FALSE,
    exported_only=FALSE, status="skipped", runtime_seconds=0,
    input_normalization_used="counts", reference_level_used="fine",
    command_run="run_music.R", warnings="",
    error_message="MuSiC not installed (see benchmarks/envs/install_external_tools.sh)",
    output_path=""))
  cat("MuSiC: skipped (not installed)\n"); quit(status=0)
}
# Real run would go here, reading prepared_inputs and writing out_pred(method).
write_meta(method, list(method=method, version=as.character(packageVersion("MuSiC")),
  executed=FALSE, imported=FALSE, exported_only=FALSE, status="not_implemented_here",
  runtime_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")),
  input_normalization_used="counts", reference_level_used="fine",
  command_run="run_music.R", warnings="execution stub: installed but runner not wired",
  error_message="", output_path=""))
cat("MuSiC: installed; runner stub (wire MuSiC::music_prop here)\n")

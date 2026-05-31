# Common helpers for external R bulk runners. Fail gracefully, write metadata.
suppressWarnings(suppressMessages({}))
prep_dir <- function() file.path("benchmarks","outputs","prepared_inputs")
out_pred <- function(method) {
  d <- file.path("benchmarks","outputs","bulk","predictions"); dir.create(d, recursive=TRUE, showWarnings=FALSE)
  file.path(d, paste0(method, ".tsv"))
}
write_meta <- function(method, fields) {
  d <- file.path("benchmarks","outputs","bulk","method_metadata"); dir.create(d, recursive=TRUE, showWarnings=FALSE)
  con <- file(file.path(d, paste0(method, ".json")), "w")
  cat("{\n", file=con)
  keys <- names(fields)
  for (i in seq_along(fields)) {
    v <- fields[[i]]
    vs <- if (is.logical(v)) tolower(as.character(v)) else if (is.numeric(v)) as.character(v) else paste0('"', gsub('"','\\\\"', v), '"')
    cat(sprintf('  "%s": %s%s\n', keys[i], vs, if (i<length(fields)) "," else ""), file=con)
  }
  cat("}\n", file=con); close(con)
}
have_pkg <- function(p) requireNamespace(p, quietly=TRUE)

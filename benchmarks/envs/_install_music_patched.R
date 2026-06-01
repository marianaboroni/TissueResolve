# Install MuSiC on R 4.1.2 by relaxing the TOAST version constraint.
# MuSiC's TOAST (>= 1.10.1) requirement needs Bioconductor 3.15 / R >= 4.2;
# TOAST is only used by MuSiC2, not classic music_prop. We relax it to the
# available TOAST 1.8.3 so the real music_prop can run. Documented workaround.
.libPaths(c("benchmarks/envs/Rlib", .libPaths()))
options(repos = "https://cloud.r-project.org")
td <- tempfile(); dir.create(td)
tb <- file.path(td, "music.tar.gz")
download.file("https://github.com/xuranw/MuSiC/archive/refs/heads/master.tar.gz",
              tb, quiet = TRUE)
untar(tb, exdir = td)
pkg <- list.files(td, pattern = "^MuSiC", full.names = TRUE)
pkg <- pkg[dir.exists(pkg)][1]
desc <- file.path(pkg, "DESCRIPTION")
d <- readLines(desc)
d <- gsub("TOAST \\(>= [0-9.]+\\)", "TOAST (>= 1.8.3)", d)
writeLines(d, desc)
cat("patched TOAST constraint in", desc, "\n")
install.packages(pkg, lib = "benchmarks/envs/Rlib", repos = NULL, type = "source")
ok <- requireNamespace("MuSiC", quietly = TRUE)
cat("MuSiC installed=", ok, "\n")
if (ok) {
  cat("VERSION=", as.character(packageVersion("MuSiC")), "\n")
  cat("eset_api=", "bulk.eset" %in% names(formals(MuSiC::music_prop)), "\n")
}

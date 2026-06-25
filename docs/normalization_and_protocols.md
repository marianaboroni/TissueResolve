# Normalization and protocol awareness

Deconvolution accuracy depends on the input scale and on protocol/library
compatibility between the reference and the query. TissueResolve and the
benchmark harness make these explicit rather than assuming them.

## Normalization detection

`benchmarks/shared/normalization.py` (and the benchmark runners) classify an
input matrix as one of: `counts`, `cpm`, `log_normalized`, `scaled`, or
`unknown`, using integer-ness, value range, row sums, and zero fraction. For
AnnData, raw counts are preferred from `layers['counts']` or `.raw` before `X`.

Each benchmark run records:

- the detected normalization status of every input,
- the normalization each method received and whether a conversion was applied,
- whether raw counts were available,
- a warning when a method receives non-ideal input (e.g. a counts-based method
  fed normalized values).

Conversions available: `convert_counts_to_cpm`, `convert_counts_to_log_cpm`.
Nothing is converted silently — the decision is written to the report.

## Protocol and library detection

`benchmarks/shared/protocol_detection.py` detects the reference protocol/library
type (scRNA 3'/5', single-nucleus, SMART-seq, mixed, unknown) and the query
protocol (bulk RNA-seq, Visium fresh-frozen/FFPE, unknown) from `obs`/`uns`/`var`
metadata (`assay`, `suspension_type`, `library_type`, feature types, …) and the
file structure, with a confidence score.

If the protocol is unknown, TissueResolve does **not** guess silently: it uses
conservative defaults, issues a CAUTION, and lets you override with
`--reference-protocol` / `--query-protocol` (benchmark) or the protocol-aware
gene weighting in the core tool.

The benchmark compares protocol-aware vs protocol-naive runs and writes
`protocol_detection.tsv`, `protocol_compatibility.tsv`, and
`library_type_detection.tsv`.

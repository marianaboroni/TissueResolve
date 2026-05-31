# Input Formats

## Single-cell reference (`.h5ad`)

TissueResolve accepts a single-cell or single-nucleus reference in the
`AnnData` `.h5ad` format.

Required elements:

- `adata.obs` contains a cell-type label column such as `cell_type`
- `adata.var` contains gene identifiers
- raw counts or expression matrices for reference building

The reference may also be a saved `ReferenceSignature` directory.

## Bulk counts table

Bulk query data should be a table with genes in rows and samples in columns.
Common formats:

- TSV with a gene identifier column and sample columns
- CSV with a header row and gene names in the first column

The tool automatically inspects the file and infers the query type.

## Visium input

Supported spatial query formats:

- a Visium `.h5ad` file with spatial coordinates and count data
- a Space Ranger-style folder with `filtered_feature_bc_matrix` and `spatial/`

## Optional marker gene file

A text file with one gene name per line may be supplied for marker selection
in spatial workflows.

## Gene identifiers

Gene-name harmonization is critical. Common identifier styles are:

- gene symbols (e.g. `TP53`)
- Ensembl IDs (e.g. `ENSG00000141510`)

If your reference and query use different conventions, map them before running
TissueResolve.

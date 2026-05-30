"""
I/O routines for TissueResolve.

Submodules (implemented in Stages 3–4):

``bulk``
    Read bulk count matrices from TSV, CSV, or h5ad.

``spatial``
    Load 10x Visium data from SpaceRanger output directories or h5ad files.
    Save spatial deconvolution results to h5ad.

``reference``
    Load single-cell / single-nucleus reference data from h5ad or CSV.

``validation``
    Gene-overlap checks and array coordinate validation.
"""

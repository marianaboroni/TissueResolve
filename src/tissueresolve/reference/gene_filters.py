"""
Protocol-aware gene filter lists for TissueResolve.

All gene lists are loaded lazily on first access and cached per-instance.
There are NO module-level mutable singletons — instantiate a fresh
:class:`GeneFilterSet` per pipeline run.

Gene lists bundled with TissueResolve
--------------------------------------
``intronic_dominant_{genome}.tsv``
    Genes where ≥ 30 % of expected signal originates from intronic reads
    in nuclear preparations.  Source: GTEx exon/intron split counts.
    Relevant when the reference is snRNA-seq.

``dissociation_stress_{genome}.tsv``
    Genes artificially up-regulated by enzymatic cell dissociation during
    scRNA-seq preparation (FOS, JUN, HSP70 family).
    Source: van den Brink et al. 2017 (Nat Methods); Machado et al. 2021.
    Relevant when the reference is scRNA-seq.

``length_biased_3prime_{genome}.tsv``
    Genes > 5 000 bp that are systematically under-detected by 3′-end UMI
    capture relative to full-length sequencing.
    Source: hg38/mm10 transcript length annotations.

``safe_universal_{genome}.tsv``
    Genes that pass all three filters above.  Fallback when protocol is
    unknown.

``hypervariable_inflammatory.tsv``
    Hypervariable inflammatory genes (IFN/TNF response, stress genes).
    Not genome-specific.

``protein_coding_{genome}.tsv``
    Protein-coding gene identifiers.  Used as an optional hard filter to
    exclude pseudogenes and ncRNAs.

``blacklist_prefixes.txt``
    Gene name prefixes to exclude (MT-, RPL, RPS, etc.).  Not
    genome-specific.

Notes
-----
All files are optional.  When a file is absent, the corresponding filter
returns an empty set with a ``logging.warning``.  The caller always
receives explicit information about what was (or was not) applied.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

__all__ = ["GeneFilterSet"]

logger = logging.getLogger("tissueresolve.reference.gene_filters")

# Default blacklist prefixes (used when the file is absent).
_DEFAULT_BLACKLIST_PREFIXES: list[str] = [
    "MT-", "RPL", "RPS", "MRPL", "MRPS",
    "SNHG", "LINC", "AC0", "AL0",
]


def _package_data_dir() -> Path:
    """Return path to ``src/tissueresolve/data/gene_lists/``."""
    return Path(__file__).parent.parent / "data" / "gene_lists"


class GeneFilterSet:
    """Instance-based gene filter set.

    Each instance has an independent cache so that two separate pipeline
    runs with different genomes cannot contaminate each other.

    Parameters
    ----------
    genome:
        Reference genome assembly (``"hg38"`` or ``"mm10"``).  Determines
        which genome-specific gene lists are loaded.
    data_dir:
        Override the path to the gene list directory.  Defaults to
        ``src/tissueresolve/data/gene_lists/``.
    """

    def __init__(
        self,
        genome: str = "hg38",
        data_dir: Optional[Path] = None,
    ) -> None:
        self.genome = genome
        self._data_dir: Path = Path(data_dir) if data_dir else _package_data_dir()
        self._cache: dict[str, frozenset] = {}
        self._prefix_cache: Optional[list[str]] = None

    # ------------------------------------------------------------------
    # Public accessors
    # ------------------------------------------------------------------

    def intronic_dominant(self) -> frozenset:
        """Genes with intronic retention bias (snRNA reference concern)."""
        return self._load_set(f"intronic_dominant_{self.genome}.tsv")

    def dissociation_stress(self) -> frozenset:
        """Dissociation-stress genes (scRNA reference concern)."""
        return self._load_set(f"dissociation_stress_{self.genome}.tsv")

    def length_biased_3prime(self) -> frozenset:
        """Long genes under-detected by 3′-UMI capture."""
        return self._load_set(f"length_biased_3prime_{self.genome}.tsv")

    def safe_universal(self) -> frozenset:
        """Genes that pass all three protocol filters."""
        return self._load_set(f"safe_universal_{self.genome}.tsv")

    def hypervariable_inflammatory(self) -> frozenset:
        """Hypervariable inflammatory genes (not genome-specific)."""
        return self._load_set("hypervariable_inflammatory.tsv")

    def protein_coding(self) -> frozenset:
        """Protein-coding gene identifiers."""
        return self._load_set(f"protein_coding_{self.genome}.tsv")

    def blacklist_prefixes(self) -> list[str]:
        """Gene-name prefixes to exclude (MT-, RPL, …)."""
        if self._prefix_cache is not None:
            return self._prefix_cache

        path = self._data_dir / "blacklist_prefixes.txt"
        if not path.exists():
            logger.warning(
                "blacklist_prefixes.txt not found at %s.  "
                "Using built-in defaults: %s.",
                path, _DEFAULT_BLACKLIST_PREFIXES,
            )
            self._prefix_cache = list(_DEFAULT_BLACKLIST_PREFIXES)
            return self._prefix_cache

        prefixes: list[str] = []
        with path.open(encoding="utf-8") as fh:
            for line in fh:
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    prefixes.append(stripped)

        if not prefixes:
            logger.warning(
                "blacklist_prefixes.txt is empty.  Using built-in defaults."
            )
            self._prefix_cache = list(_DEFAULT_BLACKLIST_PREFIXES)
        else:
            self._prefix_cache = prefixes
            logger.debug("Loaded %d blacklist prefixes.", len(prefixes))

        return self._prefix_cache

    # ------------------------------------------------------------------
    # Filter helpers
    # ------------------------------------------------------------------

    def filter_by_blacklist_prefixes(
        self,
        gene_names: list[str],
    ) -> tuple[list[str], list[str]]:
        """Partition *gene_names* into (kept, removed) by blacklist prefixes.

        Parameters
        ----------
        gene_names:
            Candidate gene identifiers.

        Returns
        -------
        tuple (kept, removed)
            Both are lists preserving the order of *gene_names*.
        """
        prefixes = self.blacklist_prefixes()
        kept, removed = [], []
        for g in gene_names:
            if any(g.startswith(p) for p in prefixes):
                removed.append(g)
            else:
                kept.append(g)
        if removed:
            logger.debug(
                "Blacklist prefix filter: removed %d / %d genes.",
                len(removed), len(gene_names),
            )
        return kept, removed

    def genes_removed_report(
        self, gene_names: list[str]
    ) -> dict[str, list[str]]:
        """Return a dict mapping filter name → genes removed by that filter.

        This provides an explicit record of what was removed so the caller
        can log or inspect it.  A gene may appear under multiple filters.

        Parameters
        ----------
        gene_names:
            Candidate gene list.

        Returns
        -------
        dict[str, list[str]]
            Keys: filter category names.
            Values: genes removed (may be empty).
        """
        gene_set = set(gene_names)
        return {
            "blacklist_prefix": [
                g for g in gene_names
                if any(g.startswith(p) for p in self.blacklist_prefixes())
            ],
            "not_protein_coding": (
                sorted(gene_set - self.protein_coding())
                if self.protein_coding()
                else []
            ),
            "hypervariable_inflammatory": sorted(
                gene_set & self.hypervariable_inflammatory()
            ),
            "intronic_dominant": sorted(gene_set & self.intronic_dominant()),
            "dissociation_stress": sorted(gene_set & self.dissociation_stress()),
            "length_biased_3prime": sorted(gene_set & self.length_biased_3prime()),
        }

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _load_set(self, filename: str) -> frozenset:
        """Load a gene set from a TSV file; return from cache if available."""
        if filename in self._cache:
            return self._cache[filename]

        path = self._data_dir / filename
        if not path.exists():
            logger.warning(
                "Gene list '%s' not found at %s.  "
                "Filter for this category is disabled for genome='%s'.",
                filename, path, self.genome,
            )
            self._cache[filename] = frozenset()
            return self._cache[filename]

        genes: set[str] = set()
        with path.open(encoding="utf-8") as fh:
            for line in fh:
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    genes.add(stripped)

        self._cache[filename] = frozenset(genes)
        logger.debug("Loaded %d genes from '%s'.", len(genes), filename)
        return self._cache[filename]

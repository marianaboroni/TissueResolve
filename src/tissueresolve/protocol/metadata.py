"""
Protocol metadata for TissueResolve.

Two conceptually distinct uses of "protocol" are handled by this module:

1. **Bulk + reference protocol compatibility** (bulk deconvolution)
   How a bulk RNA-seq library was prepared (polyA vs ribodepleted) and how the
   single-cell reference was generated (scRNA vs snRNA, 10x 3′ vs full-length,
   exonic-only vs exonic+intronic counts) determines which gene-level biases
   affect the deconvolution.  These are assessed by :mod:`protocol.risk`.

2. **Spatial platform metadata** (spatial deconvolution)
   The spatial technology used to generate a dataset (Visium, Visium HD,
   Slide-seq, etc.) and the modality of the reference used with it.
   Spatial-specific scale-factor correction is handled by
   :mod:`protocol.mismatch`.

All metadata classes are JSON-serialisable and round-trip through
:meth:`ProtocolMetadata.to_dict` / :meth:`ProtocolMetadata.from_dict`.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any

__all__ = [
    # Bulk enums
    "BulkProtocol",
    "RefModality",
    "RefCapture",
    "RefCounting",
    # Spatial enum
    "SpatialPlatform",
    # Metadata containers
    "ProtocolMetadata",
]


# ---------------------------------------------------------------------------
# Bulk-specific enumerations
# ---------------------------------------------------------------------------


class BulkProtocol(str, Enum):
    """Sequencing library preparation for bulk RNA-seq."""

    POLYA = "polyA"
    """polyA-selected library.  Long-gene 3′-bias present."""
    RIBODEP = "ribodepleted"
    """Ribosomal-RNA-depleted library.  Partial intronic signal expected."""
    UNKNOWN = "unknown"


class RefModality(str, Enum):
    """Single-cell / single-nucleus reference modality."""

    SCRNA = "scRNA"
    """Dissociated single-cell RNA-seq.  Dissociation-stress genes present."""
    SNRNA = "snRNA"
    """Single-nucleus RNA-seq.  Intronic reads retained; cytoplasmic RNA low."""
    UNKNOWN = "unknown"


class RefCapture(str, Enum):
    """Capture chemistry used for the single-cell reference."""

    UMI_3PRIME = "10x_3prime"
    """10x Genomics 3′-end UMI capture.  Under-detects genes > 5 000 bp."""
    UMI_5PRIME = "10x_5prime"
    """10x Genomics 5′-end UMI capture.  Less length bias than 3′."""
    FULL_LENGTH = "full_length"
    """Full-length capture (Smart-seq2, etc.).  No systematic length bias."""
    DROPSEQ = "dropseq"
    """Drop-seq.  3′-end bias similar to 10x 3′."""
    UNKNOWN = "unknown"


class RefCounting(str, Enum):
    """How reads were counted for the reference matrix."""

    EXONIC = "exonic"
    """Only exonic reads counted.  Compatible with both scRNA and snRNA refs."""
    EXONIC_INTRONIC = "exonic_intronic"
    """Exonic + intronic reads counted.  Amplifies intronic retention in snRNA."""
    UNKNOWN = "unknown"


# ---------------------------------------------------------------------------
# Spatial-specific enumeration (new in TissueResolve)
# ---------------------------------------------------------------------------


class SpatialPlatform(str, Enum):
    """Spatial transcriptomics platform.

    Used to tag spatial datasets and select platform-appropriate
    preprocessing and graph-building strategies.
    """

    VISIUM = "visium"
    """10x Genomics Visium (55 μm spot diameter, hexagonal array)."""
    VISIUM_HD = "visium_hd"
    """10x Genomics Visium HD (8 μm bins, square lattice)."""
    SLIDESEQ = "slideseq"
    """Slide-seq / Slide-seq v2 (10 μm puck)."""
    MERFISH = "merfish"
    """MERFISH / FISH-based in-situ transcriptomics."""
    XENIUM = "xenium"
    """10x Genomics Xenium in-situ platform."""
    UNKNOWN = "unknown"


# ---------------------------------------------------------------------------
# Unified metadata container
# ---------------------------------------------------------------------------


@dataclass
class ProtocolMetadata:
    """Describes the sequencing protocols of a bulk or spatial dataset and its
    single-cell reference.

    All fields default to ``UNKNOWN`` so that partial information is valid.
    Always provide as much information as available — more detail enables
    more targeted protocol-risk filtering.

    Attributes
    ----------
    bulk_protocol:
        Library preparation of the bulk RNA-seq data.
    ref_modality:
        Single-cell vs single-nucleus reference.
    ref_capture:
        Capture chemistry of the reference.
    ref_counting:
        Counting strategy of the reference (exonic vs exonic+intronic).
    spatial_platform:
        Spatial technology.  Only relevant for spatial workflows;
        leave as ``UNKNOWN`` for bulk deconvolution.
    """

    bulk_protocol: BulkProtocol = BulkProtocol.UNKNOWN
    ref_modality: RefModality = RefModality.UNKNOWN
    ref_capture: RefCapture = RefCapture.UNKNOWN
    ref_counting: RefCounting = RefCounting.UNKNOWN
    spatial_platform: SpatialPlatform = SpatialPlatform.UNKNOWN

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_strings(
        cls,
        *,
        bulk_protocol: str = "unknown",
        ref_modality: str = "unknown",
        ref_capture: str = "unknown",
        ref_counting: str = "unknown",
        spatial_platform: str = "unknown",
    ) -> "ProtocolMetadata":
        """Construct from plain strings (case-insensitive).

        Parameters
        ----------
        bulk_protocol:
            One of: ``polyA``, ``ribodepleted``, ``unknown``.
        ref_modality:
            One of: ``scRNA``, ``snRNA``, ``unknown``.
        ref_capture:
            One of: ``10x_3prime``, ``10x_5prime``, ``full_length``,
            ``dropseq``, ``unknown``.
        ref_counting:
            One of: ``exonic``, ``exonic_intronic``, ``unknown``.
        spatial_platform:
            One of: ``visium``, ``visium_hd``, ``slideseq``, ``merfish``,
            ``xenium``, ``unknown``.

        Raises
        ------
        ValueError
            When any string does not match a valid enum value.
        """
        return cls(
            bulk_protocol=_parse_enum(BulkProtocol, bulk_protocol),
            ref_modality=_parse_enum(RefModality, ref_modality),
            ref_capture=_parse_enum(RefCapture, ref_capture),
            ref_counting=_parse_enum(RefCounting, ref_counting),
            spatial_platform=_parse_enum(SpatialPlatform, spatial_platform),
        )

    # ------------------------------------------------------------------
    # Predicates
    # ------------------------------------------------------------------

    def is_fully_unknown(self) -> bool:
        """True when all fields are UNKNOWN."""
        return (
            self.bulk_protocol == BulkProtocol.UNKNOWN
            and self.ref_modality == RefModality.UNKNOWN
            and self.spatial_platform == SpatialPlatform.UNKNOWN
        )

    def is_bulk_fully_unknown(self) -> bool:
        """True when all bulk-relevant fields are UNKNOWN.

        If *any* of ``bulk_protocol``, ``ref_modality``, or ``ref_capture``
        is known, partial risk assessment is possible and this returns False.
        """
        return (
            self.bulk_protocol == BulkProtocol.UNKNOWN
            and self.ref_modality == RefModality.UNKNOWN
            and self.ref_capture == RefCapture.UNKNOWN
        )

    def is_spatial(self) -> bool:
        """True when a spatial platform has been specified."""
        return self.spatial_platform != SpatialPlatform.UNKNOWN

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def to_dict(self) -> dict[str, str]:
        """Return a plain-string dictionary (JSON-serialisable)."""
        return {
            "bulk_protocol": self.bulk_protocol.value,
            "ref_modality": self.ref_modality.value,
            "ref_capture": self.ref_capture.value,
            "ref_counting": self.ref_counting.value,
            "spatial_platform": self.spatial_platform.value,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "ProtocolMetadata":
        """Reconstruct from a plain-string dictionary.

        Tolerates extra keys and missing keys (uses ``"unknown"`` as default).
        """
        return cls.from_strings(
            bulk_protocol=d.get("bulk_protocol", "unknown"),
            ref_modality=d.get("ref_modality", "unknown"),
            ref_capture=d.get("ref_capture", "unknown"),
            ref_counting=d.get("ref_counting", "unknown"),
            spatial_platform=d.get("spatial_platform", "unknown"),
        )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _parse_enum(enum_cls: type, value: str) -> Enum:
    """Case-insensitive enum lookup.  Raises ``ValueError`` on no match."""
    value_lower = value.lower()
    for member in enum_cls:
        if member.value.lower() == value_lower:
            return member
    valid = [m.value for m in enum_cls]
    raise ValueError(
        f"'{value}' is not a valid {enum_cls.__name__}.  "
        f"Valid values: {valid}."
    )

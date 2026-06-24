# Pelka CRC Dataset — Verification (Stage 1)

Primary-source verification of the named second-tissue candidate **before** any
acquisition. Verified via primary repositories (GEO, HCA Data Explorer, Broad
Single Cell Portal), not mirrors. **Outcome: the named Pelka atlas FAILS the
acquisition gate in this environment; stopping for a decision (no silent
substitution).** Nothing downloaded or committed.

## Identity disambiguation (important)
Two distinct datasets share almost the same title:

- **Pelka et al. 2021, *Cell* 184:4734** — "Spatially organized multicellular
  immune hubs in human colorectal cancer." **371,223 cells, 41,364 genes, 62
  patients (28 MMRp + 34 MMRd)**, deep cTNI broad+fine labels, matched Visium in
  the paper. Hosted on **Broad Single Cell Portal SCP1162**; raw sequencing on
  **EGA study phs002407**. *(This is the atlas the plan recommended.)*
- **Masood/Indiana Univ. — GEO GSE200997 / HCA project 4d9d56e4** — "Refining
  Colorectal Cancer Classification … Single-Cell Atlas." **49,589 cells, 16 CRC +
  7 normal (~23 samples)**, open (CC BY 4.0), raw UMI count matrix + annotations +
  donor IDs. *A different, smaller study — must NOT be silently used as "Pelka".*

## Requirements table — Pelka et al. 2021 *Cell* atlas (the named candidate)
| Requirement | Required | Verified | Evidence | Decision |
|---|---|---|---|---|
| Raw / count matrices openly available | yes | **NO** | SCP1162 download = "Please sign in to download data"; raw FASTQs on EGA phs002407 (controlled) | **FAIL** |
| Anonymous / scriptable acquisition (no credentials) | yes | **NO** | SCP requires Broad/Terra sign-in; EGA requires a Data Access Committee application | **FAIL** |
| Multiple donors | yes | yes | 62 patients (28 MMRp + 34 MMRd) | pass |
| Broad labels | yes | likely | published cTNI broad compartments | pass (verify on access) |
| Fine labels | yes | likely | published fine cell programs / hubs | pass (verify on access) |
| Donor IDs | yes | likely | per-patient design | pass (verify on access) |
| Complex TME | yes | yes | epithelial + stromal + immune, MMRd/MMRp | pass |
| Reuse allowed (open license) | yes | **UNCONFIRMED/NO** | SCP terms not stated openly; EGA controlled | **FAIL** |
| Donor-disjoint evaluation possible | yes | yes (if acquired) | 62 donors | pass-if-acquired |

**Three critical requirements fail** (open raw counts, scriptable/no-credential
access, confirmed reuse license). Per the Stage-1 rule ("Stop if any critical
requirement fails"), the named Pelka atlas is **not acquirable in this
environment** without a credentialed/Data-Access-Committee process.

## Why I did not silently substitute
GSE200997 (Masood) is genuinely open and anonymously downloadable, **but it is a
different, smaller study** (~23 samples vs 62 donors; fine-label depth uncertain —
HCA lists some cell types as "Unspecified"). The rules forbid selecting a
replacement silently and forbid assuming annotations are suitable without audit.
The substitution is a **decision for the committee**, not an automatic swap.

## Options (for decision)
1. **Apply for Pelka access** (Broad SCP sign-in for the count matrices, and/or EGA
   phs002407 DAC for raw) — restores the *named* 62-donor atlas, but needs
   credentials I do not have and cannot place in code (rule). Out of scope for an
   automated session.
2. **Open CRC substitute — GSE200997 (Masood)**: anonymously downloadable raw UMI
   matrix (~68 MB) + annotations + donor IDs, CC BY 4.0. Smaller (~23 samples) and
   a different cohort; fine-annotation depth must be audited (Stage 3). Acceptable
   as an *open* CRC second tissue if the committee approves the substitution and
   accepts the reduced donor count / annotation audit.
3. **Other tissue** — NSCLC (Lambrechts E-MTAB-6149, ArrayExpress; check open raw
   counts) or HGSOC (Vázquez-García; often Synapse-gated). Each needs the same
   open-raw-count + scriptable-access verification.

## Recommendation
**Do not proceed to download/Task C.** The named Pelka atlas fails the open-access
gate. Recommend the committee choose: (a) authorise a credentialed Pelka/EGA
acquisition (outside this automated session), or (b) approve **GSE200997 as the
open CRC substitute** (smaller; annotation audit required), or (c) switch to a
verified-open NSCLC dataset. I will build the Stage-2 acquisition scaffold + run
Stages 3–N **only against a confirmed-open, anonymously-acquirable target** to
honour the no-credentials / no-silent-substitution / no-fabrication rules.

## Sources
- GEO GSE200997: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200997
- HCA project (GSE200997): https://explore.data.humancellatlas.org/projects/4d9d56e4-610d-4748-b57d-f8315e3f53a3
- Broad SCP1162 (Pelka 2021): https://singlecell.broadinstitute.org/single_cell/study/SCP1162
- EGA controlled study: https://www.ega-archive.org/studies/phs002407
- Pelka et al. 2021, *Cell* 184:4734–4752.

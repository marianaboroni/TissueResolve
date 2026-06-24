# Second-Tissue Dataset Candidates (Task B)

Prepared in parallel with Task A. **No data downloaded or fabricated** (rule 15);
all figures are best-known estimates to be **verified on acquisition** (record
exact values in a `download_manifest.json` then). Required properties: raw counts ·
multiple donors · broad+fine labels · adequate cells/subtype · epithelial/stromal/
immune diversity · donor-disjoint splits feasible · benchmark-compatible license.

## Candidates (priority order)

### 1. Colorectal cancer — Pelka et al. 2021, *Cell* (CRC atlas)  ★ top pick
- **Source/accession:** GEO GSE178341; also on CELLxGENE Census. *(verify)*
- **Tissue:** colorectal carcinoma + normal mucosa.
- **Donors:** ~60 patients *(verify)*. **Cells:** ~370k *(verify)*.
- **Counts:** raw UMI (10x 3'). **Labels:** broad + fine (cTNI immune/stromal/
  epithelial hierarchy) — well suited to a fine→broad mapping.
- **scRNA.** **License:** typically open (verify on CELLxGENE).
- **Strengths:** rich TME (epithelial+stromal+immune), many donors, deep fine
  labels, matched Visium exists in the literature → bulk *and* spatial validation.
- **Limitations:** large (subsample needed); fine labels are author-specific.
- **Decision:** **PREFERRED** — best match to the breast-atlas structure.

### 2. Ovarian cancer (HGSOC) — Vázquez-García et al. 2022, *Nature* (MSK SPECTRUM)
- **Source:** CELLxGENE / Synapse (syn25569736) *(verify)*.
- **Donors:** ~40+ patients. **Cells:** ~270k *(verify)*. **Counts:** raw UMI.
- **Labels:** broad + fine immune/stromal/epithelial. **scRNA.** License: verify.
- **Strengths:** multi-site, many donors, strong TME; distinct from breast.
- **Limitations:** treatment heterogeneity; access may require Synapse terms.
- **Decision:** **STRONG ALTERNATE**.

### 3. Lung cancer (NSCLC) — Lambrechts et al. 2018 / extended NSCLC atlases
- **Source:** CELLxGENE; ArrayExpress E-MTAB-6149 *(verify)*.
- **Donors:** ~8–35 depending on cohort. **Counts:** raw UMI. **Labels:** broad+fine.
- **Strengths:** independent tissue, immune-rich. **Limitations:** fewer donors in
  the original; consider a pooled NSCLC atlas. **Decision:** acceptable.

### 4. Adipose — Emont et al. 2022, *Nature* (human white adipose atlas)
- **Source:** CELLxGENE / Single Cell Portal *(verify)*.
- **Donors:** ~adipose from many subjects. **Counts:** raw (snRNA). **Labels:**
  broad+fine. **snRNA** (tests sc-vs-sn robustness — a useful stressor).
- **Strengths:** very different biology; snRNA cross-modality test. **Limitations:**
  fewer immune cells; snRNA mismatch with sc reference assumptions. **Decision:**
  good *secondary* stressor, not primary.

### 5. Kidney — KPMP / Lake et al. 2023, *Nature*
- **Source:** CELLxGENE / KPMP atlas *(verify)*. **Counts:** raw (sn/sc). **Labels:**
  detailed nephron + immune/stromal. **Strengths:** rich epithelial diversity.
  **Limitations:** snRNA-heavy; specialized epithelial hierarchy. **Decision:**
  alternate.

## Recommendation
Acquire **Pelka CRC (1)** first (closest structural analogue: TME diversity, many
donors, fine hierarchy, matched bulk+spatial), with **HGSOC (2)** as the alternate.
On acquisition: harmonise to gene symbols, define a fine→broad mapping, build
`download_manifest.json` (name, ID, version, date, donors, labels, protocol,
license, filters), and run the Stage-1 identifiability ceiling + the soft-gating
prospective gates (per `SECOND_TISSUE_SELECTION_AND_VALIDATION_PLAN.md`).

**Until a second tissue is acquired + validated, soft gating and any new signature
method remain experimental and unpromoted.** No further breast-only tuning while
selection is unresolved (per instruction).

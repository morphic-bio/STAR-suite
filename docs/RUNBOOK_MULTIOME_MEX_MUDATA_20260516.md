# Multiome MEX to MuData Runbook

Date: 2026-05-18

Status: implemented STAR-side multiome path. STAR-suite now has a reusable
MuData builder, a STAR/Chromap local MEX boundary, a remote post-MEX
CellBender/MuData handoff, and a manifest-driven JAX_Multiome01 production
wrapper. Earlier mentions of MuData in handoffs should be treated as
aspirational references that are superseded by the scripts listed here.

## Goal

Create a reproducible test workflow that turns the PBMC 3K multiome test set
into paired RNA and ATAC `AnnData` modalities inside one `.h5mu` object before
running the same shape on a production multiome set.

The RNA side should follow the UCSF/MSK h5ad pattern:

- raw-backed object plus filtered/default usable-cell views where possible
- `X` as count matrix, `layers["counts"]` when useful for downstream tools
- `.obs` carrying cell-calling/QC masks and provenance
- `.var` carrying feature ids, names, feature type, and genomic coordinates
- `.uns` carrying source paths and pipeline metadata

The ATAC side should follow the scverse/SnapATAC2 style:

- rows are cells/barcodes, columns are peaks
- `.X` is a sparse peak count matrix
- `.var_names` are stable peak ids such as `chr:start-end`
- `.var` includes `chrom`, `chromStart`, `chromEnd`, `feature_types = "Peaks"`
- `.obs` includes fragment/QC/evidence metrics such as total fragments, peak
  fragments, peak cut sites, peak fraction, and ATAC module call
- fragment and peak source paths are recorded in `.uns`

The final container should be:

```python
mudata.MuData({"rna": rna_adata, "atac": atac_adata})
```

with a shared observation axis (`axis=0`) and barcode names in one canonical
namespace.

## Primary References

Local STAR-suite:

- `scripts/build_multiome_mudata.py` - production-oriented RNA+ATAC `.h5mu`
  builder from RNA h5ad or GEX MEX plus ATAC peak MEX.
- `scripts/run_star_multiome_lane_smoke.sh` - STAR/Chromap multiome sample
  runner. Preferred production mode stops after local STAR cell calling,
  Velocyto packaging, Chromap BAM/sidecar materialization, and ATAC peak-MEX
  construction.
- `scripts/run_remote_multiome_post_mex_rsync.sh` - preferred post-MEX remote
  boundary. It stages RNA MEX plus local ATAC peak MEX/metrics to the GPU
  server, runs RNA downstream/CellBender and MuData there, then syncs only
  downstream outputs back.
- `scripts/run_jax_multiome01_production.sh` - manifest-driven production
  wrapper for the nine JAX_Multiome01 samples. It runs the local STAR/Chromap
  boundary sample-by-sample, queues serialized remote post-MEX jobs, and can
  invoke per-sample Globus large-file transfer/cleanup.
- `scripts/upload_jax_multiome01_large_files_globus.sh` - Globus uploader for
  raw FASTQs, generated Y/noY FASTQs, and BAMs. Generated local BAM/FASTQ files
  are deleted only after the Globus task succeeds; raw input FASTQs are
  preserved.
- `scripts/normalize_multiome_atac_barcode_fastq.py` - legacy fallback for
  materializing normalized ATAC barcode FASTQs. Production uses native
  `--chromapAtacReadFormat` instead.
- `scripts/run_multiome_mudata_smoke.sh` - smoke/production wrapper that splits
  ARC MEX, verifies or regenerates GeneFull and Velocyto MEX packaging, uses
  the remote downstream/CellBender path when configured, and writes MuData
  outputs.
- `scripts/package_star_genefull_mex.py` - packages STARsolo
  `Solo.out/GeneFull` raw/filtered MEX into CellRanger-style
  `outs/raw_feature_bc_matrix` and `outs/filtered_feature_bc_matrix`.
- `scripts/run_scrna_downstream_gene_full_velocyto.sh` - current downstream h5ad
  wrapper for UCSF/MSK-style RNA objects.
- `docs/RUNBOOK_VELOCYTO_CR_COMPAT_POLICY_20260519.md` - planned
  GeneFull-count-preserving Velocyto policy for production velocity layers where
  `spliced + unspliced + ambiguous` must equal the GeneFull count matrix.
- `scripts/run_remote_cellbender_rsync.sh` - remote-only CellBender runner for
  an existing downstream h5ad directory.
- `scripts/build_gene_full_velocyto_h5ad.py` - current GeneFull plus Velocyto
  h5ad builder pattern.
- `scripts/integrate_feature_library.py` - feature-library h5ad integration and
  barcode namespace handling.
- `scripts/extract_cr_feature_type_mex.py` - splits a Cell Ranger-style combined
  MEX by `features.tsv` feature type.
- `scripts/run_multiome_cell_call_from_arc.sh` - experimental ARC-derived
  GEX/ATAC evidence and multiome cell-calling driver.
- `scripts/run_multiome_cell_call_external_gex_from_arc.sh` - same path with
  GEX MEX extraction plus SimpleED.
- `core/features/libscrna/tools/scrna_*` - ARC table, GEX evidence, ATAC
  evidence, and combiner tools.
- `core/features/libchromap_contract/tools/star_multiome_atac_peak_mex.cpp` -
  native binary-sidecar peak caller, peak-MEX builder, and ATAC metrics writer.
- `core/features/libscrna/src/AtacEvidenceFromPeaks.cc` - lower-level
  binary-sidecar to per-barcode ATAC evidence implementation used by the native
  materialization path.
- `docs/LIBCHROMAP_CONTRACT.md` - STAR `--chromapAtac*` contract and validation.
- `docs/HANDOFF_STAR_LIBCHROMAP_MACS3_INTEGRATION_20260425.md` - prior
  handoff, including the now-addressed ATAC evidence-from-peaks direction.

Local Chromap-suite:

- `/mnt/pikachu/Chromap-suite/README.md` - PBMC 3K multiome benchmark and
  sidecar/peak output summary.
- `/mnt/pikachu/Chromap-suite/mcp_server/config.yaml` - dataset id
  `pbmc_unsorted_3k_multiome`.
- `/mnt/pikachu/Chromap-suite/mcp_server/workflows/chromap_atac_bam_fragments.yaml`
  - standalone ATAC BAM/fragments/peaks workflow schema.
- `/mnt/pikachu/Chromap-suite/mcp_server/recipes/registry.yaml` - handoff
  artifacts for fragments, binary sidecar, peaks, and ATAC evidence.

External model docs:

- MuData docs: modalities are stored as `AnnData` objects under `.mod`, with
  shared observations represented through obs maps in `.h5mu`.
  https://mudata.readthedocs.io/stable/io/spec.html
- MuData quickstart: `MuData({"A": adata, "B": adata2})`, `mdata.update()`,
  and `mdata.write(...)` are the expected construction/writer pattern.
  https://mudata.readthedocs.io/stable/notebooks/quickstart_mudata.html
- SnapATAC2 multiome tutorial: ATAC peak matrices are represented as AnnData
  with peak features and `chrom`, `chromStart`, `chromEnd` annotations before
  RNA/ATAC integration.
  https://scverse.org/SnapATAC2/version/2.2/tutorials/modality.html
- SnapATAC2 import/peak matrix APIs: `pp.import_data` imports fragments and
  computes basic ATAC QC; `pp.make_peak_matrix` builds a cell-by-peak matrix
  from peaks.
  https://scverse.org/SnapATAC2/version/2.5/api/_autosummary/snapatac2.pp.import_data.html
  https://scverse.org/SnapATAC2/version/2.5/api/_autosummary/snapatac2.pp.make_peak_matrix.html

## Test Datasets

Use the 100K fixture first:

```bash
TEST_ROOT=/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k
ARC_OUT=${TEST_ROOT}/pbmc_unsorted_3k_100k_arc/outs
CHROMAP_INDEX=${TEST_ROOT}/chromap_index
OUT=tests/multiome_mudata_smoke_output_$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "${OUT}"
```

Known useful fixture files:

- ARC combined MEX: `${ARC_OUT}/raw_feature_bc_matrix` and
  `${ARC_OUT}/filtered_feature_bc_matrix`
- ARC barcode metrics: `${ARC_OUT}/per_barcode_metrics.csv`
- ARC ATAC whitelist/translation support:
  `${CHROMAP_INDEX}/737K-arc-v1_atac.txt` and
  `${CHROMAP_INDEX}/atac2gex.tsv`
- STAR/Chromap 100K outputs exist under the same fixture root and can be used
  after the ARC-only packaging prototype is stable.

The full public PBMC 3K multiome set remains useful for regression and parity
checks:

```bash
FULL_ROOT=/mnt/pikachu/atac-seq/10xMultiome/pbmc_unsorted_3k
FULL_ARC_OUT=${FULL_ROOT}/pbmc_unsorted_3k_full_arc_live/outs
```

For new schema work, do not start with production data. The 100K fixture has
both ARC MEX and ARC per-barcode metrics, so it remains the safest schema and
parity target. For the current production implementation, the JAX one-lane smoke
is the closest end-to-end rehearsal because it exercises native ATAC barcode
parsing, low-memory Chromap, remote post-MEX, and MuData assembly.

## GEX Reference Policy

Use the same STAR GEX reference used by the UCSF and MSK perturb-seq runs:

```bash
GEX_GENOME_DIR=/storage/autoindex_110_44/bulk_index
```

This is the shared UCSF/MSK STAR index and is the required `--genomeDir` for
the expression side of the multiome smoke and JAX KOLF-2 production runs. Do
not substitute the Cell Ranger ARC multiome reference or a Chromap ATAC index
for STAR GEX alignment.

The local index records a Cell Ranger-style GRCh38 reference surface:

- `--cellrangerRefRelease 2024-A`
- `sjdbGTFfile=/storage/autoindex_110_44/bulk_index/cellranger_ref/genes.gtf`
- `genomeFastaFiles=/storage/autoindex_110_44/bulk_index/cellranger_ref/genome.fa`

Before running Phase 2 or Phase 3, validate:

```bash
test -d "${GEX_GENOME_DIR}"
grep -E 'cellrangerRefRelease|sjdbGTFfile|genomeFastaFiles' \
  "${GEX_GENOME_DIR}/genomeParameters.txt"
```

ATAC chrom sizes and peak coordinates should remain GRCh38-compatible with this
GEX reference. For SnapATAC2 prototypes, `snap.genome.hg38` is acceptable only
after confirming contig names match the Chromap fragments and peak file.

Current STAR/Chromap production defaults:

```bash
GEX_GENOME_DIR=/storage/autoindex_110_44/bulk_index
GEX_WHITELIST=/mnt/pikachu/atac-seq/10xMultiome/pbmc_unsorted_3k/open_source_full_20260424_015259/refs/737K-arc-v1_gex.txt
CHROMAP_REF=/storage/autoindex_110_44/bulk_index/cellranger_ref/genome.fa
CHROMAP_INDEX=/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/genome.index
ATAC_WHITELIST=/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt
ATAC_TO_GEX=/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/atac2gex.tsv
```

The JAX ARC ATAC barcode read uses bases 9-24 of the ATAC i5/read2 file,
reverse-complemented. Production should express this natively as:

```text
--chromapAtacReadFormat bc:8:23:-
--chromapAtacBarcodeTranslate <ATAC_TO_GEX>
--chromapAtacBarcodeTranslateFromFirst 1
```

Do not pre-normalize the ATAC barcode FASTQ unless debugging the legacy fallback
path explicitly.

## Implemented Smoke Driver

The primary entry point is:

```bash
scripts/run_multiome_mudata_smoke.sh \
  --arc-outs "${ARC_OUT}" \
  --out-dir "${OUT}"
```

This ARC-only mode writes split RNA and ATAC peak MEX directories and then
attempts ARC-only `.h5mu` outputs. The Python environment must include
`mudata`; without it, the split MEX step can still be tested with:

```bash
scripts/run_multiome_mudata_smoke.sh \
  --arc-outs "${ARC_OUT}" \
  --out-dir "${OUT}" \
  --skip-arc-only \
  --skip-star-hybrid
```

For STAR RNA plus ARC ATAC, use a STAR run that was produced with the GEX
reference and Y-removal policy above:

```bash
scripts/run_multiome_mudata_smoke.sh \
  --arc-outs "${ARC_OUT}" \
  --star-run-dir /path/to/sample/run \
  --remote-host 10.159.4.53 \
  --remote-root /path/on/remote/local_disk/multiome_cellbender \
  --cellbender-gpu \
  --cellbender-cpu-cores 24 \
  --out-dir "${OUT}" \
  --force
```

The tested 100K path used an existing local downstream h5ad directory, then
offloaded only CellBender to the remote GPU host:

```bash
scripts/run_scrna_downstream_gene_full_velocyto.sh \
  --run-dir <sample>/run \
  --output-dir <sample>/downstream_genefull_velocyto_smoke_qc \
  --min-genes 1 \
  --max-genes 100000 \
  --mt-pct-cutoff 100 \
  --adaptive-filter

scripts/run_multiome_mudata_smoke.sh \
  --arc-outs "${ARC_OUT}" \
  --star-run-dir <sample>/run \
  --sample-dir <sample> \
  --downstream-dir <sample>/downstream_genefull_velocyto_smoke_qc \
  --remote-host 10.159.4.53 \
  --remote-root /home/lhhung/multiome_remote_downstream_smoke \
  --cellbender-gpu \
  --cellbender-cpu-cores 8 \
  --skip-arc-only \
  --out-dir "${OUT}" \
  --force
```

The relaxed downstream thresholds above are for the tiny 100K fixture only; the
production path should keep the UCSF/MSK-style filtering policy unless a
dataset-specific QC decision says otherwise.

Remote host/root can also be supplied through:

```bash
export MULTIOME_REMOTE_HOST=10.159.4.53
export MULTIOME_REMOTE_ROOT=/path/on/remote/local_disk/multiome_cellbender
```

When `--star-run-dir` is provided, the wrapper checks
`outs/raw_feature_bc_matrix`, `outs/filtered_feature_bc_matrix`,
`outs/raw_velocyto_feature_bc_matrix`,
`outs/filtered_velocyto_feature_bc_matrix`, and
the GeneFull/Velocyto manifests. If any required GeneFull MEX file is missing,
it automatically runs:

```bash
python3 scripts/package_star_genefull_mex.py --run-dir <star_run_dir>
```

Velocyto MEX is produced by STAR's internal post-Solo writer. If the required
Velocyto MEX files are missing, treat the STAR run as incomplete or rerun with a
current binary. The legacy `prepare_velocyto_mex.py` path is a repair/backfill
tool only and is not part of the preferred production boundary.

When STAR does not write `Solo.out/Velocyto/filtered/barcodes.tsv`, the
Velocyto packager uses `Solo.out/GeneFull/filtered/barcodes.tsv` as the
filtered barcode axis and subsets the raw Velocyto layers accordingly.

Preferred post-MEX boundary:

- local STAR/Chromap produces the RNA GeneFull MEX, Velocyto MEX, ATAC BAM,
  binary ATAC fragment sidecar, sidecar-derived peaks/summits, ATAC peak MEX,
  and ATAC metrics;
- the remote GPU server receives only post-MEX inputs, runs RNA downstream
  h5ad/CellBender and MuData assembly, then rsyncs downstream h5ad and `.h5mu`
  outputs back;
- STAR, Chromap, libchromap, libMACS3, and the native sidecar reader stay on
  the local STAR/Chromap host. This is the preferred boundary for production
  now, but the downstream side may evolve as the MuData schema matures.

Remote CellBender behavior for older ARC/hybrid wrappers:

- if no downstream dir is supplied, the wrapper uses
  `scripts/run_remote_scrna_downstream_rsync.sh` with `--run-cellbender`;
- if an existing downstream dir is supplied but lacks
  `cellbender/cellbender_counts.h5`, it uses
  `scripts/run_remote_cellbender_rsync.sh`;
- local downstream without remote CellBender is refused unless
  `--allow-local-downstream` is passed explicitly.

## Current STAR/Chromap Driver

Use `scripts/run_star_multiome_lane_smoke.sh` for one sample. FASTQ arguments
may be comma-separated lane lists; this is how the production wrapper passes
the multi-lane JAX libraries.

For production-shaped local materialization only:

```bash
scripts/run_star_multiome_lane_smoke.sh \
  --gex-r1 <gex_R1.fastq.gz[,lane2...]> \
  --gex-r2 <gex_R2.fastq.gz[,lane2...]> \
  --atac-r1 <atac_R1.fastq.gz[,lane2...]> \
  --atac-barcode <atac_R2_i5.fastq.gz[,lane2...]> \
  --atac-r2 <atac_R3.fastq.gz[,lane2...]> \
  --out-dir <sample_out> \
  --threads 16 \
  --chromap-threads 16 \
  --chromap-low-mem \
  --chromap-macs3-frag-low-mem \
  --chromap-start-mode concurrent \
  --stop-after-local-mex
```

This writes `LOCAL_MEX_READY.txt` only after all local boundary artifacts are
present:

- STAR GeneFull raw/filtered MEX and Velocyto raw/filtered MEX under
  `star_sample/run/outs/`;
- Y/noY BAM and FASTQ outputs from the expression run;
- Chromap ATAC sorted BAM, binary fragment sidecar, and `.chroms.tsv`;
- sidecar-derived `atac_peaks.narrowPeak` and `atac_summits.bed`;
- native C++ ATAC peak MEX under `atac/peak_mex/`;
- per-barcode ATAC metrics in `atac/atac_metrics.tsv`.

Post-MEX remote handoff:

```bash
scripts/run_remote_multiome_post_mex_rsync.sh \
  --sample-dir <sample_out>/star_sample \
  --remote-host 10.159.4.53 \
  --remote-root /home/lhhung/jax_multiome_remote_downstream_production \
  --output-name downstream_genefull_velocyto_cellbender \
  --run-cellbender \
  --adaptive-filter \
  --cellbender-cpu-cores 24 \
  --cellbender-gpu
```

The remote host receives only RNA MEX, Velocyto MEX, ATAC peak MEX, ATAC
metrics, and the scripts needed for downstream processing, including the
adaptive n_genes/MT QC helpers. It runs
`run_scrna_downstream_gene_full_velocyto.sh`, CellBender, and
`build_multiome_mudata.py`, validates required RNA Velocyto layers plus ATAC
`counts`, syncs h5ad and `.h5mu` outputs back, excludes the remote
`mudata_venv/`, and writes `REMOTE_POST_MEX_READY.txt`.

For sparse smoke tests only, the filtered RNA/ATAC barcode intersection may be
empty. The remote wrapper passes `--allow-empty-barcode-intersection` for the
filtered `.h5mu` so the path can be validated; production samples should treat
zero filtered observations as a QC warning that requires review.

## Barcode Namespace Policy

Use the GEX/combined Cell Ranger barcode namespace as the MuData `.obs_names`.

Rules:

- strip or preserve `-1` consistently across both modalities. ARC-only outputs
  preserve ARC's `-1` suffixes. STAR-hybrid outputs strip the suffix while
  building MuData because STARsolo GeneFull barcodes are unsuffixed and ARC MEX
  barcodes are suffixed.
- if Chromap fragments are in ATAC-barcode namespace, translate to GEX namespace
  using `atac2gex.tsv` before building the final ATAC modality, or run Chromap
  with `chromapAtacBarcodeTranslate` so downstream evidence already uses the
  GEX namespace.
- fail on duplicate barcodes after normalization.
- write both `barcode_raw` and `barcode_canonical` columns in each modality's
  `.obs` for traceability.

## RNA Production Path Parity

The ARC-only Phase 1 below does not run STAR, so it cannot exercise Y-removal
or STAR Velocyto output. It is only a schema smoke test.

For any Phase 2 or Phase 3 run where STAR produces the expression side, use the
same downstream h5ad surface as the UCSF/MSK perturb-seq work:

- run STAR with `--soloFeatures GeneFull Velocyto`.
- package GeneFull MEX with `scripts/package_star_genefull_mex.py`, producing:
  - `outs/raw_feature_bc_matrix/`
  - `outs/filtered_feature_bc_matrix/`
  - `outs/gene_full_feature_bc_matrix_manifest.json`
- STAR internally packages Velocyto MEX, producing:
  - `outs/raw_velocyto_feature_bc_matrix/`
  - `outs/filtered_velocyto_feature_bc_matrix/`
  - `outs/velocyto_feature_bc_matrix_manifest.json`
- build the RNA AnnData through `scripts/run_scrna_downstream_gene_full_velocyto.sh`,
  not by hand-reading `Solo.out/GeneFull` directly.
- carry the RNA object into MuData from the downstream h5ad output, preserving
  `layers["spliced"]`, `layers["unspliced"]`, and `layers["ambiguous"]` when
  present.

For Y-removal, follow the UCSF path rather than the MSK path. MSK does not use
Y-removal, but the intended production set is the JAX KOLF-2 line and does
require Y-removal. The STAR expression run should therefore enable integrated
Y/noY outputs:

```text
--outSAMtype BAM Unsorted
--emitNoYBAM yes
--emitYNoYFastq yes
--emitYNoYFastqCompression gz
```

For this multiome smoke run and the JAX KOLF-2 production run, keep these
Y-removal flags enabled. Record `y_removal_enabled = true` in both the run
manifest and `mdata.uns`.

## Phase 1: ARC MEX Prototype

This phase proves the MuData writer and schemas without relying on a new STAR
run.

Split ARC combined MEX into GEX and peak-only matrices:

```bash
python3 scripts/extract_cr_feature_type_mex.py \
  --input-mex-dir "${ARC_OUT}/raw_feature_bc_matrix" \
  --feature-type "Gene Expression" \
  --out-dir "${OUT}/raw_gex_mex"

python3 scripts/extract_cr_feature_type_mex.py \
  --input-mex-dir "${ARC_OUT}/filtered_feature_bc_matrix" \
  --feature-type "Gene Expression" \
  --out-dir "${OUT}/filtered_gex_mex"

python3 scripts/extract_cr_feature_type_mex.py \
  --input-mex-dir "${ARC_OUT}/raw_feature_bc_matrix" \
  --feature-type "Peaks" \
  --out-dir "${OUT}/raw_atac_peak_mex"

python3 scripts/extract_cr_feature_type_mex.py \
  --input-mex-dir "${ARC_OUT}/filtered_feature_bc_matrix" \
  --feature-type "Peaks" \
  --out-dir "${OUT}/filtered_atac_peak_mex"
```

Build a prototype `pbmc_100k_arc_unfiltered.h5mu` from `raw_*` and
`pbmc_100k_arc_filtered.h5mu` from `filtered_*`.

Implementation notes for the first writer script:

- read each MEX with `scanpy.read_10x_mtx(..., var_names="gene_ids")`.
- preserve sparse CSR matrices.
- copy raw counts to `layers["counts"]` for both RNA and ATAC.
- RNA `.var` should preserve Ensembl IDs as index, `gene_symbols`, and
  `feature_types`.
- ATAC `.var` should parse peak ids into `chrom`, `chromStart`, and `chromEnd`.
- join `per_barcode_metrics.csv` by the canonical barcode and add at least:
  `gex_umis_count`, `gex_genes_count`, `atac_fragments`,
  `atac_peak_region_fragments`, `atac_peak_region_cutsites`,
  `atac_peak_fraction`, `arc_is_cell`, `arc_excluded_reason`.
- set `.obs["is_cell"]` from the filtered barcode set for the filtered object
  and from ARC `is_cell` for ARC reference views.
- construct `MuData({"rna": rna, "atac": atac})`, run `mdata.update()`, then
  write `.h5mu`.

The implemented writer is `scripts/build_multiome_mudata.py`. Keep
notebook-only exploration out of production.

## Phase 2: STAR GEX + ARC ATAC Hybrid

This phase checks whether the STAR expression side can be swapped in while the
ATAC side remains ARC-derived and therefore easy to validate.

Inputs:

- STAR run directory with CellRanger-compatible `outs/` packaging:
  - `outs/raw_feature_bc_matrix`
  - `outs/filtered_feature_bc_matrix`
  - `outs/gene_full_feature_bc_matrix_manifest.json`
  - `outs/raw_velocyto_feature_bc_matrix`
  - `outs/filtered_velocyto_feature_bc_matrix`
  - `outs/velocyto_feature_bc_matrix_manifest.json`
- downstream RNA h5ad outputs from `scripts/run_scrna_downstream_gene_full_velocyto.sh`
- ARC peak-only MEX from Phase 1
- ARC `per_barcode_metrics.csv` for reference labels and ATAC metrics

Expected work:

1. Run or reuse a STAR expression run with `--soloFeatures GeneFull Velocyto`
   and the Y/noY flags from "RNA Production Path Parity". Pass
   `--genomeDir "${GEX_GENOME_DIR}"` using the reference from "GEX Reference
   Policy".
2. Run `scripts/package_star_genefull_mex.py --run-dir <star_run_dir>` if the
   packaged GeneFull MEX surface is not already present.
3. Confirm STAR wrote `outs/raw_velocyto_feature_bc_matrix`,
   `outs/filtered_velocyto_feature_bc_matrix`, and
   `outs/velocyto_feature_bc_matrix_manifest.json`.
4. Run `scripts/run_scrna_downstream_gene_full_velocyto.sh --run-dir <star_run_dir>
   --output-dir <downstream_dir>` using the same options selected for the
   UCSF/MSK-style h5ad path, including adaptive n_genes/MT QC and CellBender
   when the test goal includes it. Y-removal remains enabled independently of
   this downstream h5ad wrapper.
5. Build `rna` for MuData from the downstream h5ad, preserving counts and
   Velocyto layers instead of re-reading only GeneFull MEX.
6. Build `atac` from ARC peak MEX.
7. Inner-join observations for the strict filtered `.h5mu`; keep a raw/unfiltered
   `.h5mu` if barcode sets differ.
8. Validate RNA count sums, filtered barcode counts, and Velocyto layer shapes
   against STAR logs, packaged MEX, and the downstream h5ad summary.
9. Validate Y-removal artifacts when enabled:
   - `Aligned.out_Y.bam`
   - `Aligned.out_noY.bam`
   - `y_separated/*.fastq.gz`
10. Validate ATAC feature/cell counts against ARC MEX.

This isolates expression-side packaging from Chromap/SnapATAC2 peak-matrix
generation.

## Phase 3: STAR/Chromap ATAC Path

This phase replaces ARC ATAC with Chromap-suite outputs.

Preferred production-shaped ATAC artifacts from STAR:

- `atac_possorted.bam`
- `atac_fragments.bin` plus `atac_fragments.bin.chroms.tsv` binary sidecar

Preferred local post-materialization ATAC artifacts:

- `atac_peaks.narrowPeak`
- `atac_summits.bed`
- `atac/peak_mex/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}`
- `atac/atac_metrics.tsv`

STAR must be built cleanly with Chromap support before debugging this path:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR WITH_CHROMAP=1
```

STAR parameter shape:

```text
--chromapAtacEnable 1
--chromapAtacStartMode concurrent
--chromapAtacOutputFormat BAM
--chromapAtacOutputFragments <out>/atac_possorted.bam
--chromapAtacSecondaryFragments <out>/atac_fragments.bin
--chromapAtacSortBam 1
--chromapAtacReadFormat bc:8:23:-
--chromapAtacBarcodeTranslate <refs>/atac2gex.tsv
--chromapAtacBarcodeTranslateFromFirst 1
--chromapAtacLowMem 1
--chromapAtacLowMemRam 0
--chromapAtacMacs3FragLowMem 1
--chromapAtacTempDir <sample>/star_sample/chromap_tmp
```

After STAR/Chromap materializes the sidecar, build peaks, ATAC peak MEX, and
metrics locally:

```bash
core/features/libchromap_contract/star_multiome_atac_peak_mex \
  --sidecar <out>/atac_fragments.bin \
  --barcode-translate <refs>/atac2gex.tsv \
  --barcode-translate-from-first \
  --call-peaks-from-sidecar \
  --peaks <out>/atac_peaks.narrowPeak \
  --summits-out <out>/atac_summits.bed \
  --out-dir <sample>/atac/peak_mex \
  --metrics-tsv <sample>/atac/atac_metrics.tsv \
  --threads 16 \
  --temp-dir <sample>/star_sample/chromap_tmp
```

Do not use the removed `chromapAtacMacs3FragPeaksSource` STAR flag. The old
file-source path could spill/re-read a misleading `.tsv.gz` fragments file and
is not the production boundary. Low-memory production runs should spill through
Chromap's native low-memory path and the binary sidecar, then read that sidecar
once in `star_multiome_atac_peak_mex`.

The current sample wrapper uses fixed thread budgets for STAR and Chromap
(`--threads 16` and `--chromap-threads 16` in production) rather than dynamic
permit sharing. Keep those budgets explicit when using concurrent mode.

ATAC AnnData construction options:

1. SnapATAC2 route, only when an explicit text fragment export exists:

   ```python
   import snapatac2 as snap

   adata_frag = snap.pp.import_data(
       "<exported_atac_fragments.tsv.gz>",
       chrom_sizes=snap.genome.hg38,
       file="atac_fragments.h5ad",
       sorted_by_barcode=False,
       whitelist="final_barcodes.txt",
       shift_left=0,
       shift_right=0,
   )
   atac = snap.pp.make_peak_matrix(
       adata_frag,
       peak_file="atac_peaks.narrowPeak",
       file="atac_peaks.h5ad",
   )
   ```

2. Direct MEX route, recommended for the current STAR/Chromap production path:
   read with `scanpy.read_10x_mtx`, then apply the same ATAC `.var` and `.obs`
   schema as Phase 1.

SnapATAC2 import expects fragment rows with at least:

```text
chrom  start  end  barcode  count
```

Chromap-suite text fragments match this downstream shape when explicitly
exported. The current production wrapper does not export a text fragments file
for ATAC AnnData construction; it reads the native ATAC peak MEX produced from
the binary sidecar.

## ATAC Evidence and Cell Calls

There are two evidence surfaces:

1. ARC-derived evidence, built by:

   ```bash
   scripts/run_multiome_cell_call_from_arc.sh
   ```

   This produces `gex_evidence.tsv`, `atac_evidence.tsv`,
   `multiome_calls.tsv`, and comparison outputs against ARC.

2. STAR/libchromap sidecar-derived evidence, written by the native ATAC
   peak-MEX materializer as `atac/atac_metrics.tsv`.

Current STAR/libchromap ATAC metrics include:

```text
barcode
atac_peak_region_cutsites
atac_peak_region_fragments
atac_fragments
atac_peak_fraction
```

Before feeding this into `scrna_multiome_combine`, add a small normalization
step or extend the toolchain to derive:

```text
atac_module_call = atac_peak_region_cutsites >= 1
atac_low_targeting = atac_peak_fraction < <threshold> when threshold > 0
atac_source = chromap_sidecar_peak_mex
```

Until that adapter exists, keep ARC-derived multiome calls as the validation
baseline and attach STAR/libchromap evidence to `.obs` as metrics only. The
optional STAR flag `--chromapAtacEvidenceFromPeaksOutput` belongs to the older
in-STAR peak/evidence path and is not part of the current production wrapper.

## Target MuData Schema

Write two outputs for the smoke:

- `unfiltered_multiome.h5mu`: all barcodes present in both modalities after
  namespace normalization.
- `filtered_multiome.h5mu`: final usable cells only. For Phase 1 this can be
  ARC `is_cell`; for the current STAR/Chromap production path it is the
  intersection of the downstream RNA filtered h5ad and ATAC peak MEX barcodes,
  with `cell_call_source = star_downstream_filtered_h5ad_chromap_atac`. Once
  the multiome combiner is promoted, this should move to
  `multiome_calls.final_is_cell`.

Global `.obs` columns:

- `barcode_raw`
- `barcode_canonical`
- `is_cell`
- `cell_call_source`
- `call_provenance`
- `gex_module_call`
- `atac_module_call`
- `effective_atac_module_call`
- `gex_rescue_eligible`
- `atac_low_targeting`
- `arc_is_cell` when ARC labels are available

`rna` modality:

- `X`: GeneFull counts
- `layers["counts"]`: raw count copy
- `layers["spliced"]`, `layers["unspliced"]`, `layers["ambiguous"]` when the
  expression side came from a STAR Velocyto run
- `.obs`: RNA QC and cell-call columns (`gex_umis_count`, `gex_genes_count`,
  `is_cell`)
- `.var`: `gene_symbols`, `feature_types`, and optional genomic coordinates
- `.uns`: `gene_expression_source`, `gene_expression_feature_kind`, STAR/ARC
  source metadata, Velocyto source metadata, and Y-removal mode

`atac` modality:

- `X`: peak counts
- `layers["counts"]`: raw peak count copy
- `.obs`: `atac_fragments`, `atac_peak_region_fragments`,
  `atac_peak_region_cutsites`, `atac_peak_fraction`, ATAC module call fields
- `.var`: `feature_types`, `chrom`, `chromStart`, `chromEnd`
- `.uns`: `fragments_source`, `peaks_source`, `evidence_source`,
  `peak_caller`, `tn5_shift_mode`

## Validation Gates

Run these before trying the full PBMC 3K set:

1. MEX split checks:
   - raw ARC 100K features split into 36,601 Gene Expression features and 2,847
     Peaks in the inspected fixture.
   - filtered ARC 100K barcodes: 2,002.
   - raw ARC 100K barcodes: 36,354.

2. MuData structural checks:

   ```python
   import mudata as md

   m = md.read("unfiltered_multiome.h5mu", backed=True)
   assert set(m.mod.keys()) == {"rna", "atac"}
   assert m.mod["rna"].n_obs == m.mod["atac"].n_obs
   assert "chrom" in m.mod["atac"].var
   assert "gene_symbols" in m.mod["rna"].var
   m.file.close()
   ```

3. Count checks:
   - RNA matrix sum equals source Gene Expression MEX sum.
   - ATAC matrix sum equals source peak MEX sum or SnapATAC2 peak matrix sum.
   - no duplicate `.obs_names` in either modality.
   - Phase 2/3 RNA objects include `spliced`, `unspliced`, and `ambiguous`
     layers with the same shape as `rna.X`.

4. Cell set checks:
   - Phase 1 filtered cells match ARC filtered MEX barcode set exactly.
   - Phase 2 STAR RNA filtered cells match STAR filtered barcode set exactly.
   - Current Phase 3 filtered cells match the barcode intersection between the
     downstream RNA filtered h5ad and the ATAC peak MEX.
   - Future promoted multiome-combiner runs should instead match
     `multiome_calls.tsv` where `final_is_cell == 1`.

5. Y-removal checks when enabled:
   - `Aligned.out_Y.bam` and `Aligned.out_noY.bam` exist.
   - `y_separated/` contains paired Y/noY FASTQs.
   - the manifest records `y_removal_enabled = true`.

6. Evidence checks:
   - `atac_peak_fraction == atac_peak_region_fragments / atac_fragments` when
     `atac_fragments > 0`.
   - `atac_module_call` threshold is recorded in `.uns`.
   - ARC comparison summaries are stored under the smoke output root.

7. Reload checks:
   - `mudata.read_h5mu(..., backed=True)` succeeds.
   - each modality can be read independently from the `.h5mu`.
   - filtered `.h5mu` with zero observations is acceptable only for sparse smoke
     tests that passed `--allow-empty-barcode-intersection`; production output
     should be reviewed before handoff if this happens.

## Tested STAR-Side Smoke

Completed on 2026-05-16 in:

```bash
tests/multiome_mudata_smoke_output_star_end2end_20260516T045156Z
```

Key artifacts:

- exact STAR command:
  `tests/multiome_mudata_smoke_output_star_end2end_20260516T045156Z/RUN_STAR_GEX.sh`
- STAR run:
  `tests/multiome_mudata_smoke_output_star_end2end_20260516T045156Z/star_sample/run`
- downstream h5ad with remote CellBender:
  `tests/multiome_mudata_smoke_output_star_end2end_20260516T045156Z/star_sample/downstream_genefull_velocyto_smoke_qc`
- MuData outputs:
  - `star_arc_unfiltered_multiome.h5mu`
  - `star_arc_filtered_multiome.h5mu`

STAR-side checks:

- clean rebuilt `core/legacy/source/STAR` before running.
- `--genomeDir /storage/autoindex_110_44/bulk_index`.
- `--soloFeatures GeneFull Velocyto`.
- Y-removal enabled with `--emitNoYBAM yes`, `--emitYNoYFastq yes`, and
  `--emitYNoYFastqCompression gz`.
- `Aligned.out_Y.bam` contained 242 reads.
- `Aligned.out_noY.bam` had 0 `chrY` alignments.
- `y_separated/` contained paired Y/noY FASTQs.

RNA/downstream checks:

- `filtered_counts.h5ad`: 2,693 cells x 38,606 genes.
- `final_counts.h5ad`: 736,320 barcodes x 38,606 genes.
- RNA layers present after remote CellBender:
  `spliced`, `unspliced`, `ambiguous`, and `denoised`.
- The PBMC 100K fixture produced zero nonzero Velocyto counts in current STAR
  outputs, matching older local UCSF 100K smoke artifacts. The layer files and
  shapes are still validated. For production, treat unexpectedly zero Velocyto
  counts as a dataset/reference issue to investigate.

Remote CellBender checks:

- ran on `10.159.4.53` with `biodepot/cellbender:0.3.2`.
- wrote `cellbender/cellbender_counts.h5`.
- no `CELLBENDER_FAILED.txt` remained after the successful rerun.
- propagated `layers["denoised"]` into `unfiltered_counts.h5ad`,
  `filtered_counts.h5ad`, and `final_counts.h5ad`.

MuData checks:

- `star_arc_unfiltered_multiome.h5mu`: 36,354 obs, 38,606 RNA vars, 2,847 ATAC vars.
- `star_arc_filtered_multiome.h5mu`: 2,693 obs, 38,606 RNA vars, 2,847 ATAC vars.
- modalities: `rna`, `atac`.
- RNA layers: `ambiguous`, `counts`, `denoised`, `spliced`, `unspliced`.
- ATAC layers: `counts`.
- global and modality `.obs` columns include `is_cell`, `cell_call_source`,
  `call_provenance`, `gex_module_call`, `atac_module_call`,
  `effective_atac_module_call`, `gex_rescue_eligible`, `atac_low_targeting`,
  and `arc_is_cell`.
- `mdata.uns["multiome"]["y_removal_enabled"] == "true"`.

## Tested JAX One-Lane STAR/Chromap Smoke

Completed on 2026-05-17 in:

```bash
tests/jax_multiome_lane_smoke_20260517T052512Z
```

This is the first tested STAR-side Multiome path that uses STAR GEX and
in-process Chromap ATAC from raw FASTQs rather than ARC MEX for the ATAC side.

Inputs:

- GEX lane: `25E113-L13_GT25-09244_TGTCCCAACG-TGGACATCGA_S113_L008_R1/R2_001.fastq.gz`.
- ATAC lane: `25E113-L1_GT25-09222_SI-NA-G2_S8_L005_R1/R2/R3_001.fastq.gz`.
- GEX reference: `/storage/autoindex_110_44/bulk_index`.
- ATAC barcode window: bases 9-24 of the ATAC R2/i5 read, reverse-complemented
  by native Chromap read-format support (`--chromapAtacReadFormat bc:8:23:-`).
- ATAC-to-GEX barcode translation:
  `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/atac2gex.tsv`.

Command shape:

```bash
scripts/run_star_multiome_lane_smoke.sh \
  --gex-r1 /mnt/pikachu/JAX_Multiome01/raw/25E113-L13_GT25-09244_TGTCCCAACG-TGGACATCGA_S113_L008_R1_001.fastq.gz \
  --gex-r2 /mnt/pikachu/JAX_Multiome01/raw/25E113-L13_GT25-09244_TGTCCCAACG-TGGACATCGA_S113_L008_R2_001.fastq.gz \
  --atac-r1 /mnt/pikachu/JAX_Multiome01/raw/25E113-L1_GT25-09222_SI-NA-G2_S8_L005_R1_001.fastq.gz \
  --atac-barcode /mnt/pikachu/JAX_Multiome01/raw/25E113-L1_GT25-09222_SI-NA-G2_S8_L005_R2_001.fastq.gz \
  --atac-r2 /mnt/pikachu/JAX_Multiome01/raw/25E113-L1_GT25-09222_SI-NA-G2_S8_L005_R3_001.fastq.gz \
  --out-dir tests/jax_multiome_lane_smoke_$(date -u +%Y%m%dT%H%M%SZ) \
  --threads 16 \
  --chromap-threads 16 \
  --chromap-low-mem \
  --chromap-macs3-frag-low-mem \
  --chromap-start-mode concurrent \
  --remote-host 10.159.4.53 \
  --remote-root /home/lhhung/jax_multiome_remote_downstream_smoke \
  --remote-output-name downstream_genefull_velocyto_cellbender \
  --cellbender-cpu-cores 24 \
  --cellbender-gpu
```

STAR/Chromap checks:

- `--soloFeatures GeneFull Velocyto`.
- `--soloInlineHashMode no`; inline-hash mode skipped standard MEX emission in
  this checkout and should stay off for this smoke.
- Y-removal enabled with `--emitNoYBAM yes`, `--emitYNoYFastq yes`, and
  `--emitYNoYFastqCompression gz`.
- `Aligned.out_Y.bam`, `Aligned.out_noY.bam`, and paired `y_separated/` FASTQs
  were emitted.
- Chromap consumed the raw ATAC barcode FASTQ directly; the Python barcode
  normalizer is retained only behind `--normalize-atac-barcode` for fallback.
- Chromap wrote BAM plus binary sidecar, then local sidecar post-processing
  wrote peaks/summits:
  - `atac_possorted.bam`
  - `atac_fragments.bin` plus `.chroms.tsv`
  - `atac_peaks.narrowPeak`
  - `atac_summits.bed`
- ATAC peak MEX was generated by the native C++ builder
  `core/features/libchromap_contract/star_multiome_atac_peak_mex` from the
  binary sidecar; the Python/bedtools helper is not compatible with this
  production boundary.
- ATAC peak MEX:
  - 70,940 peaks
  - 241,858 barcodes
  - 25,640,564 nonzero peak-barcode entries

Remote downstream checks:

- remote host: `10.159.4.53`.
- Python helper backend used the `biodepot/scrna-matrices:latest` Docker image
  because the remote host did not have `anndata`/`scanpy` installed in host
  Python.
- CellBender ran in `biodepot/cellbender:0.3.2` and wrote
  `cellbender/cellbender_counts.h5`.
- RNA h5ads include `layers["denoised"]`.

MuData checks:

- `star_chromap_unfiltered_multiome.h5mu`: 241,858 obs, 38,606 RNA vars,
  70,940 ATAC vars.
- `star_chromap_filtered_multiome.h5mu`: 4,943 obs, 38,606 RNA vars,
  70,940 ATAC vars.
- modalities: `rna`, `atac`.
- RNA layers: `ambiguous`, `counts`, `denoised`, `spliced`, `unspliced`.
- ATAC layers: `counts`.
- `mdata.uns["multiome"]["y_removal_enabled"] == "true"`.

Native adapter status:

- STAR exposes `--chromapAtacReadFormat` and passes it through
  `ChromapAtacConfig.read_format` to Chromap. For this JAX ARC lane the
  production-shaped setting is `bc:8:23:-`, equivalent to 1-based bases 9-24
  with reverse-complement.
- `scripts/run_star_multiome_lane_smoke.sh` now defaults to the native barcode
  path and avoids materializing a normalized barcode FASTQ.
- `core/features/libchromap_contract/star_multiome_atac_peak_mex` replaces the
  Python/bedtools ATAC peak-MEX helper for production-shaped smoke runs.

## JAX_Multiome01 Production Recipe

Inputs:

- raw FASTQs: `/mnt/pikachu/JAX_Multiome01/raw`
- metadata workbook:
  `/mnt/pikachu/DPC_metadata_template_Multiome1-complete.xlsx`
- output root:
  `/mnt/pikachu/JAX_Multiome01_processed/star_multiome_<timestamp>`
- remote downstream root:
  `/home/lhhung/jax_multiome_remote_downstream_production` on `10.159.4.53`
- Globus destination endpoint: `61fb8b9a-9b52-456e-928c-30c0fb0140bf`
- Globus destination root: `/JAX_Multiome01_processed/large_files`

Production launch shape:

```bash
scripts/run_jax_multiome01_production.sh \
  --threads 16 \
  --chromap-threads 16 \
  --chromap-low-mem \
  --chromap-macs3-frag-low-mem \
  --chromap-start-mode concurrent \
  --globus-upload-large-files \
  --no-sync-images
```

Resume an existing output root after a completed sample:

```bash
scripts/run_jax_multiome01_production.sh \
  --output-root /mnt/pikachu/JAX_Multiome01_processed/star_multiome_prod_globus_20260517T183219Z \
  --threads 16 \
  --chromap-threads 16 \
  --chromap-low-mem \
  --chromap-macs3-frag-low-mem \
  --chromap-start-mode concurrent \
  --globus-upload-large-files \
  --no-sync-images \
  --skip-build \
  --start-at <sample_label_or_slug>
```

Operational notes:

- the wrapper builds STAR with Chromap support once unless `--skip-build` is
  supplied. Do not rebuild while an active production run is using the binary.
- each sample runs local STAR/Chromap through `--stop-after-local-mex`, then a
  background post-MEX worker runs the remote downstream/MuData handoff and,
  when requested, the Globus large-file upload.
- remote post-MEX workers are serialized with
  `logs/remote_post_mex.lock`, so the GPU server is not oversubscribed.
- the background post-MEX subshell must keep stdin redirected from `/dev/null`.
  Do not remove the `< /dev/null` on that block: SSH or child processes can
  otherwise consume the manifest stream and make the wrapper stop after an
  early sample.
- Globus upload waits for task success before deleting generated local BAM and
  Y/noY FASTQ files. Raw input FASTQs listed in the manifest are uploaded but
  not deleted.
- fragment sidecars are not uploaded by default. Add `--include-fragments` to
  `scripts/upload_jax_multiome01_large_files_globus.sh` only if the handoff
  needs generated ATAC sidecars outside the local result tree.

Per-sample completion markers:

- local boundary: `<sample_out>/LOCAL_MEX_READY.txt`
- remote downstream/MuData boundary: `<sample_out>/REMOTE_POST_MEX_READY.txt`
- final MuData:
  `<sample_out>/mudata/star_chromap_unfiltered_multiome.h5mu` and
  `<sample_out>/mudata/star_chromap_filtered_multiome.h5mu`
- Globus state:
  `<run_root>/logs/globus_large_files/upload_state.tsv`

## Production Readiness Criteria

Keep the production path on this shape only while all are true:

- `scripts/build_multiome_mudata.py --help` and
  `scripts/run_star_multiome_lane_smoke.sh --help` work in the production
  environment.
- `scripts/run_remote_multiome_post_mex_rsync.sh --help` works on the local
  host, and the remote host can run `mudata` directly or create the per-job
  `mudata_venv`.
- The script accepts either combined ARC-style MEX or separate STAR/Chromap
  artifacts.
- Barcode namespace normalization is explicit and tested through
  `--chromapAtacReadFormat` plus `chromapAtacBarcodeTranslate`.
- Phase 2/3 expression-side tests preserve STAR Velocyto layers and validate
  the selected Y-removal mode.
- The ATAC peak-MEX adapter for STAR/libchromap output is native C++ and has
  parity checks against the Python/bedtools smoke output.
- The 100K fixture writes valid `unfiltered_multiome.h5mu` and
  `filtered_multiome.h5mu`.
- The JAX one-lane and production samples write valid `.h5mu` files with
  required RNA layers (`counts`, `spliced`, `unspliced`, `ambiguous`, and
  `denoised` when CellBender ran) plus ATAC `counts`.
- The production wrapper keeps the post-MEX background subshell stdin redirected
  from `/dev/null` so manifest iteration cannot be drained by `ssh`.
- Generated outputs stay under `tests/multiome_mudata_smoke_output_*`,
  `tests/jax_multiome_lane_smoke_*`,
  `/mnt/pikachu/JAX_Multiome01_processed/star_multiome*`, `/tmp`, or another
  untracked artifact root documented in `tests/ARTIFACTS.md`.

## Open Implementation Tasks

1. Add postflight summaries similar to UCSF/MSK `summary.txt`, including modality
   shapes, count sums, cell counts, and source artifacts.
2. Add a regression harness for the production wrapper manifest loop so
   background remote workers cannot inherit and consume the loop stdin.
3. Promote the STAR/libchromap ATAC metrics adapter into the multiome combiner
   once the desired final GEX/ATAC cell-call policy is agreed.
4. Add a compact post-MEX validation report for each `.h5mu`, including layer
   names, modality shapes, barcode intersections, and Y-removal provenance.

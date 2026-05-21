# Handoff: Integrate the Adaptive mt% QC Guard into the scRNA Downstream Pipeline

**Date:** 2026-05-18
**Checkout:** `/mnt/pikachu/STAR-suite`
**New files committed/staged:** see "Files added" below — **untracked, not yet committed**
**Integration target:** `scripts/run_scrna_downstream_gene_full_velocyto.sh`

## Summary

The scRNA downstream cell-QC filter uses a fixed `mt_pct <= 5%` cutoff. For
samples with a biologically elevated baseline mt% (MSK 30KO ES: median mt% of
called singlets ≈ 5.07%) the fixed cutoff bisects a healthy population —
measured: ~93% of dropped ES cells fail `mt% > 5`, only ~10% fail the n_genes
cap.

A replacement **combined strict-floor + MAD soft guard** has been written,
tested, and documented as standalone STAR-suite helper scripts:

```
mt_threshold = max(mt_floor 5%, median(mt_pct) + n_mad·MAD(mt_pct))   # per sample
keep cell  <=>  mt_pct <= mt_threshold                                 # combined guard
```

Method rationale, the 2.5-vs-3 nMADs discussion, raw-vs-scaled MAD, and the
sample-level flag are all in the runbook:
`docs/RUNBOOK_SCRNA_MT_ADAPTIVE_FILTER_20260518.md`.

**This handoff is the integration task** — wire the new scripts into the
pipeline orchestrator and run the one-time conversion of existing releases.
The scripts themselves are done; do not rewrite them.

## Files added (done — tested, do not modify without reason)

| File | Role |
| --- | --- |
| `scripts/scrna_mt_adaptive.py` | Shared core module: `compute_mt_threshold`, `mt_sample_flag`, `qc_filter_mask`, `apply_mt_adaptive_filter`, `merge_threshold_json`. Imported by the two CLI scripts below — keeps both paths identical. |
| `scripts/apply_adaptive_mt_filter.py` | NEW-IMPLEMENTATION CLI. Recomputes `filter` / `singlet_filtered` obs on the unfiltered downstream h5ad; merges mt keys into `adaptive_qc_threshold.json`. |
| `scripts/convert_h5ad_mt_adaptive.py` | CONVERSION CLI. Retrofits an already-built sample directory. |
| `scripts/generate_qc_histogram_mt_adaptive.py` | QC histogram with two mt rejection lines (floor + adaptive). Replaces `generate_qc_histogram.py`. |
| `docs/RUNBOOK_SCRNA_MT_ADAPTIVE_FILTER_20260518.md` | Method + procedures. |

All four scripts `py_compile` clean. `convert_h5ad_mt_adaptive.py --dry-run`
and `generate_qc_histogram_mt_adaptive.py` were run against real MSK archive
h5ads and verified (ES → threshold 7.69%, +12,433 cells; PP1 → floored to 5%,
0 delta; S6_1 → threshold 10.01%, flagged).

## Verified behaviour (dry-run, MSK archive)

| Sample | mt_threshold | filter cells strict-5% → adaptive | flag |
| --- | ---: | --- | --- |
| ES | 7.69% | 18,187 → 30,620 | no |
| PP1 | 5.00% (floored) | 23,708 → 23,708 | no |
| S6_1 | 10.01% | 18,902 → 24,376 | **yes** |

## Task 1 — Integrate into the pipeline orchestrator

File: `scripts/run_scrna_downstream_gene_full_velocyto.sh`. Three edits. **Do
not modify `combineFilters.py`** — it still runs with its strict-5% cutoff;
the apply step overrides the `filter` column afterward and preserves the old
one as `filter_strict_mt5`.

### Edit 1 — declare the new script paths

Near the existing `GENERATE_QC_HISTOGRAM=...` declaration (~line 13), add:

```bash
APPLY_ADAPTIVE_MT="${REPO_ROOT}/scripts/apply_adaptive_mt_filter.py"
GENERATE_QC_HISTOGRAM_MT="${REPO_ROOT}/scripts/generate_qc_histogram_mt_adaptive.py"
```

### Edit 2 — apply the adaptive mt guard after combineFilters.py

Immediately **after** the `combineFilters.py` Docker call and **before** the
QC-histogram block, insert (gated on `ADAPTIVE_FILTER`):

```bash
if [[ "${ADAPTIVE_FILTER}" == "1" ]]; then
  run_py "${APPLY_ADAPTIVE_MT}" \
    --input-h5ad     "${UNFILTERED_H5AD}" \
    --threshold-json "${ADAPTIVE_QC_JSON}" \
    --mt-floor       "${MT_PCT_CUTOFF}" \
    --n-mad          "${N_MAD}"
fi
```

- Runs on the host backend via `run_py` (needs `anndata`/`numpy` — same as
  `compute_adaptive_qc_threshold.py` / `postprocess_downstream_filters.py`).
- Reads `min_genes` and `effective_max_genes` from `${ADAPTIVE_QC_JSON}`
  (written earlier by `compute_adaptive_qc_threshold.py`) — so this step must
  run after that one. It already will, given placement.
- Recomputes `filter` / `singlet_filtered` in `${UNFILTERED_H5AD}` in place,
  so the existing `postprocess_downstream_filters.py` call needs **no change**
  — it builds the filtered views from the recomputed `filter`.
- The existing `MT_PCT_CUTOFF` env var (default 5) is reused as the strict
  floor; `N_MAD` (default 3) as the MAD multiplier.

### Edit 3 — swap the QC histogram script

Replace the `generate_qc_histogram.py` invocation inside the
`if [[ "${ADAPTIVE_FILTER}" == "1" ]]; then` histogram block with the
two-line version. Old call:

```bash
  docker "${DOCKER_ARGS[@]}" "${DOCKER_IMAGE}" \
    python3 "${GENERATE_QC_HISTOGRAM}" \
      --input-h5ad "${UNFILTERED_H5AD}" --output-dir "${OUTPUT_DIR}" \
      --min-genes "${MIN_GENES}" --max-genes "${EFFECTIVE_MAX_GENES}" \
      --mt-pct-cutoff "${MT_PCT_CUTOFF}" --n-mad "${N_MAD}" \
      --raw-adaptive-max "${RAW_ADAPTIVE_MAX_GENES}"
```

New call (same Docker wrapper; new arg set — reads everything from the JSON):

```bash
  docker "${DOCKER_ARGS[@]}" "${DOCKER_IMAGE}" \
    python3 "${GENERATE_QC_HISTOGRAM_MT}" \
      --input-h5ad     "${UNFILTERED_H5AD}" \
      --output-dir     "${OUTPUT_DIR}" \
      --threshold-json "${ADAPTIVE_QC_JSON}"
```

- The new script needs `mt_pct_floor` / `mt_pct_raw_threshold` /
  `mt_pct_threshold` / `mt_pct_n_mad` in the JSON — Edit 2 writes those, so
  Edit 3 must come after Edit 2's step at runtime (it does — histogram block
  is already after combineFilters.py).
- The script always writes `gene_quantile_histogram.{png,html}`. When
  `mt_pct_flag` is `true` in the JSON it **additionally** writes
  `mt_quantile_histogram.{png,html}` (a standalone mt% distribution for
  human review). Flagged samples therefore have one extra artifact — make
  sure any staging/upload step globs both `*_histogram.{png,html}` rather
  than the single `gene_quantile_histogram` filename.
- Confirm the Docker image (`biodepot/scrna-matrices:latest`) has `plotly` +
  `kaleido` (it ran the old histogram script, so it should). If `run_py`'s
  backend also has plotly, running it via `run_py` instead of Docker is fine
  too.
- Leave `generate_qc_histogram.py` in the repo for the non-adaptive
  (`ADAPTIVE_FILTER=0`) path, or point both paths at the new script — your
  call; the new script only requires the JSON to exist.

### Ordering after integration

```
compute_adaptive_qc_threshold.py     # n_genes adaptive max  -> adaptive_qc_threshold.json
combineFilters.py (Docker)           # mt_pct/n_genes/singlet + legacy strict-5% filter
apply_adaptive_mt_filter.py          # NEW: overrides filter/singlet_filtered with mt guard
generate_qc_histogram_mt_adaptive.py # NEW: QC graph, two mt lines
postprocess_downstream_filters.py    # builds filtered_counts / default_singlet (unchanged)
```

No change is needed in `scripts/run_msk30ko_pipeline_from_manifest.py` — it
only passes `--adaptive-filter`, which sets `ADAPTIVE_FILTER=1`.

### Module-path note

`apply_adaptive_mt_filter.py` and `generate_qc_histogram_mt_adaptive.py` do
`import scrna_mt_adaptive`. That resolves because the module sits in the same
`scripts/` directory (script dir is on `sys.path[0]`). If the Docker call for
the histogram mounts only individual files, ensure `scrna_mt_adaptive.py` is
present in `scripts/` inside the container mount — the existing
`-v "${OUTPUT_DIR}:${OUTPUT_DIR}"` plus repo mount should already cover it;
verify the repo `scripts/` dir is mounted for the Docker histogram call.

## Task 2 — One-time conversion of existing releases

Separate from the pipeline edit. Retrofits the already-shipped MSK 30KO
release. Conversion needs each sample's all-barcodes `final_counts.h5ad`
(production-sample archive), **not** the delivery-tree `filtered_counts.h5ad`.

```bash
# preview
python3 scripts/convert_h5ad_mt_adaptive.py --dry-run \
  --sample-dir /mnt/pikachu/msk30ko-production-sample-archive/MSK-05-13-26-large-files/{ES,DE,PP1,PP2,S5_1,S5_2,S6_1,S6_2}/downstream_genefull_velocyto_cellbender

# convert to a fresh tree
python3 scripts/convert_h5ad_mt_adaptive.py \
  --sample-dir <dirs...> --output-dir /mnt/pikachu/msk30ko-mt-adaptive-rebuild

# regenerate QC graphs per converted sample
python3 scripts/generate_qc_histogram_mt_adaptive.py \
  --input-h5ad <out>/final_counts.h5ad --output-dir <out> \
  --threshold-json <out>/adaptive_qc_threshold.json
```

Then re-stage `filtered_counts.h5ad`, `default_singlet_filtered_counts.h5ad`,
`adaptive_qc_threshold.json`, and `gene_quantile_histogram.{png,html}` into the
delivery tree and re-upload. Flagged samples (`mt_pct_flag = true`) also have
`mt_quantile_histogram.{png,html}` — stage it too. Full procedure in the
runbook §"Conversion procedure".

**`DE_GemX`** is the 9th sample and is **not** in the 8-sample production
archive — it was processed in an earlier separate trial run and not re-run.
Convert it directly from that run directory (it has `final_counts.h5ad` +
`adaptive_qc_threshold.json`):

```
/storage/MSK-perturb-comparison/msk30ko_DE_GemX_velocyto_downstream_trial_20260513_190549/samples/DE_GemX/downstream_genefull_velocyto_cellbender/
```

## Validation / acceptance criteria

- Pipeline: run one MSK sample end-to-end; confirm `adaptive_qc_threshold.json`
  gains the `mt_pct_*` keys, `final_counts.h5ad` obs gains `filter_strict_mt5`
  / `singlet_filtered_strict_mt5`, and `gene_quantile_histogram.png` shows two
  mt lines.
- For every sample: `filter_cells_mt_adaptive >= filter_cells_strict_mt5`
  (the guard never drops cells the 5% cutoff kept).
- `mt_pct_threshold_was_floored` is `true` for clean samples (PP1/PP2),
  `false` for ES/S5_*/S6_*.
- Review any sample with `mt_pct_flag = true` (expected: S6_1) against its
  `mt_quantile_histogram.{png,html}` — emitted only for flagged samples —
  before release.
- Re-run the provider cell-call concordance after conversion; ES recovery
  should rise toward the provider call (some overshoot is expected and OK).

## Risks & caveats

- `combineFilters.py` runs unchanged and still computes its own strict-5%
  `filter`; `apply_adaptive_mt_filter.py` overwrites it. The legacy column is
  preserved as `filter_strict_mt5` — do not drop it (rollback + audit).
- The apply/convert steps rewrite the full unfiltered h5ad (≈2M barcodes).
  Expect a multi-hundred-MB rewrite per sample; ensure disk headroom.
- `convert_h5ad_mt_adaptive.py` writes `*.pre_mt_adaptive` backups unless
  `--no-backup`; with `--output-dir` it writes to a fresh tree and leaves
  originals untouched — prefer `--output-dir` for the release retrofit.
- The threshold is computed over `obs['singlet']` (STAR-called singlets),
  matching the n_genes adaptive convention. If `singlet` is absent the apply
  step raises `KeyError` — combineFilters.py must have run first.

## Suggested commit plan

Current branch `feature/native-multiome-atac-barcode` is unrelated. Branch
from `master`:

```
git checkout master && git checkout -b feature/scrna-mt-adaptive-filter
git add scripts/scrna_mt_adaptive.py scripts/apply_adaptive_mt_filter.py \
        scripts/convert_h5ad_mt_adaptive.py scripts/generate_qc_histogram_mt_adaptive.py \
        docs/RUNBOOK_SCRNA_MT_ADAPTIVE_FILTER_20260518.md \
        docs/HANDOFF_SCRNA_MT_ADAPTIVE_FILTER_INTEGRATION_20260518.md
# commit 1: the helper scripts + runbook + handoff
# commit 2: the run_scrna_downstream_gene_full_velocyto.sh integration (Task 1)
```

No `Co-Authored-By` trailers (repo policy).

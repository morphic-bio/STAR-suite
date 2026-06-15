# CITE-seq MEX Validation Runbook

Date: 2026-06-15

## Claim

STAR-suite supports CITE-seq-style feature-barcode assays through the MEX stage:

- ADT / protein tags are emitted as `Antibody Capture`.
- HTO / CMO / sample-hash tags are emitted as `Multiplexing Capture`.
- Mixed ADT+HTO references can split protein and hash MEX outputs while keeping
  the merged 10x-style feature-barcode matrix compatible with downstream tools.

This runbook validates the raw feature-barcode processing contract. It does not
claim downstream protein normalization, dsb/CLR, WNN, clustering, or biological
interpretation.

**TotalSeq-C validation warning:** the formal public 10x benchmark in this
runbook is the PBMC 1k **TotalSeq-B** dataset. STAR-suite emits the generic
10x `Antibody Capture` MEX contract used by both TotalSeq-B and TotalSeq-C
feature-barcode libraries, but we have not yet run STAR-suite against a
TotalSeq-C raw-read dataset. Do not claim TotalSeq-C dataset-level validation
until that stress test is added and recorded.

## Tiers

| Tier | Purpose | Command |
|---|---|---|
| A | No-download synthetic regression | `bash tests/run_citeseq_validation_suite.sh --quick` |
| B | Public DOGMA-HIV real HTO/ADT smoke, using local staged files | included in `--quick`; skips if local data are absent |
| C | Public 10x PBMC 1k TotalSeq-B ADT benchmark | `bash tests/run_citeseq_pbmc1k_pf_benchmark.sh --download` |

Tier C is the formal public 10x benchmark because it has raw FASTQs, feature
reference, and Cell Ranger MEX outputs available from 10x. TotalSeq-C is
deferred until a usable raw-read benchmark is selected.

## Public Benchmark: PBMC 1k TotalSeq-B v3.1

Dataset:

- 10x page: <https://www.10xgenomics.com/datasets/1-k-human-pbm-cs-stained-with-a-panel-of-total-seq-b-antibodies-single-indexed-3-1-standard-4-0-0>
- Sample ID: `SC3_v3_NextGem_SI_PBMC_CSP_1K`
- Cell Ranger version: 4.0.0
- Chemistry: Chromium Next GEM Single Cell 3' v3.1, single-indexed
- Assay: GEX + 17 TotalSeq-B cell-surface protein tags

Pinned assets:

```text
https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_CSP_1K/SC3_v3_NextGem_SI_PBMC_CSP_1K_fastqs.tar
https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_CSP_1K/SC3_v3_NextGem_SI_PBMC_CSP_1K_feature_ref.csv
https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_CSP_1K/SC3_v3_NextGem_SI_PBMC_CSP_1K_filtered_feature_bc_matrix.tar.gz
https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_CSP_1K/SC3_v3_NextGem_SI_PBMC_CSP_1K_raw_feature_bc_matrix.tar.gz
```

The FASTQ tar is about 6.3 GB. Keep it outside the repository.

### One-command benchmark

```bash
bash tests/run_citeseq_pbmc1k_pf_benchmark.sh --download
```

The wrapper:

1. stages the public fixture under `/mnt/pikachu/citeseq_pbmc1k_totalseqb_v31`
   unless `CITESEQ_PBMC1K_FIXTURE` or `--fixture` points elsewhere;
2. runs `assignBarcodes --output-mode adt_mex` on the ADT FASTQs only;
3. uses the official Cell Ranger filtered barcodes as the filtered cell
   universe, materialized as a plain TSV with Cell Ranger `-1` suffixes removed
   for `assignBarcodes`;
4. normalizes the Cell Ranger TRU barcode list into the NXT assignment namespace
   used by the ADT FASTQs in this v3.1 fixture, then translates STAR-suite
   barcodes back to TRU for comparison;
5. compares STAR-suite `Antibody Capture` MEX counts to the official Cell
   Ranger filtered feature-barcode MEX.

Default local assignment whitelist:

```bash
export CITESEQ_NXT_WHITELIST=/storage/scRNAseq_output/whitelists/3M-february-2018_NXT.txt
```

Override it when running on another host. The wrapper defaults to
`CITESEQ_PBMC1K_SOURCE_NAMESPACE=TRU` and
`CITESEQ_PBMC1K_TARGET_NAMESPACE=NXT`. This is fixture-specific: the official
Cell Ranger filtered MEX reports TRU barcodes, while the antibody FASTQs match
the NXT-split whitelist. The NXT/TRU transform is the reversible two-base
STAR-suite namespace transform, not a biological barcode correction.

### Manual staging

```bash
bash scripts/download_10x_citeseq_pbmc1k_fixture.sh \
  --outdir /mnt/pikachu/citeseq_pbmc1k_totalseqb_v31
```

Then run without network access:

```bash
bash tests/run_citeseq_pbmc1k_pf_benchmark.sh \
  --fixture /mnt/pikachu/citeseq_pbmc1k_totalseqb_v31
```

Expected benchmark outputs:

```text
/mnt/pikachu/citeseq_pbmc1k_pf_benchmark_*/
  assignBarcodes.log
  compare.log
  citeseq_pbmc1k_compare.json
  citeseq_pbmc1k_feature_totals.tsv
  BENCHMARK_SUMMARY.txt
```

The comparison gate records:

- selected `Antibody Capture` feature count;
- filtered barcode overlap;
- total ADT UMI delta;
- per-feature total Pearson;
- per-barcode ADT total Pearson;
- sparse entry-level absolute deltas on common feature/barcode axes.

Observed local result after namespace normalization:

```text
features: Cell Ranger=32 STAR=32 common=32
barcodes: Cell Ranger=1206 STAR=1206 common=1206
total ADT UMIs: Cell Ranger=4,080,957 STAR=4,096,584 rel_delta=0.0038
feature-total Pearson=0.9999989
barcode-total Pearson=0.9987676
```

Default pass thresholds are intentionally looser than the observed result, to
allow for cross-tool UMI handling differences while still catching broken MEX
contracts:

```text
feature-total Pearson >= 0.95
barcode-total Pearson >= 0.90
total relative delta <= 0.25
missing ADT features from STAR = 0
```

Tighten thresholds only if repeated runs on the same fixture remain stable
across clean builds.

## Public Real Dataset: DOGMA-HIV

DOGMA-HIV is a public real mixed ADT+HTO dataset. The current STAR-suite smoke
uses a local staged copy of the public data rather than auto-downloading it. It
starts from raw ADT/HTO FASTQs and exercises native STAR-suite HTO demux:

```bash
bash tests/multi_feature/test_hiv_dogma_hto_demux_smoke.sh
```

The script uses local staged FASTQ paths by default and skips cleanly when they
are absent. It starts from raw ADT/HTO FASTQs, but it is not yet the complete
end-to-end STAR-suite run from all DOGMA libraries; that remains a separate
validation after the PBMC benchmark.

Current local smoke output on the staged YW8 filtered barcode surface:

```text
hash_features=6
protein_features=165
singlet=3112
doublet=113
negative=5856
```

The script still prints non-fatal notes against older demux guidance targets;
those notes are not acceptance criteria for this CITE-seq MEX validation. Use
DOGMA here as the real mixed-library contract smoke. The public 10x PBMC 1k
fixture is protein-only; DOGMA covers the mixed protein+hash path until the full
end-to-end DOGMA validation is run.

## Synthetic Regression Gates

Run the quick tier before touching feature-barcode code:

```bash
bash tests/run_citeseq_validation_suite.sh --quick
```

This runs:

```bash
bash core/features/process_features/tests/test_adt_mex.sh
bash core/features/process_features/tests/test_hash_demux_mex.sh
bash tests/multi_feature/test_adt_protein_multifeature_arm.sh
bash tests/multi_feature/test_hto_cmo_hash_demux_arm.sh
bash tests/multi_feature/test_hiv_dogma_hto_demux_smoke.sh
```

The first four are synthetic and should always pass in a normal build
environment. The DOGMA smoke may skip when local data are unavailable.

## TotalSeq-C Status

Current 10x protein-tag datasets are no longer always labeled "CITE-seq"; they
are presented as Cell Surface Protein / TotalSeq / Feature Barcode products.
That naming does not change STAR-suite's MEX contract: ADT/protein rows should
still emit as `Antibody Capture`.

Do not list TotalSeq-C as validated yet. The current validation set covers:

- public TotalSeq-B PBMC raw-read ADT assignment;
- synthetic hash-only and mixed ADT+HTO contracts;
- DOGMA-HIV real mixed ADT+HTO smoke.

The next TotalSeq-C validation should be treated as an added Tier D stress test,
not as part of the current acceptance set. It needs raw feature FASTQs, a
feature reference, and a Cell Ranger MEX oracle.

## Acceptance Criteria

- Synthetic ADT and HTO tests pass.
- DOGMA-HIV raw ADT FASTQ smoke either passes or skips because local data are
  absent; it must not silently fail.
- Public PBMC 1k benchmark emits STAR-suite ADT MEX with `Antibody Capture`
  rows and passes the comparison thresholds.
- Documentation states clearly that no TotalSeq-C raw-read dataset-level
  validation has been completed yet.
- No generated FASTQs, MEX outputs, logs, or downloaded 10x assets are committed.
- Any retained benchmark output path is documented in `tests/ARTIFACTS.md`.

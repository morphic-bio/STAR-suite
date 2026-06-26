# STAR Suite vs `nf-core/rnaseq` star_salmon — output parity check

A small, self-contained harness that runs **both** pipelines on the **same real
human reads** and the **same reference**, emits **every intermediate**, and shows
that STAR Suite's outputs **match** the reference chain — per-transcript and
per-gene quant concordance, plus BAM and read-count agreement. The point here is
**parity**: that swapping in STAR Suite changes the *speed*, not the *answer*.

**Parity vs speed — two separate questions, two separate artifacts.** This harness
answers parity, and it does so on a fast `chr22` slice because parity is
depth-robust (both pipelines see the same reads, so agreement reflects the method,
not coverage). **Speed** is a full-dataset question — see STAR Suite's published
full-genome benchmarks, which anyone can re-run; timing is reported here too but is
**indicative only** (chr22 is too small to represent the integrated-pass advantage).

We know maintainer time is scarce, so this is built to be quick to read, quick to
run, and easy to point at your own data.

## What it compares

- **A — reference chain:** Trim Galore → STAR (2.7.x) → Salmon — the exact tools
  `nf-core/rnaseq --aligner star_salmon` runs. Its Salmon `quant.sf` is the baseline.
- **B — STAR Suite (integrated):** one binary doing trim → align → quant, emitting
  Salmon-format `quant.sf` + `quant.genes.sf` (tximport-ready).
- **B' — STAR Suite BAM → Salmon:** STAR Suite's transcriptome BAM handed to *your*
  Salmon, for exact `star_salmon`-equivalent numbers.

You get speed (A vs B) and two honest concordance views (B is STAR Suite's own
quant; B' isolates the alignment by using the same Salmon).

## Run it

```bash
# fast: real human reads aligned to a chr22-only reference (minutes)
./run_compare.sh --mode chr22

# full GRCh38 (the production number; longer)
./run_compare.sh --mode full
```

Outputs land in `runs/<mode>_<timestamp>/`; the Markdown report is printed at the
end and saved to `compare/report.md`. Every intermediate (trimmed FASTQs, genome
and transcriptome BAMs, all `quant.sf`) is kept on disk for inspection.

## Current chr22 TranscriptVB Parity

The 2026-06-26 parity artifact is:

```bash
benchmarks/nfcore_rnaseq_compare/runs/chr22_20260626_003522/
```

The cached real-read fixture is small enough to keep locally for smoke testing:

```bash
benchmarks/nfcore_rnaseq_compare/fixtures/chr22/
```

Latest fixed-parameter checks on that fixture. The auto-detect row is a
representative pass of `tests/run_transcriptvb_chr22_parity_smoke.sh`; exact last
digits can vary slightly because STAR's threaded online FLD consumes alignment
chunks in thread-completion order.

| Comparison | NumReads Pearson | TPM Pearson | NumReads Spearman | total NumReads | half-L1 moved |
|---|---:|---:|---:|---:|---:|
| STAR auto `A`, no error vs Salmon fixed `ISR`, no error, `-p 1` | `0.999990519803` | `0.999990675522` | `0.999107439749` | `12771.994` vs `12771.999` | `13.034` |
| STAR fixed `ISR`, no error vs Salmon fixed `ISR`, no error, `-p 1` | `0.999995116948` | `0.999996544483` | `0.999999425373` | `12772.004` vs `12772.005` | `9.351` |
| STAR fixed `ISR`, no error, no-FLD/no-effective-correction vs Salmon same flags | `0.999999999990` | `0.999999999998` | `1.000000000000` | `12771.996` vs `12771.998` | `0.010` |

`half-L1 moved` is the approximate number of fractional read-equivalents assigned
to a different transcript after multimapper allocation. In the representative
auto-detect smoke, only seven transcripts differ by more than one read-equivalent,
with a maximum single-transcript delta of `3.819`.

### Salmon Thread-Order Caveat

Salmon's alignment-mode online fragment-length distribution is thread-order
sensitive on this small chr22 fixture. Salmon `-p 1` and Salmon `-p 32` disagree
materially here (`NumReads` Pearson about `0.9957`) and show the same large outlier
pattern that initially looked like a STAR Suite gap.

The cause is ordering, not a different BAM. Salmon updates its alignment-mode FLD
online while consuming alignment records; with multiple worker threads, early
alignment observations can reach the FLD estimator in a different order. On this
small fixture that early-order difference is large enough to alter effective-length
weights and VB allocation for a small number of multimapping transcripts. STAR
Suite's integrated quantification also builds an online FLD from threaded
alignment chunks, so changing `--runThreadN` can move a small number of fractional
read-equivalents in this tiny smoke. The integration target is therefore fixed as
STAR at 32 threads, matching the retained chr22 artifact, versus Salmon `-p 1`.
`run_compare.sh --mode chr22` defaults Salmon quantification to one thread while
leaving STAR threaded; override with `--salmon-threads` or `SALMON_THREADS` when
intentionally studying Salmon's threaded behavior.

### Retained Artifact And Smoke Test

Keep the 2026-06-26 chr22 run as a local retained artifact because it contains the
diagnostic intermediate outputs used to explain the final residual differences.
For routine integration, use the same-BAM smoke wrapper instead of checking in the
full exploratory run tree.

The same-BAM smoke wrapper is:

```bash
tests/run_transcriptvb_chr22_parity_smoke.sh
```

It reuses the cached chr22 fixture/reference, runs STAR TranscriptVB auto-detect
with 32 STAR threads by default, runs Salmon on the exact STAR transcriptome BAM
with `-p 1`, and gates on the correlation/read-movement thresholds above.

## Point it at your own data / tools

All paths are CLI- or env-overridable — nothing is hardwired beyond defaults:

```bash
./run_compare.sh --mode chr22 \
  --reads-dir /path/with/SAMPLE_1.fastq.gz_and_2 \
  --threads 16 --spots 3000000
STAR_SUITE_BIN=... STAR_UPSTREAM_BIN=... SALMON_BIN=... ./run_compare.sh --mode full
```

The default fast reads are a public human PE bulk run (GSE88509 / GSM2344101 /
SRR4422207); the chr22 reference is built from a self-consistent GRCh38 genome+GTF.

## Notes on fairness

- Same reference, same reads, same thread count for both pipelines.
- In `chr22` mode, Salmon quantification defaults to one thread because the small
  fixture exposes Salmon's online-FLD thread-order sensitivity. This is a stability
  choice for parity, not a performance comparison; full mode defaults Salmon back
  to the requested thread count.
- `chr22` mode downsamples by restricting the reference to chromosome 22 — real
  human reads, only chr22 alignments — so it runs fast while staying real. Use
  `--mode full` for the production-scale number.
- Trimming: STAR Suite trims internally (cutadapt-compatible); the reference chain
  uses Trim Galore. The report's "input reads into alignment" row is the post-trim
  count for each — close agreement indicates equivalent trimming, and the trimmed
  FASTQs are kept so you can diff them directly.
- We're happy to adjust parameters to match your pipeline's exact settings.

## Version matching (important for a fair comparison)

The reference baseline uses **upstream STAR 2.7.11b** — the exact version STAR Suite
is built on (`STAR --upstream-version`) **and** the version `nf-core/rnaseq` pins
(`bioconda::star=2.7.11b`). Matching it means any alignment difference reflects STAR
Suite's changes, not an incidental upstream version bump. The harness checks this and
refuses to run on a mismatch; set `STAR_UPSTREAM_BIN` to your 2.7.11b binary if your
PATH `STAR` is a different version.

## Requirements

Upstream **STAR 2.7.11b** (on PATH or via `STAR_UPSTREAM_BIN`), the STAR Suite `STAR`
binary, `salmon`, `trim_galore`/`cutadapt`, `samtools`, `gffread`, and `sra-tools`
(`fastq-dump`) for the default public reads. `python3` + `numpy` for the report.

# Handoff: UCSF GEX-only no-BAM benchmark wrapper (planning agent)

**Date:** 2026-03-26 (updated after unblock validation)  
**Repo:** `STAR-suite`  
**Script:** `scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh`

Use this document when updating overnight runbooks, assigning benchmark agents, or reconciling “internal gzip only” policy with the historical `7a7fb08` arm.

**Stale checkout warning:** Symptom triage matches this handoff but the run still fails the way an *older* wrapper would (pre-created `tmp`, missing `zcat` on historical, or silent wrong mode)—the agent may be on a **branch/worktree that never picked up the wrapper fixes**, or edited the handoff without patching the script. The canonical script lives only at `scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh` in this repo; **do not assume** the benchmark-only clone matches without `git diff` / verification below.

---

## Verify live script matches this handoff (mandatory for agents)

Run from repo root (line numbers in docs **drift**; use `grep`):

```bash
f=scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh
grep -q 'MODE_FLAG_COUNT' "$f"
grep -qE 'MODE_FLAG_COUNT.*-ne 1' "$f"
grep -q 'rm -rf "${OUTDIR}/tmp"' "$f"
# Exactly one non-comment code line with --readFilesCommand zcat (historical CMD only)
test "$(grep -cE '^[[:space:]]*--readFilesCommand zcat' "$f")" -eq 1
grep -q 'HANDOFF_UCSF_GEXONLY_BENCHMARK_WRAPPER_20260326' "$f"   # invariants header in script
```

**Expected implementation (conceptual, not line numbers):**

- After input validation: `mkdir -p "${OUTDIR}"` then **`rm -rf "${OUTDIR}/tmp"`** only (no `mkdir …/tmp`).
- **`--historical-vanilla` CMD array includes `--readFilesCommand zcat`** immediately after `--readFilesIn`.
- **`--modern-optimized` CMD has no `--readFilesCommand`.**
- Mode parsing uses **`MODE_FLAG_COUNT`** and **`if [[ "${MODE_FLAG_COUNT}" -ne 1 ]]`** then exit 1.

The script header comments document the same four invariants for human review.

---

## Executive summary (planning)

| Topic | Decision |
|--------|-----------|
| Two STAR modes | **`--historical-vanilla`** vs **`--modern-optimized`** — different CLI families; never merge flags. |
| Mode flags | **Exactly one** required; passing both or duplicating the same flag **exits 1** (no silent last-wins). |
| Historical + `.gz` | **`7a7fb08` needs `--readFilesCommand zcat`** for gzipped FASTQs; otherwise STAR reads gzip bytes as FASTQ and fails immediately at mapping. |
| Modern + `.gz` | **No zcat** — native internal gzip streaming. |
| `--outTmpDir` | Wrapper **`rm -rf "${OUTDIR}/tmp"`** before STAR; **do not** pre-create `tmp` (breaks at least `7a7fb08`: “could not make temporary directory … tmp/”). |
| Timing tables | Historical vs modern wall time is **not** I/O-comparable unless you footnote **zcat** on the historical arm or use uncompressed inputs (usually impractical). |

---

## Purpose

Single entry point for **full-sample** UCSF GEX Solo runs (no `pfMulti`, no guide FASTQs), `--outSAMtype None`.

**Compressed FASTQs**

- **`--modern-optimized`:** STAR native `.gz` streaming — **no** `--readFilesCommand zcat`.
- **`--historical-vanilla`:** `7a7fb08`-class builds treat `.fastq.gz` paths as **uncompressed** unless `--readFilesCommand` is set. The wrapper uses **`zcat`**, matching `tests/run_ucsf_solo_gex_100k_benchmark.sh`.

Overnight docs that mandate “internal gzip only” for **all** STAR runs should **carve out an explicit exception** for the historical vanilla arm, or plan to decompress to plain FASTQ (not recommended at full scale).

---

## Modes (non-interchangeable)

| Mode | Use case | Typical STAR build |
|------|-----------|---------------------|
| `historical-vanilla` | Legacy Solo surface aligned with the 100k harness | e.g. commit `7a7fb08` |
| `modern-optimized` | Bridge + inline hash + CR-compat GeneFull flags | current optimized / `004c36e9` line |

Do **not** expect matching Solo statistics between modes: historical uses **Rescue** multimappers and **Unstranded**; modern uses **Unique**, **Forward**, `--clip3pPolyG yes`, `--soloCrMultimapRescue`, bridge hash, `dynamicThread*`, etc. **`7a7fb08` does not parse** the modern CR-compat / bridge CLI; that is why the split exists.

---

## Argument handling (for automation)

- Parser increments a counter for each `--historical-vanilla` or `--modern-optimized`.
- If `count != 1` after parsing: print error and **exit 1** (covers “both flags” and “`--historical-vanilla` twice”).
- `--outdir` is still required.

---

## Authoritative references

- `docs/HANDOFF_SOLO_OPTIMIZATION_20260324.md` — `7a7fb08` worktrees, 100k/2M artifacts.
- `tests/run_ucsf_solo_gex_100k_benchmark.sh` — historical Solo flag shape + **`zcat`** for `.gz`.
- Related benchmark-process notes (if present in your branch): `docs/HANDOFF_SOLO_OVERNIGHT_BENCHMARK_BLOCKERS_20260326.md` — `parametersDefault.xxd` after `make clean`, older-SHA build pitfalls.

---

## Usage

Run **modern optimized first**; run **historical vanilla last** in the overnight
plan (after optimized GEX-only, CR9 if applicable, and the perturb no-BAM
matrix). Do **not** subject the historical vanilla invocation to a short
automation timeout — allow multi-hour wall time.

```bash
# Modern optimized (current / bridge binary) — run this before historical vanilla
./scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh \
  --outdir /storage/.../ucsf_gexonly_no_bam/star_optimized \
  --star-bin /path/to/STAR_optimized \
  --modern-optimized

# Historical vanilla (7a7fb08-class binary) — run last; no short timeout
./scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh \
  --outdir /storage/.../ucsf_gexonly_no_bam/star_vanilla_7a7fb08 \
  --star-bin /path/to/STAR_7a7fb08 \
  --historical-vanilla
```

**Required:** `--outdir` and **exactly one** mode flag.

**Optional:** `--threads N` (default 32), `--star-bin` (default `${REPO_ROOT}/core/legacy/source/STAR`).

**Env overrides**

| Variable | Default |
|----------|---------|
| `UCSF_GEXONLY_GEX_DIR` | `/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2/GEX` |
| `UCSF_GEXONLY_SOLO_WHITELIST` | `.../3M-february-2018_TRU.txt` (Cell Ranger 9 install path) |
| `UCSF_GEXONLY_GENOME_DIR` | `/storage/autoindex_110_44/bulk_index` |

---

## Historical vanilla CLI (canonical shape)

Aligned with the known-good 100k harness and validated for argv parsing + startup on **`7a7fb08`**:

- `--readFilesIn "${R2_FILES}" "${R1_FILES}"` (comma-separated lanes, `.fastq.gz`)
- **`--readFilesCommand zcat`**
- `--outFileNamePrefix`, `--outTmpDir` (see tmp note above)
- `--outSAMtype None`
- `--clipAdapterType CellRanger4`, `--alignEndsType Local`, `--chimSegmentMin 1000000`
- `--soloType CB_UMI_Simple` + CB/UMI positions + `--soloFeatures GeneFull`
- `--soloMultiMappers Rescue`, `--soloStrand Unstranded`
- `--soloUMIfiltering MultiGeneUMI_CR`, `--soloUMIdedup 1MM_CR`, `--soloCellFilter EmptyDrops_CR`, `--soloCbUbRequireTogether no`, `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts`

**Omitted** (modern-only / unsupported on historical): `--clip3pPolyG`, `--soloCrMultimapRescue`, `--soloCrGexFeature`, `dynamicThread*`, `STAR_SOLO_NONFLEX_HASH_BRIDGE`, `--soloInlineHashMode`.

---

## Modern optimized CLI (summary)

- `export STAR_SOLO_NONFLEX_HASH_BRIDGE=1`
- `--soloInlineHashMode yes`
- `--soloMultiMappers Unique`
- `--clip3pPolyG yes`, `--soloStrand Forward`
- `--soloCrGexFeature genefull`, `--soloCrMultimapRescue yes`
- `--dynamicThreadInterface 1`, `--dynamicThreadConstMapPermits`, `--dynamicThreadTelemetry 1`
- `--outSAMtype None`
- **No** `--readFilesCommand zcat`

---

## Validation status (2026-03-26 engineering check)

On host **pikachu**, with a clean worktree build of **`7a7fb08`** and EBs2_2 GEX inputs:

1. **Parameters:** historical CLI **parsed**; `Log.out` showed “Finished reading parameters” and the expected effective command line (no “unrecognized parameter” fatals from bridge/CR-compat flags).
2. **Genome:** “Finished loading the genome”.
3. **Mapping:** “started mapping” with **no** immediate read-format FATAL once **`zcat`** was enabled.
4. **Without `zcat`:** FATAL at first read — gzip bytes interpreted as FASTQ (`wrong read ID line format`).
5. **With pre-created `tmp`:** FATAL — could not make temporary directory under `--outTmpDir`.

Example artifact root used for the successful **start** smoke (not a full overnight completion):  
`/storage/ucsf_gex_hist_validate_20260326/historical_vanilla_7a7fb08_zcat`

Planning agents should still schedule **full-sample** historical runs explicitly; this check only **unblocked** the arm.

---

## Outputs

Per run directory:

- `RUN_COMMAND.sh` — reproducible quoted command
- `BENCHMARK_SUMMARY.txt` — `mode=…`, timing, cells, `read_decompress=zcat` or `native_gzip`, `outdir=…`
- Normal STAR tree: `Solo.out/GeneFull/…`, `Log.out`, etc.

### Detecting completion (agents)

The script writes `BENCHMARK_SUMMARY.txt` **only after** STAR returns. If that
file is present with sensible fields, treat the run as **finished** unless the
outdir was reused for a failed retry without cleanup. Use
`Solo.out/GeneFull/filtered/barcodes.tsv` to corroborate. Do not rely on
`Log.out` length or `pgrep` as the primary signal (see `AGENTS.md` **Benchmark
Hygiene**).

### CR9 GEX parity (GeneFull only)

After a successful GEX-only run, compare `Solo.out/GeneFull` to Cell Ranger’s
Gene Expression MEX (no `outs/crispr_analysis` on the STAR side):

```bash
CR_RUN=/storage/ucsf-full/bench_20260218_dynamic_first/cellranger_runs/cr_full_iPSC2_1_AALG2_crstar32_20260218_205804 \
OUT_REPORT=/storage/.../star_optimized.../GEX_PARITY_vs_CR9.txt \
STAR_RUN=/storage/.../star_optimized... \
  ./scripts/run_ucsf_gexonly_gex_parity_vs_cr.sh
```

Uses `report_additional_parity_metrics.py --skip-feature-call-parity` (per-barcode
and per-gene correlations; same NXT translation defaults as the paper methodology).
Override `CR_RUN` if your CR reference is a different run that matches the GEX FASTQs.

---

## Build recipe: historical binary (`7a7fb08` worktree)

After `make clean` in `core/legacy/source`:

1. Regenerate embedded defaults before linking, e.g.  
   `xxd -i parametersDefault > parametersDefault.xxd`  
   (avoids `Depend.list` / stale-default issues; see overnight blockers handoff if checked in.)
2. `make -C core/features/libscrna` (or equivalent) if `STAR` links `libscrna.a`.
3. `make -j8 STAR`

Sanity check (optional): `strings …/STAR | rg soloAddTagsToUnsorted` should be **empty** if `parametersDefault` on disk has no such token and xxd was regenerated cleanly.

---

## Planning implications

- **Runbooks:** GEX-only overnight matrix lists **two** STAR invocations with **explicit mode flags**, not one CLI for both SHAs.
- **Policy text:** “Internal gzip only” applies to **modern**; **historical** documents **`zcat`** or accepts a documented exception.
- **Parity / CR9:** GeneFull filtered/raw comparisons across modes need careful interpretation (strand, multimappers, bridge). CR reference reuse is orthogonal.
- **Automation:** Enforce **exactly one** mode flag in wrappers and CI smoke scripts.

---

## File location

- `scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh` (executable).

Forks should preserve the two command families and document any intentional drift.

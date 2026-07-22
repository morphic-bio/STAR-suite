# Handoff: UCSF GEX-only Parity Fix (2026-03-26)

## Scope

This handoff is only about the UCSF **GEX-only** parity surface:

- STAR run: `run_ucsf_gexonly_no_bam_benchmark.sh`
- comparison target: Cell Ranger 9 Gene Expression MEX
- feature surface: `Solo.out/GeneFull`

It is **not** about the perturb-seq parity table. The perturb parity numbers in
`README.md` are still the current valid ones.

## Current Branch / Checkout

- Repo: `/mnt/pikachu/STAR-suite`
- Branch: `feature/flex-optimization-using-solo-20260325`

## What Is Known Good

### UCSF GEX-only benchmark timings

Legacy STAR:

- Run: `/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_vanilla_7a7fb08_retry2/`
- Summary: `BENCHMARK_SUMMARY.txt`
- Wall: `23.92 min`
- Cells: `13,872`

Current optimized STAR, zcat:

- Run: `/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_zcat_20260326/`
- Summary: `BENCHMARK_SUMMARY.txt`
- Wall: `13.75 min`
- Cells: `13,723`

Current optimized STAR, native gzip:

- Run: `/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_retry2/`
- Summary: `BENCHMARK_SUMMARY.txt`
- Wall: `15.78 min`
- Cells: `13,728`

Conclusion:

- current optimized + `zcat` is the best validated UCSF GEX-only timing
- current vs legacy improvement is real: about `10m10s`

### The wrapper itself is not the GEX parity problem

Script:

- `/mnt/pikachu/STAR-suite/scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh`

It intentionally has two different command families:

- `--historical-vanilla`
- `--modern-optimized`

So legacy vs current are **not** apples-to-apples flag sets, even though they
share the same wrapper file.

### The helper does use GeneFull

Parity helper:

- `/mnt/pikachu/STAR-suite/scripts/run_ucsf_gexonly_gex_parity_vs_cr.sh`

This helper validates:

- `STAR_RUN/Solo.out/GeneFull/raw`
- `STAR_RUN/Solo.out/GeneFull/filtered`

So the earlier bad parity result was **not** caused by comparing the wrong STAR
feature matrix.

## What Is Wrong

The old UCSF GEX-only parity result was against the wrong CR run.

The stale pinned CR path was:

- `/storage/ucsf-full/bench_20260218_dynamic_first/cellranger_runs/cr_full_iPSC2_1_AALG2_crstar32_20260218_205804/`

That run has:

- `233,828,574` reads
- `7,325` filtered cells

from:

- `outs/metrics_summary.csv`

This is obviously not the matching CR surface for the UCSF full GEX-only STAR
run, which is:

- `444,896,731` reads
- `13,723` to `13,872` filtered cells depending on arm

Consequences:

- the old `GEX_PARITY_vs_CR9.txt` numbers are not trustworthy
- the near-zero barcode correlations are expected from the mismatched file set
- the old README GEX-only parity note was misleading

## What Was Changed

### GEX parity helper hardened

Updated:

- `/mnt/pikachu/STAR-suite/scripts/run_ucsf_gexonly_gex_parity_vs_cr.sh`

Change:

- `CR_RUN` is now required explicitly
- the helper no longer silently defaults to the stale `cr_full_iPSC2_1_AALG2...`
  path

### README note corrected

Updated:

- `/mnt/pikachu/STAR-suite/README.md`

Change:

- the old GEX-only CR parity note was replaced with a warning that parity is
  pending a correct matching UCSF `EBs2_2` Cell Ranger path

## What Was Searched

I searched `/storage` and `/mnt/pikachu` for a matching local UCSF Cell Ranger
run with roughly the correct scale.

Result:

- I did **not** find a local CR run on this host that clearly matches the
  `EBs2_2` full GEX-only benchmark surface
- the only obvious local UCSF CR run under `cellranger_runs/` is the stale one
  above

## Current Valid Numbers

### Valid

- UCSF GEX-only timing:
  - legacy `23.92 min`
  - current optimized + `zcat` `13.75 min`
- perturb timings:
  - UCSF `16.4 min`
  - MSK `25.0 min`
- perturb parity table in `README.md`

### Not valid yet

- UCSF GEX-only barcode/gene correlation vs CR9

## Next Agent: Exact Task

1. Find the actual matching UCSF `EBs2_2` Cell Ranger 9 run.
2. Verify it has roughly the expected scale:
   - around `445M` reads
   - around `13.7k` filtered cells
3. Rerun:

   ```bash
   CR_RUN=/path/to/matching/cr_run \
   STAR_RUN=/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_zcat_20260326 \
   /mnt/pikachu/STAR-suite/scripts/run_ucsf_gexonly_gex_parity_vs_cr.sh
   ```

4. If needed, also run against the legacy GEX-only STAR arm:

   ```bash
   CR_RUN=/path/to/matching/cr_run \
   STAR_RUN=/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_vanilla_7a7fb08_retry2 \
   /mnt/pikachu/STAR-suite/scripts/run_ucsf_gexonly_gex_parity_vs_cr.sh
   ```

5. Only after that, update `README.md` and `tests/ARTIFACTS.md` with GEX-only
   parity numbers.

## Files Most Relevant

- `/mnt/pikachu/STAR-suite/scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh`
- `/mnt/pikachu/STAR-suite/scripts/run_ucsf_gexonly_gex_parity_vs_cr.sh`
- `/mnt/pikachu/STAR-suite/scripts/report_additional_parity_metrics.py`
- `/mnt/pikachu/STAR-suite/README.md`
- `/mnt/pikachu/STAR-suite/tests/ARTIFACTS.md`

## One-Line Summary

The UCSF GEX-only timing numbers are good; the old GEX-only CR parity numbers
are bad because the helper was pointed at the wrong Cell Ranger run.

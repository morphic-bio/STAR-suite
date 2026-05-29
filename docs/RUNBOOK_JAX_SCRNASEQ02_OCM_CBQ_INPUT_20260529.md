# JAX scRNAseq02 OCM CBQ Input Runbook

Date: 2026-05-29

Canonical recipe lives in morphic-recipes:

```text
/mnt/pikachu/morphic-recipes/docs/RUNBOOK_JAX_SCRNASEQ02_OCM_CBQ_INPUT_20260529.md
```

STAR-suite keeps this copy because the MCP workflow schema and compatibility
launcher live here.

## Summary

The JAX scRNAseq02 OCM composite smoke recipe now supports:

```bash
--star-input-format fastq|cbq
```

CBQ mode stages one ordered paired CBQ per lane from the STAR-specific staged
FASTQs, using STARsolo mate order:

```bash
cbq_ordered_encoder --readFilesIn <R2.fastq.gz> <R1.fastq.gz> --outFile <lane>_R2_R1.cbq
```

The rendered STAR command uses:

```bash
--readFilesType Binseq PE
--readFilesIn lane1_R2_R1.cbq,lane2_R2_R1.cbq
```

CBQ mode must use:

```bash
--star-yremove no
```

The recipe exits if CBQ is combined with Y/noY FASTQ emission, which is still
FASTQ-only.

## Passing Smoke

```text
/tmp/ocm_cbq_smoke_20260529T174218Z
```

Details:

- clean STAR-suite worktree: `/tmp/star_suite_ocm_cbq`
- read pairs: `1000`
- threads: `8`
- STAR-suite version: `1.0.3`
- `STAR_COMPLETED.txt`: `completed_utc=2026-05-29T17:45:53Z`
- `Log.final.out`: `Number of input reads = 1000`, uniquely mapped reads = `964`

## FASTQ-vs-CBQ Parity Smoke

On 2026-05-29, a 1K FASTQ-vs-CBQ smoke using the same staged downsample passed:

- CBQ run: `/tmp/ocm_cbq_impl_prepare_20260529T181158Z`
- FASTQ run: `/tmp/ocm_cbq_impl_fastq_parity_20260529T181518Z`
- shared settings: `--read-pairs 1000`, `--threads 8`,
  `--star-yremove no`, `--star-out-samtype None`
- both runs: `Number of input reads = 1000`, uniquely mapped reads = `964`
- `star_composite/run/Solo.out`: byte-identical
- `star_composite/outs`: byte-identical, including native per-sample OCM MEX

The CBQ/Y-removal guard was also checked: `--star-input-format cbq
--star-yremove yes` exits before staging with the expected FASTQ-only error.

STAR-suite regression wrapper:

```bash
tests/run_cbq_ocm_composite_smoke.sh
```

The wrapper is registered in `tests/production_module_regression_manifest.tsv`
as a host-local production smoke because it requires the JAX scRNAseq02 OCM raw
files and genome index. A wrapper sanity run passed at:

```text
/tmp/star_suite_cbq_ocm_wrapper_check_20260529T182012Z
```

## 100K Parity Run

The release-level 100K FASTQ-vs-CBQ OCM run passed on 2026-05-29:

```text
/mnt/pikachu/JAX_scRNAseq02_processed/ocm_cbq_100k_20260529T182838Z
```

Settings:

- `READ_PAIRS=100000`
- `THREADS=16`
- `RUN_FASTQ_PARITY=1`
- CBQ and FASTQ both used `--star-yremove no` and `--star-out-samtype None`

STAR summaries:

| Mode | Started | Mapping Start | Finished | Input Reads | Unique Reads | Unique % |
| --- | --- | --- | --- | ---: | ---: | ---: |
| CBQ | May 29 18:28:40 | May 29 18:29:42 | May 29 18:31:03 | 100000 | 95687 | 95.69 |
| FASTQ | May 29 18:31:03 | May 29 18:32:09 | May 29 18:33:29 | 100000 | 95687 | 95.69 |

Parity:

- `fastq/star_composite/run/Solo.out` vs `cbq/star_composite/run/Solo.out`:
  byte-identical by `diff -qr`
- `fastq/star_composite/outs` vs `cbq/star_composite/outs`:
  byte-identical by `diff -qr`
- Native OCM filtered cells per tag matched:
  `OB1=5325`, `OB2=4127`, `OB3=4277`, `OB4=3757`
- CBQ/Y-removal gate still failed as expected before staging.

Input size check:

- staged FASTQ.gz bytes: `15658160`
- staged CBQ bytes: `10461811`
- output root size: about `504M`

## Reproduce

```bash
git -C /mnt/pikachu/STAR-suite worktree add --detach /tmp/star_suite_ocm_cbq origin/master
make -C /tmp/star_suite_ocm_cbq/core/legacy/source clean
make -C /tmp/star_suite_ocm_cbq/core/legacy/source -j8 STAR cbq-ordered-encoder

OUT_ROOT=/tmp/ocm_cbq_smoke_$(date -u +%Y%m%dT%H%M%SZ)

STAR_SUITE_ROOT=/tmp/star_suite_ocm_cbq \
/mnt/pikachu/morphic-recipes/scripts/run_jax_scrnaseq02_ocm_composite_smoke.sh \
  --read-pairs 1000 \
  --out-root "$OUT_ROOT" \
  --prepare \
  --run-star \
  --force \
  --threads 8 \
  --star-input-format cbq \
  --star-yremove no \
  --star-out-samtype None
```

Check:

```bash
test -f "$OUT_ROOT/STAR_COMPLETED.txt"
test -s "$OUT_ROOT/stage/star_composite_cbq/cbq_manifest.tsv"
grep -F -- '--readFilesType Binseq PE' "$OUT_ROOT/RUN_STAR_COMPOSITE.sh"
! grep -F -- '--readFilesCommand' "$OUT_ROOT/RUN_STAR_COMPOSITE.sh"
! grep -F -- '--emitYNoYFastq' "$OUT_ROOT/RUN_STAR_COMPOSITE.sh"
```

## Release-Level Follow-Up

Run a 100K FASTQ-vs-CBQ comparison with `--star-yremove no` on both sides and
compare raw `Solo.out/GeneFull` plus `Solo.out/Velocyto` outputs after any
needed canonicalization. Record wall time, output paths, and comparison status
before treating OCM CBQ as a release gate.

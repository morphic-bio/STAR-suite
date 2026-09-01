# BGZF ingest benchmarks — 2026-09-01

Protocol: cold page cache, 32 threads, full JAX SC2300771 eight-lane no-align
Flex count-only invocation, comparing the existing internal gzip path, BGZF
record-range readers, and indexed ordered CBQ range readers. Each timed run uses
a fresh external output directory and `/usr/bin/time -v`; command lines and
logs are retained in this directory.

The JAX FASTQs are structurally complete BGZF streams without the optional
canonical empty EOF member. Index integrity is established by `BC/BSIZE` hops
ending exactly at physical EOF, followed by raw inflate, ISIZE, and CRC checks.

The public 320k scFFPE dataset was skipped because
`/mnt/pikachu/tenx_320k_scFFPE` contains the download tar and `config.csv`, but
the FASTQs are not yet extracted.

Run order:

```bash
make -C core/legacy/source bgzf-input-harness
docs/benchmarks/bgzf_ingest_20260901/prepare_jax_bgzi.sh
docs/benchmarks/bgzf_ingest_20260901/run_jax_benchmark.sh gzip
docs/benchmarks/bgzf_ingest_20260901/run_jax_benchmark.sh bgzf
docs/benchmarks/bgzf_ingest_20260901/run_jax_benchmark.sh cbq
```

Timed runs must be serialized and launched only when the machine is otherwise
quiet. The BGZF index build is deliberately outside the timed region because
the `.bgzi` file is an amortized cache, matching the runbook design.

## Results

Pending quiet-host execution.

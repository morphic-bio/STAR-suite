# BGZF ingest benchmarks — 2026-09-01

Protocol: cold page cache, 32 threads, full JAX SC2300771 eight-lane no-align
Flex count-only invocation, comparing the existing internal gzip path, BGZF
on-demand range readers, and indexed ordered CBQ range readers. FASTQ and CBQ
inputs were staged on the `/home` NVMe SSD for their respective rows. Each timed
run uses a fresh external output directory and `/usr/bin/time -v`; command lines
and logs are retained in this directory. The wrapper requires the
`Flex count-only no-genome: active` marker before writing its success marker.

The JAX FASTQs are structurally complete BGZF streams without the optional
canonical empty EOF member. Stream integrity is established by `BC/BSIZE` hops
ending exactly at physical EOF, followed by raw inflate, ISIZE, and CRC checks.

The public 320k scFFPE dataset was skipped because its tar was still incomplete
(`361,491,693,568` of `533,915,095,040` bytes) and the FASTQs were not extracted.
The download was stopped throughout all timed runs.

Run order:

```bash
docs/benchmarks/bgzf_ingest_20260901/run_jax_benchmark.sh gzip
docs/benchmarks/bgzf_ingest_20260901/run_jax_benchmark.sh bgzf
docs/benchmarks/bgzf_ingest_20260901/run_jax_benchmark.sh cbq
```

Timed runs must be serialized and launched only when the machine is otherwise
quiet. BGZF readers discover and cache bounded work ranges from inline
`BC/BSIZE` fields during the timed run; no pre-index or sidecar is used.

## Results

All three successful rows produced exactly the same Flex counters:
`2,011,130,186` total, `1,681,459,858` keep, `16,111,757` deny, and
`313,558,571` miss.

| Input path | Commit | Wall | CPU | Peak RSS | Versus CBQ |
| --- | --- | ---: | ---: | ---: | ---: |
| gzip (`readFilesBgzfMode off`) | `eb564c5` | 10:42.19 | 718% | 44,437,864 KiB (42.38 GiB) | +41.48% |
| BGZF on-demand range | `b639121` | 7:38.71 | 1583% | 53,259,228 KiB (50.79 GiB) | +1.06% |
| CBQ indexed range | `eb564c5` | 7:33.91 | 1500% | 48,396,764 KiB (46.15 GiB) | baseline |

The gzip and CBQ controls were completed immediately before the two BGZF-only
performance commits (`719e863`, `b639121`). Neither control path invokes the
changed BGZF code, so the completed cold-cache controls were not duplicated.
The accepted log prefixes are:

- `jax_gzip_20260901T210716Z`
- `jax_bgzf_20260901T221610Z`
- `jax_cbq_20260901T213904Z`

An initial 13-minute gzip run was rejected because a nonzero chimeric setting
disabled the cache-only path. Two pre-acceptance BGZF profiling rows were also
excluded. Their logs and external outputs are retained outside the committed
results for diagnosis.

# BGZF reader optimization, 2026-09-04

All runs used the full JAX BGZF FASTQ dataset staged on the `/home` NVMe SSD,
32 pipeline threads, `--flexNoAlign 1`, RAM buckets, CRC checking, a cold page
cache, and a fresh output directory. Each run has its command, preflight,
completion marker, `/usr/bin/time` output, pipeline logs, output manifest, and
continuous host-activity trace in the correspondingly named directory.

| Commit | Candidate state | Wall (s) | CPU | Peak RSS (GiB) | Decision |
|---|---|---:|---:|---:|---|
| `ad652fc` | Initial persistent-inflater/batched-handoff candidate | 298.36 | 2430% | 40.6 | Refine work granularity |
| `d90cf5a` | Fine-grained producer work restored | 269.43 | 2381% | 40.8 | Keep |
| `e3d4ad3` | Compressed-size-weighted R2/R1 inflaters | 268.65 | 2406% | 40.6 | Keep; neutral alone |
| `e5cd3e4` | Direct record handoff | 258.31 | 2457% | 40.5 | Keep |
| `98b6a58` | Fixed-size FASTQ records | 248.45 | 2424% | 40.6 | Keep |
| `1b0c2d8` | Skip unused quality copy plus handoff cleanup | 248.37 | 2414% | 40.4 | Keep; neutral wall, lower CPU |
| `810d26b` | Recycle inflated buffers across line boundaries | 248.65 | 2406% | 40.7 | Revert; no wall benefit |
| `1f7387a` | Bounded libc mate-name delimiter scans | 238.69 | 2468% | 40.5 | Keep |
| `ec47f15` | Special-case unused R2 name fields | 238.55 | 2474% | 40.6 | Revert; within noise |

The retained path reduced wall time from 269.43 to 238.69 seconds (30.74
seconds, 11.4%). The final delimiter-scan change accounts for 9.68 seconds of
that reduction relative to its immediate 248.37-second baseline. Every valid
run reported 2,011,130,186 input pairs, the same Flex pipeline counters, 16
sample outputs, and manifest digest `767c1e314df547e3`. No competing compute
job appeared in the accepted run's activity trace.

The final branch includes explicit revert commits for the two neutral
experiments. Consequently, its reader source is identical to the measured
`1f7387a` state.

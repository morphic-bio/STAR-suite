# FASTQ and CBQ reader optimization, 2026-09-04

All runs used the full JAX dataset staged on the `/home` NVMe SSD, 32 pipeline
threads, `--flexNoAlign 1`, RAM buckets, a cold page cache, and a fresh output
directory. Each artifact directory contains the command, preflight, completion
marker, `/usr/bin/time` output, pipeline logs, output manifest, and continuous
host-activity trace.

## Plain gzip FASTQ

The retained fast path consumes unused FASTQ lines without materializing
`std::string` objects, scans the read-name token within the fixed input buffer,
and consumes no-align R2 quality into scratch storage. It continues to use the
existing serial `gzFile` stream.

| State | Wall (s) | CPU | Peak RSS (GiB) | Manifest |
|---|---:|---:|---:|---|
| Control | 497.82 | 848% | 39.4 | `767c1e314df547e3` |
| Fast path, run 1 | 487.35 | 842% | 39.5 | `767c1e314df547e3` |
| Fast path, run 2 | 487.91 | 842% | 39.6 | `767c1e314df547e3` |

The optimized median was 487.63 seconds, 10.19 seconds (2.0%) faster than the
control. Both optimized runs reproduced all pipeline counters, 16 sample
outputs, and the control manifest digest. The fast path is retained.

## CBQ

Three candidates were measured. The initial handoff change removed unused
read-name and quality copies by passing borrowed fixed-buffer spans. A second
candidate streamlined packed-base decoding. A third divided the global record
space into four times as many claimable ranges (135 lane fragments instead of
39 on JAX) to test for end-of-range imbalance.

| Candidate | Wall (s) | CPU | Peak RSS (GiB) | Decision |
|---|---:|---:|---:|---|
| Borrowed fixed-buffer handoff | 218.54 | 2227% | 39.8 | Neutral |
| Streamlined packed-base decode | 217.66 | 2202% | 40.1 | Keep for simpler decoding |
| Fourfold range over-partitioning | 218.54 | 2249% | 39.8 | Revert; neutral |

The established CBQ result was 218.23 seconds, so none of the candidates
demonstrated a material reader-wall improvement. The decode candidate was only
0.57 seconds (0.3%) faster, while its mapping-phase wall time was unchanged.
It is nevertheless retained because separating the unaligned prefix, complete
packed bytes, and trailing bases eliminates repeated offset arithmetic and
makes the bounds behavior explicit.
The range candidate exactly reproduced the 218.54-second handoff result.
All three runs reproduced the pipeline counters, 16 sample outputs, and
manifest digest `767c1e314df547e3`.

CBQ already reads sequences into bounded buffers and exposes the other fields
as borrowed views. The measurements therefore indicate that additional
buffering and range scheduling are not limiting this workload. The range
scheduling candidate remains explicitly reverted.

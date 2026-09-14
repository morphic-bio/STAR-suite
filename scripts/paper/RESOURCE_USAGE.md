# Benchmark resource sampler

`resource_usage.py` records concurrent workload memory and scratch filesystem
usage for benchmark wrappers. It is the unchanged sampler used for the
September 2026 H1 HDD and M6id SSD Cell Ranger, STAR 1.9.4, and cyto campaigns.

The archived campaign copy has SHA256
`af90e16999541d3cedf11a3e8f8191aba558a564a7ebcca5de0c596c76195e09`.
Its original location in the companion paper repository is
`analysis/paper_benchmark_refresh_20260910/i3_retirement_m6id_replacement_20260913/package/resource_usage.py`.
The first H1 GSE325982 Cell Ranger run preceded this sampler.

## Use from a Python wrapper

This is an importable module, not a standalone command-line program. It needs
Python 3 on Linux, readable `/proc` and pressure files, and permission to read
the workload's process status. It uses only the Python standard library.

From the STAR Suite repository root:

```python
from pathlib import Path
from scripts.paper.resource_usage import run_measured

run_dir = Path("/scratch/benchmarks/example_run")
run_dir.mkdir(parents=True, exist_ok=False)

with (run_dir / "stdout.log").open("w") as out, \
     (run_dir / "stderr.log").open("w") as err:
    result = run_measured(
        ["/usr/bin/time", "-v", "-o", str(run_dir / "time.txt"),
         "/absolute/path/to/STAR", "--version"],
        run_dir,
        scratch="/scratch",
        interval=5,
        cwd=run_dir,
        stdout=out,
        stderr=err,
    )
raise SystemExit(result.returncode)
```

Replace the example command with the full benchmark command. Keep GNU time's
output for its separate maximum-RSS and CPU-time measurements. A very short
command may finish before its memory is sampled.

The wrapper must create the output directory, use a fresh directory for each
execution, keep workloads in the foreground, and serialize benchmarks on the
measured filesystem. The sampler refuses existing resource files and passes
remaining keyword arguments to `subprocess.Popen`. It returns the command's
exit code in a `CompletedProcess`; it does not raise for a nonzero exit.

## Recorded files and interpretation

- `RESOURCE_USAGE.jsonl`: samples every five seconds by default, plus a final
  sample. Records process count, process-tree RSS and swap, host memory,
  scratch filesystem space, selected device I/O counters, CPU counters,
  load average, and CPU pressure.
- `RESOURCE_SUMMARY.json`: sampled peaks, sampling errors, before/after disk
  allocation, peak disk growth, timestamps, and workload exit code. Inspect
  `errors` and `samples` before accepting a measurement.

Memory and disk values are bytes. Divide by `2**30` for GiB.

| Measurement | Meaning |
|---|---|
| `sampled_peak.process_tree_rss_bytes` | Maximum sampled sum of RSS across concurrently running workload processes. The traversal checks children of every thread. Shared resident pages can be counted in more than one process, and short peaks between samples can be missed. |
| GNU time maximum RSS | A separate Linux per-process high-water mark that does not sum concurrent workers. Do not substitute it for process-tree RSS in comparisons. |
| `sampled_peak.cgroup_memory.peak` | When available, the cgroup memory peak, including charged filesystem page cache and kernel memory. It is not interchangeable with process RSS. |
| `sampled_peak.host_used_excluding_available_bytes` | Host-wide `MemTotal - MemAvailable`, including unrelated processes. |
| `sampled_peak.disk_used_bytes` | Total allocation on the filesystem containing `scratch`, including staged inputs, references, and earlier results. |
| `disk_peak_growth_bytes` | Peak sampled filesystem allocation minus the pre-run baseline. Includes new temporary files and outputs, excludes pre-existing allocation, and can be affected by unrelated writers or cleanup. This is not cumulative bytes written. |
| `disk_used_after_bytes - disk_used_before_bytes` | Net retained filesystem growth at exit, not a recursive measurement of the output directory. |

The sampler reads filesystem allocation counters without recursively scanning
output directories during a run. `diskstats` contains raw cumulative Linux
device counters for `md0`, `nvme*`, and `xvd*` devices. Compare start/end values
for one chosen device layer; do not add RAID and member-device counters.

## systemd scope support

For commands starting with literal `systemd-run` and an explicit
`--unit=NAME`, the sampler also reads processes and memory counters from
`/sys/fs/cgroup/system.slice/NAME.scope`. This covers the benchmark campaigns'
`systemd-run --scope` launch convention even when process ancestry changes.
User scopes, other slices, wrapped or absolute-path `systemd-run` invocations,
and arbitrary detached daemons are not automatically discovered by this version.

The sampler does not impose a memory limit or manage benchmark provenance by
itself. Save the exact command, binary identity, inputs, environment, and
completion checks in the calling wrapper. Preserve the deployed sampler with
each campaign when making future changes to this tracked copy.

# Per-shard trim QC for the sharded bulk pipeline

Date: 2026-08-07

## Why

The sharded bulk pipeline maps on N nodes and never assembles the reads
anywhere, so there is no single point at which QC can see all of them. The two
existing options both cost something:

- in-process `--trimQcReport` gives one report *per shard*, not one for the run;
- `trim_qc_replay` gives a merged report but needs an alignment file, and it
  measures aligned reads rather than everything that passed through.

## What was added

`TrimQcCollector` is a purely additive set of histograms — length, per-position
count, per-position quality sum, per-position A/C/G/T/N, and a 101-bin GC
histogram — and STAR already merges one per chunk thread. That makes it
shardable for free. Nothing needs to be recomputed and no read is revisited;
only the counters move.

- `slam/source/TrimQcShard.{h,cpp}` — binary dump/load of the collector plus the
  eight trim counters `writeTrimQcJson` needs off `Stats`. All are plain sums.
- `TrimQcCollector::restoreFromShard` — the one entry point added to the class so
  a dump can be rehydrated without loosening its private state.
- `--trimQcShardOut PATH` — STAR writes the accumulator next to its evidence
  sidecar. It also enables trim-QC collection on its own, so a sharded run does
  not have to ask for a per-shard JSON it will throw away.
- `trim_qc_merge` — merges N dumps and writes the same `trim_qc.json` /
  `trim_qc.html` the in-process path emits.
- `trim_qc_fastq --shard-out` — the same dump from the FASTQ tool, which makes
  the merge testable without a genome.

## Cost

The dump is the histograms, not the reads. On the smoke fixture the four shard
dumps are 4.7–10.5 KB each, 28 KB total; size scales with read length, not read
count. For 150 bp paired reads it is on the order of 20 KB per shard.

This is the same move as the VB/EM gather and the 101-bin GC rendezvous: send
the summary, not the data.

## The producer was missing, and is now wired

`--trimQcReport` emitted a well-formed report with `"total_reads": 0`.
`TrimQcCollector::addRead` was never called anywhere in STAR: `ReadAlign`
constructed a collector and `STAR.cpp` merged one per chunk thread, but nothing
between them fed it a read. The files were produced, so the failure was
invisible unless you opened them. Confirmed empirically -- 2,000 PE reads
against a 200 kb reference produced a report with every histogram empty -- and
`trimQc.addRead` has never appeared anywhere in this repository's history.

The producer now lives in `ReadAlign_oneRead.cpp`, in the same
`outFilterBySJoutStage != 2` block as the existing `qualHist` accumulation.
That location is forced: a few lines later mate 2 is reverse-complemented into
`Read1[0]` and the per-mate view is gone, and sharing the two-pass guard with
`qualHist` is what keeps a two-pass run from counting every read twice.
`Read1[im]` is clipped numeric sequence while `Qual0[im]` stays in original
coordinates, so the left clip is passed as the quality offset rather than
applied to the pointer.

## Correctness

`tests/run_trim_qc_merge_smoke.sh`:

- 4 shards with **different read lengths** (75/150/100/60, 4,700 reads total),
  so the merge has to grow its per-position arrays — the case a naive
  accumulator gets wrong.
- Merged output is **byte-identical** to a single pass over the concatenated
  reads.
- Merging in a different shard order is byte-identical too.

## Y removal is not affected

`emitNoYBAM` / `emitYNoY` classify reads by whether they have a Y-chromosome
alignment, so that decision needs the alignment and stays in the mapping phase.
Per-read Y artifacts remain restricted under the data-use agreement recorded in
`docs/HANDOFF_PE_BULK_YREMOVE_TRANSCRIPTVB_REGRESSION_20260310.md`; `emitYReadNames`
is the cheap path when only a name list is needed for downstream filtering.

## Usage

```bash
# per shard, during mapping
STAR ... --trimQcShardOut "$shard_dir/trim_qc.trimqc"

# once, at gather
trim_qc_merge --out-prefix "$out/run" "$run"/shard.*/trim_qc.trimqc
```

# STAR Suite v1.9.3 Release Notes

Date: 2026-09-11

STAR Suite 1.9.3 repairs feature barcode assignment at scale and removes the
storage patterns that made large runs slow to finish. Since 1.7, position
learning was switched off whenever a feature offset was given explicitly, which
is what the guide and lineage workflows do. Without learning, every read
searched the whole feature library: A375 lost 210,785 assigned feature reads,
and the MSK lineage library (245,979 LARRY barcodes) did not finish in three
hours. 1.9.3 restores learning for explicit offsets, removes repeated work in
the learning and fallback searches, and replaces the node-based barcode, UMI and
read-counter storage found by an allocation audit. It includes the merged
changes through `0143c16`.

`STAR --version` reports `1.9.3`. Debian source packaging uses `1.9.3-1`;
Ubuntu packages use `1.9.3-1~ubuntu22.04.1` and
`1.9.3-1~ubuntu24.04.1`. Upstream STAR remains `2.7.11b`, genome-index
compatibility remains `2.7.4a`, and legacy compatibility remains `2.7.1a`.
Existing indexes do not need rebuilding.

## Feature barcode assignment

- **Explicit feature offsets learn positions again.** A run that passes a
  feature offset (including offset 0, which the Perturb-seq and lineage
  wrappers use) now learns each feature's match position from the first
  100,000 reads, as it did before 1.7. A375 assigns 9,274,336 feature reads,
  matching its March result, rather than 9,063,551.
- **Less repeated work in large libraries.** Learning over a single uniform
  anchor group now visits only the features a cached offset lookup can return,
  and the broad fallback scan is skipped for reads that a complete
  Hamming prehash proves cannot match. Feature and offset ordering, distance
  ceilings, ambiguity rules and learned histograms are unchanged.
- **Expansion buffers with no permitted Ns.** Buffer sizes used a negative
  shift exponent and could be zero, and the fallback could modify the input
  sequence while storing a match position. Both are fixed.
- **Feature learning is deterministic.** Serialized learning no longer depends
  on which input lane happens to be ready first.

## Memory and run completion

An audit of container and allocation patterns
(`docs/AUDIT_ALLOCATION_PATTERNS_20260911.md`) found several places where the
suite allocated one small container per barcode, UMI or variant. Those are
replaced, with results recorded in `docs/ALLOCATION_FIXES_20260911.md`:

- Cell-barcode correction keeps packed lookup tables and one candidate array
  instead of maps of per-variant vectors, and no longer copies millions of
  whitelist strings; GEX gather borrows Solo's existing packed whitelist table.
- Single-feature barcode–UMI counters live inline in their pooled record, so a
  general hash is allocated only for real multi-feature UMIs (A375: 99.968%
  inline).
- Bridge read counters, inline barcode collision storage, occupancy grouping,
  file-based Flex matrices, feature merge and table import, and Velocyto UMI
  counters use flat or borrowed storage.

Together these remove the long single-threaded cleanup that ran after a run
reported "ALL DONE": in the MSK test the process exited about a second after
STAR finished.

## Velocyto without BAM

Velocyto requested per-read information that the direct hash bridge did not
keep, so `--outSAMtype None` runs produced empty layers. The bridge now retains
read identities for the selected gene-like feature, using the same resolved
ambiguous-barcode decisions as the count hash and replaying the final corrected
or rejected UMI status. GEX-only runs allocate no read buffers. Velocyto now
reports an error when its source lacks per-read storage instead of writing
empty layers.

## Scheduling (opt-in, off by default)

`--dynamicThreadBgzfHierarchy 1` admits BGZF decoders within their workload's
permit budget, and `--dynamicThreadBalance 1` additionally balances the
finish times of alignment and feature assignment. Both default to `0`, and
existing runs are unaffected. Feature assignment also reports its phases
(read and assign, hash merge, barcode rescue, UMI deduplication, filtering,
output, cleanup).

## Validation

Full MSK ES workload (GEX + PolyIII guides + LARRY lineage, 32 threads, no BAM),
against the accepted April 2026 run:

| Measurement | 1.9.3 candidate | April 2026 |
| --- | ---: | ---: |
| Whole workflow | 1,706 s | 1,811 s |
| LARRY assignment (245,979 features) | 1,440.3 s | 1,442.9 s |
| PolyIII assignment (30 guides) | 56.4 s | 146.4 s |
| LARRY raw matrix | identical (3,667,851 entries) | — |
| GEX barcode Jaccard vs Cell Ranger 9 | 0.9922 | 0.9822 |
| Guide calls matching Cell Ranger 9 | 98.9708% | 98.9680% |

On 1.9.2 the same workload was still in LARRY assignment after three hours.
PolyIII differs from April in 14 of 2,168,685 matrix entries, each by one UMI,
from three additional assigned reads; this follows the deterministic-learning
change. Cell calls (32,898) reflect the 1.9.2 floor separation, not this
release.

## Known issue

Guide calls are written only for OrdMag primary cells, so cells rescued by
EmptyDrops have no row in `protospacer_calls_per_cell.csv` (MSK: 2,585 cells).
This is unchanged from 1.9.2 and is tracked in
`docs/HANDOFF_CRISPR_CALLS_MISSING_TAIL_CELLS_20260911.md`.

## Compatibility

Counts change only where assignment was previously losing reads: guide and
lineage libraries gain the reads that explicit-offset runs had been missing.
GEX counting, cell calling and Flex outputs are unchanged from 1.9.2. Hosted
tarball and Debian packages retain the portable no-Chromap build; local
production builds retain the Chromap-enabled default. The paper refresh pins
new executions to 1.9.3.

# Ovarian GEX finalization audit, 2026-09-15

The observed v1.7.1-to-v1.9.5 slowdown is accounted for by increased work in
the two stages that execute the additional assignment policies. No changed
upstream accounting, unexplained spool duplication, or changed strict/hard
output was found. A separate, real robustness defect was reproduced in the
soft-policy UMI ordering. It is present in both releases and in both the
memory and spill implementations. Its impact on the ovarian data is unknown.

The user requested investigation of the slowdown as a possible correctness
canary. This audit includes a clean v1.9.5 STAR build, source comparisons,
retained full-run metadata, and small standalone resolver diagnostics. It does
not change the released binary or scientific results.

## Correction to the earlier timing interpretation

The first timing review consulted progress logs and matrix timestamps and said
that existing evidence did not separate resolver and materializer timing.
That was incomplete: `starSpatialGex.out/run_summary.tsv` already records both.
This audit supersedes that part of the earlier interpretation. The earlier
provenance is retained without rewriting its history.

| Native stage | v1.7.1, strict/hard | v1.9.5, four policies | Change |
| --- | ---: | ---: | ---: |
| Spill merge and clique/contribution construction | 138.62 s | 143.39 s | +4.77 s |
| Downstream molecule resolution | 70.66 s | 305.52 s | +234.86 s |
| Matrix construction, writing and commit | 150.31 s | 358.22 s | +207.91 s |
| All remaining elapsed time | 1,852.47 s | 1,804.46 s | -48.01 s |
| Total | 2,212.06 s | 2,611.59 s | +399.53 s |

Resolution timing includes loading the contribution shards; materialization
is the downstream interval minus that resolution timer, including setup,
matrix-run creation, text MEX output, commit and some cleanup. These are the
actual code timer scopes, not pure compute or I/O measurements.

The older v1.6.1 full ovarian primary also requested all four policies. It
spent 317.79 s resolving and 376.86 s materializing, versus 305.52 and
358.22 s now. Its total was 47:13.65. This is useful historical evidence,
not a matched-environment timing control.

## Correctness and work accounting

All three full primaries agree exactly on 22 native fields, including:

- 474,131,092 decoded reads;
- 373,655,322 joined reads and 467,830,993 candidate rows;
- 205,906,769 read cliques;
- 257,318,254 downstream contribution records in 256 shards;
- 164,293,957 strict and 201,888,275 hard molecules.

The sealed manifests match for all 18 v1.7.1 strict/hard MEX components and
all 36 v1.6.1 four-policy components. These are comparisons of the existing
sealed checksums, not a new scan of every large matrix.

The doubled matrix-run count is exactly the expected
`256 shards * 3 scales * policies`: 1,536 runs for two policies, 3,072 for four.
The recorded matrix-run bytes agree exactly with the sum of exported nonzeros
times the 16-byte record size plus each run's 72-byte header and 32-byte
trailer. Contribution bytes likewise equal 32 bytes per contribution plus
the same run overhead. There is no unexplained record or byte excess.

The decoder, capacity code, spill implementation, downstream spool and
multi-gene reconciliation sources are byte-identical between v1.7.1 and
v1.9.5. The contiguous helper block containing clique construction,
integer/soft resolution and materialization in `SpatialGexPipeline.cpp` is
also byte-identical. The changes elsewhere in that file are Flex evidence
accounting and its validation/output fields. GEX initialization/finalization
was moved into lambdas in `STAR.cpp`; the GEX finalization call remains once
at the same point after mapping. `Pipeline::finalize` rejects double calls.

A fixed diagnostic of 100,000 synthetic cliques / 249,604 contributions
compared every policy alone and strict+hard together against all policies
enabled. Every retained molecule identity and weight matched exactly, and
the input contribution records were unchanged. This passed on both releases.

## Confirmed soft-policy defect

Both `buildSoftRawSupport` and `resolveSoftContributions` sort UMI supports
using descending support, except that nearly equal supports are ordered by
UMI. Approximate equality is not transitive. With UMIs 0, 1, 2 and supports
0.5, 0.50000000004, 0.50000000008, the comparator says:

`0 before 1`, `1 before 2`, and `2 before 0`.

This is an invalid ordering for `std::sort`. The diagnostic exercises the
production functions rather than a replacement resolver. Each test clique
has two complementary candidate probabilities summing to one and consistent
hard-policy flags. A second gene creates a controlled reconciliation contest.
Additional UMIs are each at least two substitutions from every original UMI
and from one another, so they cannot change the original correction group
through Hamming-1 edges. The original contributions and probabilities stay
fixed.

| Disconnected additional UMIs | Original group's corrected UMI | Original group's gene-0 mass |
| --- | ---: | ---: |
| 0 | 2 | 0.87500000003 |
| 13 | 2 | 0.87500000003 |
| 14 | 0 | 0 |
| 29 | 2 | 0.87500000003 |

When the representative changes to UMI 0, the other gene wins reconciliation
at that UMI. Thus this can change gene counts, not merely an internal label.
The 14-extra-UMI case takes the group above the library's small-sort boundary;
the invalid comparator allows a different representative. It is reproduced
identically on v1.7.1 and v1.9.5. Both memory and spill output agree on the
defective result.

An external diagnostic source copy replaced only the two approximate-equality
sort comparisons with exact support comparisons, retaining the UMI secondary
key. The original group then stayed at UMI 2 and mass 0.87500000003 for every
0-through-64-extra-UMI case. The policy-isolation checks also passed on that
copy. This isolates the cause; it is not a released or integrated fix.

The first comparator entered in `1d41bea` (2026-07-24), and the spill copy in
`e5f24f0` (2026-07-25). The affected scope is soft_expected spatial resolution,
shared by spatial GEX and Flex. Integer strict/hard/gated correction uses
integer support ordering and does not contain this comparator defect.

There is no evidence yet that these adversarial near ties occur in the
ovarian contribution stream, or that they caused its runtime increase.
Identical historical/current matrices cannot exclude a defect shared by both.
The next focused correction should give both soft resolvers a valid,
deterministic ordering, add this nonlocal-effect case as a regression test,
and measure any soft-output changes separately from the integer policies.

## Reproduction and retained evidence

- Core diagnostic: `tests/diagnostics/spatial_gex_policy_audit.cpp`.
- Metadata summarizer: `tests/diagnostics/summarize_spatial_gex_policy_audit.py`.
- Provenance root:
  `/mnt/pikachu/visium-hd-processing-provenance/runs/visium_hd_publication/20260915_ovarian_gex_finalization_audit`.
- Exact compile/run commands and source hashes: `commands/build_v171.json`,
  `build_v195_valid.json`, `build_v195_strict_sort_control.json` and corresponding
  `run_*.json` files. The diagnostic includes the production pipeline source
  in one translation unit; do not additionally link `SpatialGexPipeline.o`.
- Full 65-case outputs: `logs/policy_audit_{v171,v195_valid,v195_strict_sort_control}.tsv`.
- External-copy-only control patch:
  `environment/strict_sort_diagnostic_control.patch`.
- Source identity/diff record: `results/source_audit.json`.
- Generated accounting and diagnostic summary: `results/audit_summary.json`;
  generator invocation and hash: `commands/summarize01.json`.

No vendor source, full-slide rerun, alignment or reference-index load is part
of the diagnostic. The first diagnostic attempt is retained too; the final
variants additionally validate consistent hard flags and memory/spill parity.

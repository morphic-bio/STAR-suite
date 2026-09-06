# Handoff: STAR-Flex vs Cell Ranger concordance fixes (2026-09-05/06)

Branch `fix/sample-tag-1mm` in worktree `/mnt/pikachu/STAR-suite-fix-sampletag-20260905`
(base `release/integrate-v1.8.4` @ `4136243`). One commit landed (`da1a765`, sample-tag fix);
the rest is an uncommitted working set of 10 files (see "What landed"). Final binary
`core/legacy/source/STAR` sha256 prefix `e29d91ed6387049d`. Nothing pushed.

Comparator: Cell Ranger 9.0.1 on the JAX 16-plex set (2,011,130,186 read pairs, 4 declared
samples BC004/BC006/BC007/BC008), `/mnt/pikachu/benchmark_cr9_flex_full/outs/per_sample_outs`.
Concordance scripts: `docs/benchmarks/jax_matrix_20260904/concordance_vs_cr.py` (raw counts,
common cells, genes matched by Ensembl ID) and `paper_protocol_concordance.py` (the manuscript's
log1p / gene-filter protocol). Every arm below is a JAX CBQ no-align run on NVMe with
`/mnt/pikachu/nvme_jax_v184/run_arm.sh ARM [extra flags]`; outputs under `/storage/jax_tagfix/ARM`,
logs `/mnt/pikachu/nvme_jax_v184/arm_*.log`, per-arm record `/mnt/pikachu/nvme_jax_v184/tagfix_results_20260905.txt`.
"hash" = sha256 manifest of every `Gene/filtered` file over all 16 tags (`manifest()` in `run_arm.sh`).

## 1. Result

| build | UMIs ours/CR (BC004/006/007/008) | per-cell Pearson, median | Spearman | Jaccard | wall (CBQ) |
|---|---|---|---|---|---|
| 1.8.4 | 0.9737 / 0.9757 / 0.9764 / 0.9731 | 0.9943 | 0.9914 | 0.9847 | 3:20 (protocol run) |
| **final** | **0.9962 / 0.9971 / 0.9972 / 0.9941** | **0.9964** | **0.9954** | **0.9871** | 3:26 (validation run) |
| cyto 0.4.7 | - | 0.9969 | 0.9953 | 0.9846 | 3:46 (protocol run) |

Paper protocol (log1p, genes >=20 counts in both and detected in >=1% of cells; mean over samples):
per-cell Pearson 0.99168 -> **0.99495** (cyto 0.99522); Spearman 0.99046 -> **0.99453** (cyto 0.99453);
Jaccard 0.9823 -> 0.9833 (cyto 0.9836). Marker genes (WT-Day-7, ours/CR): ZIC2 18,628 -> 24,075 / 24,769.

Gates on the final binary: JAX CBQ `e8cfa27ef18cd573`; BGZF FASTQ input `e8cfa27ef18cd573` (identical);
with the original cache and `--soloProbeMismatch 0 --soloSampleTagMismatch 0 --soloSampleSearchNearby no`
the four real tags are byte-identical to 1.8.4 (`B_fix_off` = `2d2fead4b43f2403`; the sample-less tags
are now empty, so the 16-tag manifest is `23ab4485e2273603`). Small-sample tag test: PASS.

## 2. Where the 2.4% went (decomposition)

Method: per-cell / per-gene / count-bin ratios on the filtered matrices (`residual_variance.py`),
Cell Ranger `molecule_info.h5` joins, a read-level join against CR's BAM (tag classes), FASTQ
scans of 20M pairs classifying every read by Hamming distance to the nearest probe with a tolerant
matcher, a replay of our real H0/H1 cache verdict over those reads, base-quality profiles by
class, and a decode of the cache itself. Files: `/mnt/pikachu/nvme_jax_v184/residual_variance_*`,
`mol_level_fix_WT-Day-8.txt`, `probe_hd_molecules_L001.txt`, `cache_verdict_molecules_L001.txt`,
`hole_genes_scan_L001.txt`, `hd_quality_L001.txt`, `cache_deny_WT-Day-8.txt`,
`residual_variance_summary_20260905.md`.

1. **Sample tag gate (~1.7 points).** We required an exact match to one of the 128 listed 8-mers;
   CR tolerates one substitution. 2.56% of CR's counted molecules had no read we accepted; of those
   nearly all were one substitution from a unique sample.
2. **Terminal-position DENY in the H1 cache (~0.4 points, plus ZIC2 -18%).** 92% of probes had every
   substitution at positions 0/1/48/49 deny-listed (598k of 894k DENY records; interior positions
   ~4%, the genuine paralog conflicts). The synthetic-read validator soft-clips a terminal mismatch
   and cannot resolve the probe, so it returned DENY. Costs every first-cycle sequencing error and
   any molecule with a real variant there: ZIC2 probe f249247 has A->G at position 0 in 78% of its
   reads and that variant was a DENY record.
3. **N in the probe window (~0.05 points).** 1.1% of reads; the encoder rejected them; CR aligns
   through an N.
4. **Reads 2-5 substitutions from a probe (~0.4 points if taken, not taken - see section 4).**
   Fixed-position haplotypes: N4BP2 1e9d19b has no HD0 reads (differs at 30/35/36/44), H2BC15 a
   4-position haplotype in 44% of reads, APOE 5d8ed20 the genotype SNP at 33 plus one error, SOX2
   40285dc a cluster at 42-48. Molecule-inherent, so they remove whole molecules.
5. **Cell sets (Jaccard 0.985 -> 0.987).** CR-only cells (51/57/16/38 per sample, median 1.5-5.6k UMIs)
   are in our raw composite matrix at 0.99x CR's counts and are removed by our EmptyDrops, which tests
   every candidate (`nSimpleCells = 0`, no OrdMag auto-pass) and defines ambient as UMI <= 100 within
   the retain window rather than CR's rank window 45k-90k. Our extra cells (20-104) are near-threshold
   tail rescues CR does not make. Separately, the twelve tags without a sample produced thousands of
   1-3-UMI "cells" (BC001: 2,293 cells holding 2,645 UMIs) through the Simple-ED fallback.
6. **Ruled out:** whitelist (byte-identical), CB correction (raw ratio = filtered ratio), UMI dedup
   (flat across count bins), residual alignment (align arm gives identical concordance; 98.5% of the
   reads it sees come back "too short"), tag gate after the fix (we now reject 4.2% of input for the
   tag, CR 5.0%), H0 deny-listing (all 53,459 probes KEEP at H0).

## 3. What landed (uncommitted unless noted)

| change | files | notes |
|---|---|---|
| Sample tag one-mismatch tolerance (committed `da1a765`, Sol) | `flex/source/SampleDetector.{h,cpp}`, `core/legacy/source/FlexPipeline.cpp`, `Parameters.cpp`, `ParametersSolo.{h,cpp}`, `parametersDefault` | Exact stage first; on a miss, a split-key triple lookup (three CSR piece tables keyed on the two unchanged 3/3/2-base pieces) finds listed sequences within one substitution; two owners -> dropped; a tie at offset 68 is final; nearby offsets +-1/+-2 exact-only. `--soloSampleTagMismatch 0|1`; Flex preset 1 and `soloSampleSearchNearby yes`. Also `call_once` on the static canonical/label tables in `loadWhitelist` and both setters: 32 workers assigned them concurrently, which segfaulted one run in three. |
| Cache generator: terminal-position H1 variants inherit the probe's gene | `core/legacy/source/FlexHashCacheGenerate.cpp` | Positions 0,1,48,49 skip the synthetic validation (which cannot resolve a soft-clipped terminal mismatch); the merge step still denies any variant shared with another gene's probe. DENY records 894,045 -> 201,612. **Requires cache regeneration** (`--runMode hashCacheGenerate --hashCacheTiers H0,H1`, 45 s; new file `/home/lhhung/jax_stage_20260903/ref/h01_cache_edgefix.bin`, the staged `h01_cache.bin` untouched). Old caches keep the artifact silently: a cache-version bump or a load-time warning belongs in the release. |
| Single N in the probe window resolved through the cache | `core/legacy/source/FlexHashScreen.{h,cpp}` (`classifyCbqH0H1Offset0SingleN`, `classifyReadH0H1Offset0SingleN`), `FlexPipeline.cpp` (both FASTQ fused sites and the CBQ site) | The four substitutions at the N are looked up in H0/H1; exactly one gene and no DENY -> Keep; anything else -> Pass (a miss, never a certified negative, so residual alignment may still try the read). CBQ supplies the N position in `n_mask` with the lane zeroed; FASTQ scans the window. |
| `--soloProbeMismatch 0|1` (default 1) | `ParametersSolo.{h,cpp}`, `Parameters.cpp`, `parametersDefault` (+ regenerated `.xxd`) | 0 = exact cache only = 1.8.4 behaviour; 1 = the single-N resolution. Logged at start. |
| Simple-ED floor | `flex/source/libflex/FlexFilter.cpp` | Simple-ED (OrdMag fallback) passers need UMI >= max(retain threshold, `umiMin` = 500); the OrdMag threshold has no floor and lands at 1 UMI (0 on failure) in a tag without a sample. Empties the sample-less tags; those junk cells had been inflating GEM occupancy (see 4.3). |

## 4. What was tried and not landed

### 4.1 Tolerant probe match on cache misses (removed - design decision)
Built in `FlexHashScreen` behind `soloProbeMismatch 2..5`: quarter-piece index (12/13/12/13 bases)
over the 53,459 included probes, candidates verified by full-window Hamming distance, kept when
the best probe is within the cap and no other-gene probe is within one of it; later a per-piece
presence bitset (20 MB) and gate-first ordering that made it cost nothing measurable.

| arm | cap / gate | hash | UMIs ours/CR | Pearson | Jaccard | note |
|---|---|---|---|---|---|---|
| G | 5, none | 609fbfb18e5102a9 | 0.997/0.996/0.990/0.997 | 0.99759 | 0.9761 | occupancy filter removed 41-78 real cells per tag (4.3) |
| H | 5, none, + floor | a8b10331e8ba067f | 1.008/1.008/1.009/1.006 | 0.99759 | 0.9873 | +0.8% over CR, uniform |
| I3 | 3 | ce59f1c9626fc0df | 1.004/1.005/1.005/1.003 | 0.99734 | 0.9874 | |
| I2 | 2 | 28fff45cd728bb72 | 1.001/1.002/1.002/1.000 | 0.99701 | 0.9870 | loses N4BP2/APOE haplotypes |
| J | 5, UMI Q20 gate | 7a0f08a944af6b1e | 1.000/1.002/1.002/0.999 | 0.99702 | 0.9868 | parity; haplotypes kept |
| P | J + prefilter | 7a0f08a944af6b1e | = J | = J | = J | 151.9M reads/10 s vs 150.1M without the matcher |

Why the ungated excess: reads 2-5 substitutions from a probe are low-quality reads - 57-61% carry a
UMI base below Q20 (exact-match reads: 5.7%), only 60-65% have a whitelisted CB (95%) - and a
mis-read UMI mints a molecule instead of joining the real one (2.4 reads per new molecule against the
library's 7.4). The Q20 gate on the UMI fixed that.

Why removed: a runtime Hamming match is not verified against the genome. The H0/H1 cache is
trusted because every entry was aligned once and shown to resolve to its own gene; a tolerant hit
gets no such check, so a read from a paralog or pseudogene transcript a few substitutions from a
probe would be assigned to that probe (H2BC15 at +70% over CR is the likely example). Verifying the
HD2 neighbourhood per probe at generation, or aligning each tolerant hit in align mode, were
discussed; the decision was to leave the step out: the validated cache is the alignment proxy, and
after UMI collapse the verified part of the gain is small (~0.4 points) while the concentrated part
(the haplotypes) is exactly what needs the genome check. The code was removed, not left as an opt-in.

### 4.2 OrdMag auto-pass in the cell filter (optional, not default)
`--soloFlexUseSimpleED yes` (existing switch: forces the Simple-ED/OrdMag bootstrap stage and unions
its above-knee barcodes with the EmptyDrops passers). Arm K on top of J: +6/0/7/25 cells, 37 of 38
being CR cells with median 6-7k UMIs; Jaccard 0.9868 -> 0.9889; counts unchanged; **+55 s wall**
(4:32 vs 3:37) from the bootstrap. Left optional. The real fix is a port of what the non-Flex path
already has (bootstrap OrdMag as primary caller, `CR_COMPAT_PARITY_RESULTS_UCSF_iPSC2_20260221.md`;
guarded rank-window ambient, `HANDOFF_UCSF_EMPTYDROPS_RERUN_20260402.md`) - explicitly left out of
`FlexFilter.cpp` at the time because EmptyDrops and sample resolution are tied per tag in Flex. Still
a todo; needs the cheap OrdMag form (no bootstrap when expected cells are known) before it can be default.

### 4.3 Occupancy filter interaction (root-caused, fixed by the floor)
With the tolerant matcher, the sample-less tags' junk "cells" grew 10-20%; a real BC007 cell whose
barcode appeared as a 1-3-UMI "cell" in 9-15 other tags then exceeded the 99.9% simulated GEM
occupancy and the occupancy guard removed it from every tag (76 of 76 lost WT-Day-8 cells were
such barcodes; median 16k UMIs). The Simple-ED floor removes the junk cells and with it the effect.

### 4.4 Diagnostics that did not run to completion
- Four-tag runs (`--soloFlexAllowedTags` restricted to the declared samples): failed at parsing
  (duplicate flag); superseded by the floor. `run_arm.sh` now takes `ALLOWED=<file>` for this.
- Residual-alignment seeding diagnostic (why 98.5% of misses sent to the aligner return "too short"
  although the Flex preset zeroes the per-length thresholds and sets `outFilterMatchNmin 25`): a
  2M-read align-mode run with `STAR_TRIM_DEBUG_N` and denser seed starts was started and stopped
  on request - the decision is not to chase the alignment stage.

## 5. Arms and hashes (all CBQ no-align unless noted; validation runs, not timing protocol)

| arm | binary | cache | hash | wall | KEEP % | DENY % |
|---|---|---|---|---|---|---|
| B_fix_off / R / R2 (1.8.4 behaviour) | various | old | 2d2fead4b43f2403 (16 tags) / 23ab4485e2273603 (sample-less tags empty) | 3:20-3:22 | 90.84 | 0.94 |
| A_fix_on (tag fix, table impl.) / Sol's triple lookup | 1632b85d / 102cd3f3 | old | 23cd7d3c0cc9a00f | 3:25 | 94.39 | 0.94 |
| C_fix_on_align (align mode) | 102cd3f3 | old | 8f95174ff0db3581 | 13:56 | | |
| D_fix_on_bgzf (BGZF input) | 102cd3f3 | old | 23cd7d3c0cc9a00f | 3:59 | | |
| E_edgefix | 83da82fc | edgefix | a5d46028ccb76c67 | 3:31 | 95.21 | 0.12 |
| F_edgefix_N | ba9553b4 | edgefix | dc26a94366ba1092 | 3:28 | 95.22 | 0.12 |
| **N_final** (shipped) | **e29d91ed** | edgefix | **e8cfa27ef18cd573** | **3:26** | 95.22 | 0.12 |
| **M2_final_bgzf** (BGZF input) | e29d91ed | edgefix | **e8cfa27ef18cd573** | 3:46 | | |
| G/H/I/J/K/L/P (tolerant variants, removed) | see 4.1 | edgefix | see 4.1 | 3:27-4:32 | 95.27-95.32 | 0.11 |

Per-arm concordance tables: `/mnt/pikachu/nvme_jax_v184/concordance_*.txt`, `paperproto_*.txt`.

## 6. Open items
1. Commit the working set on `fix/sample-tag-1mm` (single commit or split: generator / N / filter);
   no attribution trailers. Version bump and release note: the Flex preset now changes default output
   (tag tolerance, N, floor), and the cache must be regenerated.
2. Cache-format version bump or load-time warning for caches generated before the terminal-position fix.
3. Re-time the CBQ/BGZF rows (no-align and align) with the quiet-box protocol on this binary:
   `run_jax_nvme_v184.sh` must first lose its hard-coded `--soloSampleSearchNearby no`.
4. Paper accuracy numbers and Supplementary Note 1 from `N_final` (both concordance scripts).
5. OrdMag-primary + rank-window ambient port for the Flex filter (4.2).
6. Known remaining differences: interior paralog-conflict DENYs (EIF2A -30%, NNAT, H1-10, HNRNPA3,
   PA2G4 - the deliberate genome stringency); tail cell calling; H2BC15 (CR under-counts a 4-mismatch
   haplotype; cause not established without CR's BAM).

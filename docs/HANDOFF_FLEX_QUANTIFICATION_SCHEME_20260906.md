# Handoff: how STAR-Flex is quantified against Cell Ranger (the scheme used 2026-09-04..06)

Use this scheme unchanged so numbers stay comparable across sessions and agents. Scripts live in
`docs/benchmarks/jax_matrix_20260904/` (concordance) and `docs/benchmarks/jax_matrix_20260904/analysis/`
(decomposition, harnesses). Results of the sessions that used it: `HANDOFF_FLEX_CONCORDANCE_FIXES_20260906.md`.

## 0. Rules that apply to every number
- **Two kinds of runs.** *Validation runs* (correctness, hashes, concordance) may run whenever the box is
  free of other benchmarks; their wall times are indicative only (±5 s) and are never quoted as timings.
  *Timing runs* follow the written protocol (`run_jax_nvme_v184.sh`): final binary pinned by sha, staged
  NVMe inputs, page cache dropped per arm, quiet box, nothing else running. Do not mix the two.
- **Inputs on the SSD, including `/tmp` and outputs.** `/mnt/pikachu` is rotational RAID and is never
  read or written during a run.
- **Identify every run by binary sha and cache file.** `sha256sum STAR | cut -c1-16`; the cache path is
  part of the configuration (regenerated caches are new files, the staged one is never overwritten).
- **Output identity is a manifest hash**, not a file hash: sha256 over `sha256sum` of every
  `per_sample/*/Gene/filtered/*` file except `gdna_metrics.json`, sorted (`manifest()` in `run_arm.sh`).
  It covers all 16 tags, so a change confined to sample-less tags changes the hash; when that happens,
  compare the four real tags per tag (`cat matrix.mtx barcodes.tsv features.tsv | sha256sum` per tag).
- **Genes are matched by Ensembl ID** (column 0 of `features.tsv`). Cell Ranger names all but 25 genes by
  symbol in column 1; matching on symbols intersects 25 genes and silently produces meaningless
  correlations (this bug invalidated an earlier report).
- **A cell is CB16 + tag** (`cell_key` -> `CB16|BC0xx`), never CB16 alone - two tags in one GEM are two cells. Formats: CR per-sample barcodes are
  `CB16 + probe-barcode8 + "-1"`; our filtered barcodes are 16-mers; our `Solo.out/Gene/raw` composite
  barcodes are `CB16 + tag8` (the tag's listed sequence, e.g. BC007 = AAGTAGAG).
- **Report medians over the four samples and the per-sample rows**; never a single sample as "the" number.

## 1. Concordance vs Cell Ranger (what goes in the paper)
Comparator: `/mnt/pikachu/benchmark_cr9_flex_full/outs/per_sample_outs/<sample>/count/sample_filtered_feature_bc_matrix`
(CR 9.0.1, JAX defaults, 4 declared samples). Sample <-> tag map: `DEFAULT_GROUPS` in `concordance_vs_cr.py`
(WT-Day-7=BC004, PAX6-PTC-D9-Day7=BC006, WT-Day-8=BC007, PAX6-PTC-D9-Day8=BC008).

1. **Raw-count protocol** — `python3 concordance_vs_cr.py <per_sample dir> <label>` (or `--input-kind cyto`
   for cyto's `counts/` output). Per sample: cells called by each tool, **Jaccard** of the cell sets,
   and over *common cells only* with Ensembl-matched genes: per-cell **Pearson** (median, 1st percentile,
   worst cell), per-cell **Spearman** (median), **gene-total Pearson**. Cells called by one tool only are
   excluded from the correlations; the Jaccard column is where that disagreement shows. The `r>=1-1e6`
   column is a legacy of the totals-correlation era and is always 0 here.
2. **Paper protocol** — `python3 paper_protocol_concordance.py <per_sample dir> <label> [--input-kind cyto]`:
   log1p counts, genes with >=20 total counts in *both* tools and detected in >=1% of common cells in both,
   per-cell Pearson/Spearman mean and median, gene-level Pearson on log1p totals. This is the manuscript's
   stated Methods protocol; quote it alongside the raw-count protocol, never instead of it.
3. **UMI totals** — `UMIs ours/CR` = sum of our filtered matrix / sum of CR's filtered matrix per sample
   (each tool's own called cells; `run_arm.sh` prints it). Also computed over *common cells only* by
   `residual_variance.py` ("UMIs in common cells: STAR/CR"). Say which one you mean.
4. **Totals vs profiles.** Correlations of per-cell *totals* (0.9999x) are not per-cell *profile*
   correlations (0.99x); the paper's older 0.99998/0.99993 are totals correlations. Do not mix them in
   one table.

## 2. Validation harness (`analysis/run_arm.sh`)
`./run_arm.sh ARM [extra STAR flags]` runs JAX CBQ no-align (32 threads, all-NVMe) with the fix
worktree's binary, then prints: exit, wall, manifest hash, `Hash screen: KEEP/DENY/unmatched sample tag %`,
the raw-protocol concordance median row, and per sample `UMIs ours/CR` plus a few marker genes
(this/E/CR). `CACHE=<file>` and `ALLOWED=<tag file>` override the cache and the allowed-tag list.
Sentinel `$T/arm_ARM.done`; run it detached (`setsid nohup ... & disown`) and wait on the sentinel.
Never overlap two arms. For a BGZF-input identity check use the `common` array with
`--readFilesBgzfMode range --bgzfReaderThreads 32 --bgzfCrcCheck 1 --readFilesIn R2s R1s`
(see `run_edgefix_test.sh` / the M arms) - its hash must equal the CBQ arm's.

Regression gates every change must pass:
- **1.8.4 identity**: old cache + every new behaviour switched off (`--soloProbeMismatch 0
  --soloSampleTagMismatch 0 --soloSampleSearchNearby no`) -> the four real tags byte-identical to
  `B_fix_off` (16-tag manifest `2d2fead4b43f2403`, or per-tag hashes when the sample-less tags differ).
- **FASTQ/CBQ identity**: BGZF arm hash == CBQ arm hash for the same binary and flags.
- **Small tag test**: `analysis/tagtest/run_tagtest.sh` (10,000 real read pairs rewritten into five tag
  classes; expected KEEP/UNMATCH per class with the fix on and off; PASS/FAIL).

## 3. Decomposing a residual difference (`analysis/*.py`)
Run in this order; each answers one question. All read the filtered/raw MEX directly (no scanpy).

| script | question | key outputs |
|---|---|---|
| `residual_variance.py <per_sample> [molecule sample]` | Is the difference per-cell, per-gene, per-molecule? | per-cell total ratio percentiles (uniform vs CB-driven); per-gene ratio distribution and largest deficits/excesses with symbols; **count-bin table** ((cell,gene) pairs binned by CR count: a ratio that falls with count = dedup over-merging, flat = per-read loss); presence/absence pairs; cells only in one tool with their totals in the other tool's *raw* matrix (threshold effect vs real loss); optional molecule_info join |
| `mol_level.py <per_sample> <sample> <tag>` | Are the missing molecules singletons (read-level loss) or multi-read (molecule-inherent)? Which probes carry them? | read-depth profile of molecules in equal/deficit/excess pairs; per-probe share in deficit pairs; hole genes' probes |
| `probe_hd_molecules.py [N pairs]` | For real reads: Hamming distance to the nearest probe (tolerant matcher: exact 25-mer halves, then 1-mismatch halves) aggregated per molecule (CB,UMI,tag) | per-read class %, molecules with no read within HD1, GC dependence, genes with most lost molecules |
| `hole_genes_scan.py` | For named genes: are mismatches at fixed positions (variant/paralog/probe design) or random (sequencing error)? | HD histogram and top mismatch positions per probe |
| `cache_verdict_molecules.py [N pairs]` | What would *our* cache actually say for each read? | replay of H0/H1 verdicts (KEEP / DENY / miss-but-HD2-5 / N / junk) per read and per molecule |
| `hd_quality.py` | Are reads in a class trustworthy? | mean base quality of probe/UMI/CB, fraction with any UMI or CB base < Q20, CB-in-whitelist rate, per HD class |
| `cache_deny.py` | Which probes/variants does the cache deny, and how many CR molecules sit on them? | cache record counts by class, per-probe H0 status, DENY-by-position (edge vs interior), CR molecules on denied probes |

Cache decoding, for any of these: header `<magic 8s, version u2, k u2, recsize u4, n u8>` (24 bytes),
records `<lo u8, hi u8, gene u4, cls u1, neg u1, res u2>` (24 bytes); gene = `gene & 0x7FFF`; class 0 = H0,
1 = H1, 2 = DENY (`neg != 0`). Encoding of a 50-mer: base 0 in the most significant position,
`int(seq.translate(ACGT->0123), 4)`, hi = bits 64..99, lo = bits 0..63. Probe sequences:
`star_index/flex_probe_artifacts/filtered_probe_set.csv` (included=TRUE rows), read orientation.

Read layout (R2): probe window bases 0-49, constant region 50-67, sample tag 68-75. R1: CB 0-15, UMI 16-27.

Cell Ranger side: `sample_molecule_info.h5` (`barcode_idx`, `feature_idx`, `probe_idx`, `count` = reads per
molecule, `umi`; `probes/probe_id` = `ENSG|SYMBOL|hash`), `sample_raw_feature_bc_matrix` for
"is this barcode present at all", `metrics_summary.csv` for read-assignment fractions. There is no CR BAM
for the JAX run on this box (the earlier read-level tag join used an 8M-record BAM that is no longer
present); read-level questions are answered from the FASTQ with the scripts above.

Our side beyond the filtered matrix: `Solo.out/Gene/raw` (composite CB16+tag8 barcodes, all tags) for
"did we accumulate this barcode at all" (column sums via `awk` on `matrix.mtx`), `per_sample/flexfilter_summary.tsv`
(Expected, Retain, Simple_ED, Tail_Tested, ED_Pass, Occ_Rem, Final, Total_UMIs per tag), `Log.final.out`
(`FLEX HASH SCREEN` block: evaluated, KEEP, KEEP no barcode, DENY, unmatched sample tag, PASS),
`Log.out` (`PIPELINE_STATS t=..s ... delta10s=` lines; `Flex sample tag detection:` and
`H0/H1 hash screen:` lines at start), `Solo.out/Barcodes.stats` (CB match classes).

## 4. Throughput and memory
- Streaming throughput = mean of `delta10s` over the `PIPELINE_STATS` lines with the first and last dropped
  (reads per 10 s at steady state). Use it to attribute a cost to a code change; wall on a validation run
  cannot (tail and noise are ±5 s). A change is "free" when this number is within noise of the baseline.
- Wall, peak RSS and CPU% come from `/usr/bin/time -v` (`Elapsed`, `Maximum resident set size`, `Percent of CPU`).
- Timing rows for the paper come only from the protocol script, never from `run_arm.sh`.

## 5. Interpreting the pieces
- A **uniform** deficit (per-cell IQR narrow, per-gene median ~= total ratio, count-bin flat) is a
  read-acceptance difference upstream of counting; check the tag gate, the hash tiers, N handling.
- A deficit **growing with count** is dedup; a **wide per-cell** spread is barcode correction/assignment.
- Deficit molecules that are **multi-read** are molecule-inherent (variants, probe design, DENY records);
  singleton-only losses are sequencing noise.
- **Gene holes** are diagnosed by `hole_genes_scan.py`: fixed mismatch positions = variant/paralog/design;
  then `cache_deny.py`-style lookups of the observed variant tell whether the cache denies it.
- **Cells only in CR** must be looked up in our raw composite matrix before any conclusion: present at
  ~CR's count = our filter removed them; absent = never accumulated (whitelist/tag/CB path).
- **Reads accepted by a new rule** must be checked with `hd_quality.py` before trusting the molecules
  they create: low-quality reads mint spurious molecules from mis-read UMIs and read as a uniform *excess*.
- **Sample-less tags** (tags in the whitelist with no sample) can carry junk cells that distort GEM
  occupancy; check `flexfilter_summary.tsv` for every tag, not only the compared ones.

## Correction (2026-09-06): grouped-sample identity bug in the comparator scripts
`concordance_vs_cr.py` / `paper_protocol_concordance.py` (and the native port that reproduces them)
truncated every barcode to CB16 and then summed rows with equal barcodes, so in a sample that spans
more than one tag (the 320K two-tag samples) two independent cells sharing a GEM were merged into one.
Single-tag JAX validation was unaffected (every number in this document stands). **All grouped-sample
reports produced with the Python scripts or the native comparator before this fix are invalid for
publication and must be recomputed.** Fixed in the same commit as this note: a cell is keyed by
CB16|tag (`cell_key`); STAR/cyto per-tag inputs take the tag from the directory/file name, Cell
Ranger's CB16+tag8 barcodes are mapped back to the tag id through the probe-barcode whitelist
(`--tag-map`). The native comparator must adopt the same key before it is used on grouped samples.

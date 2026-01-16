# SLAM Method Technical Note: Metrics and Benchmarks

This note captures the key SLAM metrics used for interpretation and testing, and
the published inter-method correlation ranges that provide a sanity check for
cross-tool comparisons. It is intended to inform STAR-SLAM implementation and
parity testing against GRAND-SLAM.

## Core GRAND-SLAM outputs (GEDI wiki)

Source: https://github.com/erhard-lab/gedi/wiki/GRAND-SLAM

- Per-gene read count (fractional for multimappers).
- New-to-total RNA ratio (NTR): reported as MAP and mean, with 0.05/0.95
  quantiles and Beta posterior parameters (alpha, beta).
- Mismatch summaries for exonic/intronic and sense/antisense, single vs double,
  and position-specific profiles to diagnose trimming artifacts.
- Binomial mismatch rates (p_c) for new vs old RNA used by the EM model.

Preprocessing guidance emphasizes:
- Mapped reads must include `MD` (and preferably `NH`) tags.
- Avoid soft-clipping; use end-to-end mapping and trim 5' and 3' to remove
  read-end artifacts detected via mismatch-position plots.
- If reads are already fixed-length and adapter-free (e.g., SE50), disable
  adapter clipping (`--clip3pAdapterSeq` / `--clip3pAdapterMMp`) and keep
  `--alignEndsType EndToEnd`.

## Published inter-method comparison benchmarks

The GRAND-SLAM paper itself does not publish explicit "acceptable" correlation
thresholds. A relevant external benchmark is the comparison study:

- Boileau et al. 2021, "A comparison of metabolic labeling and statistical
  methods to infer genome-wide dynamics of RNA turnover." (Briefings in
  Bioinformatics)
- Analysis scripts and tables are in:
  https://github.com/dieterich-lab/ComparisonOfMetabolicLabeling

This study compares decay-rate estimates across protocols and methods (pulseR
and GRAND-SLAM). Pearson correlations are used, along with residual standard
error and RMSD, and confidence-interval overlap.

Using the published `paper/tables/tbl.xlsx` from that repo, Pearson r between
pulseR and GRAND-SLAM decay-rate estimates across 11,603 genes is:

- Timepoints 0-1-2h:
  - SLAM: r = 0.7722
  - TLS:  r = 0.8076
  - TUC:  r = 0.8360
- Timepoints 0-2-4h:
  - SLAM: r = 0.8629
  - TLS:  r = 0.8818
  - TUC:  r = 0.8776
- Timepoints 0-4-8h:
  - SLAM: r = 0.7810
  - TLS:  r = 0.8138
  - TUC:  r = 0.8104
- Timepoints 0-1-2-4-8h:
  - SLAM: r = 0.8856
  - TLS:  r = 0.8972
  - TUC:  r = 0.8938

Note: these correlations are for decay-rate estimates (not NTR). They are
useful as a sanity range for cross-method agreement but are not a direct
acceptance threshold for STAR-SLAM vs GRAND-SLAM parity.

## Controls and "gold set" references for SLAM-seq

There is **no universal public gold-standard SLAM-seq spike-in set**. The
original SLAM-seq paper validated conversion detection using controlled inputs
and reported public datasets, but not a standardized benchmark panel.

Documented controls from the literature:
- **Synthetic RNA oligos with defined s4U** (conversion ground truth):
  - Herzog et al. 2017 (Nat Methods, SLAM-seq) used let-7 RNA oligos, including
    a single-s4U position (`let-7-s4Up9`) to validate conversion behavior.
  - Sequences are in the paper’s supplementary tables (not a packaged dataset).
- **Spike-in RNA for small RNA libraries**:
  - The same paper added *Arabidopsis thaliana* total RNA as a spike-in before
    size selection for small RNA libraries (library-normalization control).
  - This is **not** a SLAM conversion gold set.
- **Public SLAM-seq datasets**:
  - Sequencing data from the SLAM-seq paper is available in GEO:
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99978
  - The paper notes that UTR/counting-window BED files are included in the GEO
    datasets, which are relevant for reproducing their quantification choices.
- **Protocol-level QC guidance**:
  - The SLAM-seq protocol paper (Methods Mol Biol 2020, PMCID: PMC7611866)
    mentions low-depth spike-in approaches for QC but does not publish a
    standard spike-in dataset.

Practical implication: if we need a true gold set, we should generate one
internally (e.g., defined mixtures of labeled/unlabeled IVT RNA or synthetic
oligos) and treat public SLAM-seq datasets as **reference benchmarks**, not
ground truth.

## Implications for STAR-SLAM testing

- For direct parity, compare NTR and conversion-related metrics on a common
  gene set with adequate coverage (e.g., nT >= 20-50).
- Expect high correlation for well-covered genes when preprocessing is matched
  (annotation, SNP mask, trimming, mapping mode).
- Systematic deviations are often concentrated in mitochondrial and ribosomal
  genes or low-complexity regions; document these explicitly.
- If the library is stranded, set `--slamStrandness Sense` (or Antisense) in
  STAR-SLAM to mirror GEDI's strand-aware model; default `Unspecific` keeps all
  reads regardless of strand.

## Debug instrumentation (STAR-SLAM)

To pinpoint divergences, use targeted debug outputs:

- `--slamDebugGeneList` and/or `--slamDebugReadList` with `--slamDebugOutPrefix`
  to emit `<prefix>.gene.tsv` and `<prefix>.reads.tsv`.
- The gene file includes per-gene weighted counts plus strand/intronic splits
  and drop counts (SNP mask, strandness).
- The read file captures per-read assignments, drop reasons, and conversion
  positions for the selected genes/reads.

## Fixture parity benchmarks (STAR-SLAM vs GRAND-SLAM)

These values come from the 100k-read human fixture and the GRAND-SLAM reference
`fixture_ref_human.tsv.gz`, using `tests/slam/compare_fixture.py` with
read-count thresholds (>=20, >=50, >=100). They provide a concrete baseline for
expected parity in this repository.

### SNP mask from BED (slamSnpBed)

This is the most direct parity mode because the reference fixture was generated
using the same SNP mask. Use a BED only when a sample-specific VCF/BED is known;
otherwise prefer internal SNP detection for production data.

| Threshold | Filter | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |
|-----------|--------|---------|-------------|--------------|--------------|---------------|
| >=20      | readcount | 384 | 0.998865 | 0.981026 | 0.999865 | 0.999489 |
| >=50      | readcount | 76  | 0.996069 | 0.986681 | 0.998412 | 0.996712 |
| >=100     | readcount | 23  | 0.994379 | 0.994066 | 0.996847 | 0.994071 |

### Internal SNP detection (slamSnpDetect 1)

This mode is expected to diverge slightly from the BED-based reference (different
SNP calling), but still provides strong agreement. This is the recommended
default when no external SNP set is available.

#### Auto-estimated SNP mismatch threshold

When `--slamSnpDetectFrac` is not set or `<= 0` (default), STAR-SLAM automatically
estimates the mismatch fraction threshold using a knee/elbow detection algorithm:

1. **Build histogram**: For each genomic position with coverage >= 10, compute
   `f = mismatches / coverage` and bin into 100 bins across [0, 1].

2. **Cumulative distribution**: Compute `N(>= f)` from high to low bins.

3. **Knee detection (Kneedle-style)**:
   - Apply `log1p(N)` transformation for stability.
   - Normalize x (bin centers) and y (log counts) to [0, 1].
   - Find the bin maximizing distance from diagonal `y = 1 - x`.

4. **Guardrails**:
   - Minimum 1,000 eligible sites required, otherwise fallback.
   - Knee strength must exceed epsilon (0.02), otherwise fallback.
   - Result clamped to [0.10, 0.60] to avoid pathological thresholds.

5. **Fallback**: If auto-estimation fails, uses 0.22 (GEDI-parity value).

For reproducibility, set an explicit threshold: `--slamSnpDetectFrac 0.22`.
The log output reports the mode used (`auto`, `auto_fallback`, or `explicit`)
along with the estimated value and knee bin when auto-estimation succeeds.

| Threshold | Filter | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |
|-----------|--------|---------|-------------|--------------|--------------|---------------|
| >=20      | readcount | 384 | 0.989764 | 0.967831 | 0.991470 | 0.990457 |
| >=50      | readcount | 76  | 0.996314 | 0.980295 | 0.998116 | 0.994169 |
| >=100     | readcount | 23  | 0.993684 | 0.991098 | 0.996284 | 0.988142 |

Note: k/nT correlations consistently reach the ~0.999 level under the BED mask,
which reflects very tight agreement in the raw conversion signal. NTR remains
the primary biological endpoint, with correlation values interpreted in the
context of coverage thresholds and SNP handling mode.

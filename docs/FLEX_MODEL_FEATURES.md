# Flex modeling features and filtered exports

The tag-aware caller uses the complete probe count matrix for OrdMag,
candidate selection, ambient estimation, and EmptyDrops. Deprecated features
remain separate `DEPRECATED_` IDs. The per-sample filtered MEX contains only
the included active genes supplied through `--soloFlexFilteredGeneList`.

For the count-only H1X2 path, extend an existing complete probe cache using
the same public probe panel and the exact gene-ID list used to build that
cache:

```bash
python3 flex/scripts/extend_model_probe_cache.py \
  --cache /path/to/active_h01x2_cache.bin \
  --gene-ids /path/to/active_gene_ids.txt \
  --probe-csv /path/to/probe_set.csv \
  --out-dir /path/to/new_model_cache
```

The output directory must be fresh. The builder requires Python and NumPy.
It preserves existing gene indices, appends missing active and deprecated
features, and adds their exact H0 parents and H1X2 variants. Multiple-parent
global variants become DENY entries, including collisions with existing
variants. A unique new exact parent takes priority over an old non-exact
cache verdict. Conflicting exact gene/region assignments fail explicitly.
The manifest records source hashes, feature counts, and cache changes.
No genome index, alignment, Cell Ranger code, or external count matrix is
used to generate these probe records.

Use these resources with the existing Flex count-only command:

```text
--soloHashScreenFile /path/to/new_model_cache/model_h01x2_cache.bin
--soloProbeList /path/to/new_model_cache/model_gene_ids.txt
--soloRemoveDeprecated No
--soloFlexFilteredGeneList /path/to/new_model_cache/included_gene_ids.txt
--soloRunFlexFilter yes
--soloFlexCellCaller tag-aware
--flexNoAlign 1
--outSAMtype None
--soloFlexDecisionSidecar -
```

Keep the existing sample whitelist, barcode whitelist, tag offset, UMI, and
input-reader settings. The explicit `--soloRemoveDeprecated No` is needed
when using presets that otherwise remove deprecated IDs. The gene-ID list
and cache must always be used together; sorting or replacing only one can
change the meaning of numeric cache gene indices. The raw MEX retains all
modeling features; the export allowlist only affects the per-sample filtered
MEX.

The tag-aware defaults are a fixed ambient rank window of 45,000–90,000 per
tag, 100,000 EmptyDrops simulations, BH FDR 0.01, and the existing configured
500-UMI primary/candidate floor. `--soloFlexEdNiters` can explicitly override
simulation depth; zero chooses the caller default. The legacy caller keeps
its 10,000-simulation default. The standalone `scrna_simpleed --flex-tag-aware`
uses the same fixed ambient and simulation defaults; an explicit
`--ambient-umi-target-per-tag` remains available for diagnostics.

After calling all sample groups, occupancy counts distinct observed TAG8
calls per CB16 GEM, fits the positive observations using a zero-truncated
Poisson, and removes GEMs above the ordinary Poisson 0.999 quantile. Empty
tags contribute no calls. Full CB16+TAG8 identities are preserved when
exporting and comparing cells.

The shared-library floor test, grouped caller fixture, observed-occupancy
test, and `tests/emptydrops/test_model_probe_cache.py` cover the relevant
boundaries. Full CBQ/BGZF read validation and concordance artifacts for the
2026-09-09 integration are in
`/mnt/pikachu/star_suite_paper/analysis/full320k_star_deprecated_20260909/`;
consult its completion records and
[full benchmark report](HANDOFF_FULL320K_STAR_DEPRECATED_20260909.md) for
measured results. The complete STAR-read implementation improves all eight
samples, with pooled cell Jaccard 0.954180 → 0.973961 and exact CBQ/BGZF raw
count parity. The tested executable is installed at
`core/legacy/source/STAR`; the prior executable is preserved in the artifact
directory.

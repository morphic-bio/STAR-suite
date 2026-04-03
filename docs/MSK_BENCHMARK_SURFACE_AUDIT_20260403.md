# MSK Benchmark Surface Audit (2026-04-03)

## Purpose

Separate the MSK historical raw-matrix EmptyDrops isolation result from the
actual guarded end-to-end benchmark surface so future reviews do not treat the
`0.997556` isolation Jaccard as the expected production target.

## Surfaces

| Surface | Artifact | Cells | Barcode Jaccard vs CR9 | Gene Pearson vs CR9 | Notes |
| --- | --- | ---: | ---: | ---: | --- |
| Historical raw GEX-only ambient isolation | `docs/MSK_AMBIENT_RESCUE_ISOLATION_20260402.md` | `32,303` | `0.997556` | n/a | Same historical raw MEX, same `libscrna` MC engine, ambient-only swap to legacy `45k-90k` rank window |
| Modern GEX-only guarded bridge rerun | `/tmp/msk_gex_multimap_unique_guarded_20260403/` | `33,092` | `0.974260` | `0.994554` | Same bridge / `Unique` surface as the README scRNA benchmark, but with guarded ambient fallback |
| Modern full perturb guarded rerun | `/storage/MSK-perturb-comparison/paper_bench_emptydrops_guarded_20260403_2/` | `33,095` | `0.974112` | `0.994554` | README-ready 3-library benchmark surface; `27.1` min wall, CRISPR exact-match `97.68%` |

## Key Point

The `0.997556` result is valid, but it answers a narrower question:

- "If the historical raw GEX matrix is held fixed and only the ambient pool is
  corrected, how much Cell Ranger 9 overlap can be recovered?"

It does **not** answer:

- "What should the modern end-to-end MSK benchmark achieve after the ambient
  fix?"

That second question is answered by the guarded modern reruns above, which land
near `0.974`, not `0.997`.

## Knob Audit On The Modern GEX-only Surface

These probes all used the same guarded EmptyDrops policy on the same modern MSK
GEX-only raw-count surface.

Baseline artifact:

- `/tmp/msk_gex_multimap_unique_conventional_guarded_20260403/`

| Surface | Cells | Jaccard vs CR9 | Gene Pearson vs CR9 | Interpretation |
| --- | ---: | ---: | ---: | --- |
| `Unique`, `clip3pPolyG yes`, score `0` | `33,123` | `0.973646` | `0.994539` | Baseline modern GEX-only guarded surface |
| `EM`, `clip3pPolyG yes`, score `0` | `33,077` | `0.974940` | `0.985952` | Slight cell Jaccard gain, clear gene-parity loss |
| `Unique`, `clip3pPolyG yes`, score `30` | `33,117` | `0.973763` | `0.994538` | Effectively negligible |
| `Unique`, `clip3pPolyG no`, score `0` | `33,115` | `0.973941` | `0.965570` | Cell effect is tiny; gene Pearson collapses |

Additional overlap results:

- `score=30` vs `score=0`: Jaccard `0.999034` (`13` score30-only, `19`
  score0-only barcodes)
- `clip3pPolyG no` vs `yes`: Jaccard `0.996022` (`62` no-polyG-only, `70`
  polyG-only barcodes)
- bridge vs conventional `Unique` guarded runs:
  `/tmp/msk_gex_multimap_unique_guarded_20260403/` vs
  `/tmp/msk_gex_multimap_unique_conventional_guarded_20260403/`
  Jaccard `0.997255`

## Conclusions

- The ambient fallback bug was real and explains the large jump from the old
  modern MSK surface (`~0.942`) to the guarded reruns (`~0.974`).
- The historical raw-matrix isolation result (`0.997556`) should remain
  documented, but it is **not** the right benchmark target for the modern
  production surface.
- No single remaining upstream flag recovers `.99+` Jaccard on the modern MSK
  surface.
- `--soloMultiMappers EM` improves barcode overlap only slightly and costs real
  gene-level parity.
- `--outFilterScoreMin 30` is negligible on MSK GEX.
- `--clip3pPolyG` matters strongly for gene Pearson, but not for barcode
  Jaccard.

The remaining gap appears to be a distributed raw-matrix surface difference
between the historical CellGENI-style configuration and the modern CR-compat /
`Unique` counting surface, not one hidden leftover switch.

# Split-pool barcode lists, v3 chemistry (short linkers)

Lists for `--soloSplitseqWhitelist` with the barcode-read layout
`--soloSplitseqCBposition 50_8 30_8 10_8 --soloSplitseqUMIposition 0_10`
(see `docs/SPLITSEQ_RNA.md`).

## Source

The kallisto-bustools documentation repository of the Pachter lab,
<https://github.com/pachterlab/kb_docs>, commit
`a4a4d985a4491264e02aa29807dff06f31e09fec`, licence BSD-2-Clause, fetched
2026-11-16. The two files are copied unchanged into `source/`, with the
repository's licence text as `source/LICENSE.kb_docs`.

| File in kb_docs | Git blob | Copy here |
| --- | --- | --- |
| `docs/barcodes/splitseqv3_barcodes.txt` | `99dd06e6cdb838257bc8a52a0f72cc925d39698a` | `source/splitseqv3_barcodes.txt` |
| `docs/barcodes/splitseqv3_replace.txt` | `9b65031521d1596989459ec0283bf3635a4f27b5` | `source/splitseqv3_replace.txt` |
| `LICENSE` | `598e1ac67da50e5528780302b27b3d2f65e889f8` | `source/LICENSE.kb_docs` |

The on-list has three columns in read order (round 3, round 2, round 1: 96,
96 and 208 sequences). The replacement list has 106 lines "SOURCE *TARGET";
the kallisto-bustools documentation describes it as converting the
random-hexamer barcodes into their oligo-dT counterparts.

## Generated files

```bash
python3 scripts/make_splitseq_lists.py --onlist source/splitseqv3_barcodes.txt \
  --replace source/splitseqv3_replace.txt --wells 96 --out-dir 96_wells
python3 scripts/make_splitseq_lists.py --onlist source/splitseqv3_barcodes.txt \
  --replace source/splitseqv3_replace.txt --wells 48 --out-dir 48_wells
```

- `96_wells/`: round 1 has 202 lines (96 oligo-dT, 96 random-hexamer and 10
  further entries that the replacement list maps to a well). The six round-1
  entries of the on-list that have no replacement line are left out.
- `48_wells/`: the round-1 entries of wells 1 to 48 (rows A to D), 99 lines.
  Use it for a library that used the first 48 round-1 wells; a list with wells
  that were not used lets a sequencing error be corrected into an unused well.
- `round1.tsv` has two columns (observed sequence, sequence it is counted as);
  `round2.tsv` and `round3.tsv` have one. `crosswalk.tsv` gives primer type,
  well and well index of every sequence.

## What was checked, and on what

One public library: GEO GSE352893, run SRR41469160 (97,458,666 read pairs;
the record names the Evercode WT v3 kit; 48 round-1 wells).

- Read off the barcode reads with no list (reads with both linkers exact, the
  most frequent word at each position): 96 words per round. Rounds 2 and 3
  equal the on-list columns; the 96 round-1 words are the oligo-dT and
  random-hexamer entries of wells 1 to 48.
- Which round-1 words share cells, from the first 30,000,000 reads alone:
  counting round 2 and round 3 pairs with at least 20 reads, or at least 50,
  the best partner of all 96 words is the other entry of the same well in the
  replacement list. The result depends on that limit: with at least 5 reads
  it holds for 85 of the 96.
- Well order: cell labels "%02d_%02d_%02d" built from `crosswalk.tsv` well
  indices match the cell names of the deposited matrix
  (`GSE352893_iPSC_K562_mixed_integrated.h5mu`): all 16,818 deposited cells
  are among the barcodes counted, and 16,419 of the 16,818 barcodes with the
  highest counts are deposited cells (about 620 when the well indices are
  shuffled). So the line order of the two files is the plate order, row-major.
- The same 208 round-1 sequences, with the same wells and primer types, are in
  the barcode table `bc_data_n208_R1_v3_3.csv` of
  <https://github.com/mortazavilab/parse_pipeline> (commit
  `b156bbd6`; no licence file; used for comparison only, not copied).

Not checked: wells 49 to 96 of round 1 on data (the `96_wells` list rests on
the source files for those wells); any other library; any other kit version.

## SHA256

| File | SHA256 |
| --- | --- |
| `source/splitseqv3_barcodes.txt` | `a8452cc302ef720c7c1742b60f7aee6857571761aedd50964f401398c74f89b7` |
| `source/splitseqv3_replace.txt` | `303a9170b82626146db1ea25856b7c28f8d8690e032577d24e175dd327481321` |
| `source/LICENSE.kb_docs` | `61a194af4bc41f610c251ca1a0e4f80d7ed70361158055126a817bed566f1893` |
| `96_wells/round1.tsv` | `0626a7c350e10e20e0ed35aea6f7fb39d9b47c3bac0ff4b2bc871260ad6c19d9` |
| `96_wells/round2.tsv` | `346ecc6bf50e72c91e0ec115ca22918abbfa8d3dae7be0dab9e602a1564de199` |
| `96_wells/round3.tsv` | `35e18a9c0b63befe55320c25384232cecabbed7211a56454f6d77eef5e192b34` |
| `96_wells/crosswalk.tsv` | `3d8df4e3e6b303d054df9dde4284c3d49b37c0876d09ebb689479aa03ee4b528` |
| `48_wells/round1.tsv` | `761c5071c0f1651966fd7350106e255b141b9edbadb080c74f7d48bc5030f4c8` |
| `48_wells/round2.tsv` | `346ecc6bf50e72c91e0ec115ca22918abbfa8d3dae7be0dab9e602a1564de199` |
| `48_wells/round3.tsv` | `35e18a9c0b63befe55320c25384232cecabbed7211a56454f6d77eef5e192b34` |
| `48_wells/crosswalk.tsv` | `6006ff1d97dd145f9bd6863bb8e0ca2a8fbf2105390940a41e6c37821d6c3e77` |

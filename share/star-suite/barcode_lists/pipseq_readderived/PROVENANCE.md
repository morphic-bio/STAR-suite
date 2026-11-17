# PIP tier lists derived from reads

The four lists for `--soloPipseqWhitelist`, in read order (8, 6, 6 and 8
bases; 96 sequences each). See `docs/PIPSEQ_RNA.md`.

## Source

Derived from public reads alone with `scripts/make_pipseq_tier_lists.py`: the
96 most frequent sequences at each tier position in the first 500,000 barcode
reads of GEO GSE291758, run SRR32666640, among reads in which the three
constant linkers match exactly at one phase (454,620 reads). No reference list
enters the derivation. The last listed sequence has 90, 246, 36 and 2.8 times
the reads of the next one in tiers 1 to 4.

Cross-check only: the cell labels of five public GEO `barcodes.tsv.gz` files
(GSM9890436 to GSM9890440) contain the same 96 sequences per tier. The read
counts alone separate the tier 4 list weakly (2.8 times); that the set is
complete for tier 4 rests on this cross-check.

## Scope

These are the sequences read off one public library. They are not a
manufacturer's list, and nothing here shows that they hold for another kit
version or chemistry. For another kit, or to check a library, derive the
lists from its own reads with `scripts/make_pipseq_tier_lists.py` and look at
the margins in its `summary.json`.

Derived 2026-11-15; details in `tests/ARTIFACTS.md` ("Newer PIP barcode tier
lists").

## SHA256

| File | SHA256 |
| --- | --- |
| `readderived_tier1.tsv` | `15c54f6fc2e0538c473c5c27fbb6584476f65912792d8b7fcbed629a96be38a5` |
| `readderived_tier2.tsv` | `6d7cce1aeae0f88ee624b53522cebdc5533e3b0ccdb3ce6675ff3a492858c3f4` |
| `readderived_tier3.tsv` | `7d1ff6846475bb5c61b7b4eeeae35bb66a247725519f1a68f66c334753ab75a0` |
| `readderived_tier4.tsv` | `c65bb07d2dbd80cd5b82cb294e3bbadc793587a4f2df93e66aa94ff780018514` |

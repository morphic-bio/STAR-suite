# Provenance — `catatac_crispri_guide_capture.csv`

This feature reference is the CAT-ATAC CRISPRi guide library (GSE288996, *Cell Rep
Methods* 2025) reproduced **verbatim from the study's deposited supplementary**
(`mmc2.xlsx`). It is intentionally **not** curated or de-duplicated.

## The expected `duplicate feature name` warning is by design

The deposited library lists the **HIC2** gene's two guides twice, with identical
sequences:

```
HIC2_1,TCGGCCCGCGACTCCTGTT,CRISPR Guide Capture   # lines 12-13
HIC2_2,CGGCTGTGAGCGGCGCTCG,CRISPR Guide Capture
...
HIC2_1,TCGGCCCGCGACTCCTGTT,CRISPR Guide Capture   # lines 56-57 (duplicate in source)
HIC2_2,CGGCTGTGAGCGGCGCTCG,CRISPR Guide Capture
```

This duplicate is present in the **source**, confirmed in every derived file
(`catatac_gse288996/guide_ref/ps_ref_crispri.tsv`,
`catatac_gse288996/guide_ref/catatac_guide_library_full.tsv`). STAR's
`process_features` ingests the file unmodified, keeps the first occurrence, drops the
later duplicate, and logs:

```
Warning: duplicate feature name 'HIC2_1' ... ; ignoring later definition and keeping the first
Warning: duplicate feature name 'HIC2_2' ... ; ignoring later definition and keeping the first
```

**Why we keep it this way.** Feeding the deposited library byte-for-byte and letting
STAR dedupe transparently (with a logged warning) means the pipeline provably uses the
**published library exactly as deposited** — there is no opportunity to be accused of
silently curating, re-ordering, or altering the guide set. The warning *is* the audit
trail. Editing this CSV to remove the duplicate would forfeit that guarantee.

**Effect on results: none.** The two HIC2 entries are byte-identical, so de-duplication
loses no guide. The library resolves to **54 unique guides** (27 target genes × 2
guides), which is the count the verifier asserts and the merged MEX carries.

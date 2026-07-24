# Hash-partitioned spatial GEX reconciliation

Date: 2026-07-24

Branch: `performance/gex-hash-reconcile-20260724`

Base: `c0c3d59a0662b869a34367137b6d5347e2aafd4b`

## Purpose

The spatial GEX resolver emits provisional support in feature order so it can
bound UMI-correction state by feature. Final `MultiGeneUMI_CR` reconciliation
instead groups by `(umi_mode, product, candidate, corrected_umi)`. The first
implementation required a global GNU text sort between those two stages. That
is correct but unnecessarily expensive for full-slide provisional support.

`molecule_first_gex_reconcile_hash` replaces that global sort with a two-pass,
bounded hash operation:

1. validate feature order and hash each complete support row by the final
   reconciliation key into a fixed number of private partitions;
2. load and reconcile one partition at a time, retaining feature encounter
   order within each key;
3. remove each partition immediately after it is processed.

The hash is fixed FNV-1a, output traversal is partition order followed by
first-seen key order, and no result depends on `unordered_map` iteration.
Changing the requested partition count can change row order, but not molecule
content. Repeating with the same inputs and options is byte deterministic.

The scientific policy is unchanged. Integer products call the existing
`MultiGeneUmiCr::resolve` implementation. Soft expected-count reconciliation
uses the same corrected-support winner, tie, and original-support dominance
rules as `molecule_first_gex_reconcile`.

## Interface

```text
molecule_first_gex_reconcile_hash \
  --input provisional/gex_provisional_support.tsv \
  --input-feature-sorted \
  --resolved-source provisional \
  --out-dir final \
  --tmp-dir tmp \
  --hash-partitions 256 \
  --partition-buffer-kb 1024
```

The input and output directories must be distinct. The output must be new or
empty. Temporary partitions are private to the process and are deleted after
successful processing unless `--keep-partitions` is explicitly supplied for
debugging.

## Regression gates

The tracked CLI regression constructs feature-sorted input that is
deliberately not final-key sorted. It checks:

- scientific content parity with the original streaming reconciler;
- byte-identical repeat output;
- `MultiGeneUMI_CR` accept, tie, and original-dominance cases;
- soft expected-count parity;
- row conservation and required hard-link ancestry;
- rejection of nonempty output and non-feature-sorted input;
- removal of processed temporary partitions.

On the existing ovarian 100K fixture, two 16-partition runs were byte
identical. Materialization of the hash result was byte-for-byte identical to
the previous sorted result for strict, soft expected-count, hard, and gated
hard products at 2, 8, and 16 micrometres. The 93,418,150-byte provisional
table contained 570,795 rows and the measured hash reconciliation wall time
was 1.34 seconds with 25,692 KiB peak RSS.

Fixture inputs and generated parity outputs remain untracked under:

```text
/mnt/pikachu/star-spatial/gex_sidecar_tests/
  20260722_ovarian_100k_bounded_gex_parity_v2/
  20260724_ovarian_100k_hash_parity_v1/
```

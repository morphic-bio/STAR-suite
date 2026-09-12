# Eight-byte FLEX RAM records with count overflow

Implemented as candidate R11, commit `7938428`, after the user selected an
explicit overflow mechanism. Clean build and ASan/UBSan boundary, merge and
RAM/spill tests passed; physical-prefix, L004 and full CBQ/BGZF validation passed. The earlier design
calculation below is retained as context, with the implementation decisions
stated here.

- Counts 0–62 are inline; 63 is an escape to a full 30-bit count stored in a
  sorted per-run `OverflowCount` vector. Normal records never search it.
- Background merges remap the exceptional records' positions into the new
  output run, preserving all counts and probe-region flags.
- The production collapse consumes compact runs directly. Reusable legacy
  readers still return full records. RAM-to-spill conversion reconstructs
  the unchanged 12-byte schema/checksum format.
- Buckets with more than 4,096 CB indices use the old 12-byte representation.
  Exceptional counts alone never cause truncation or a whole-run fallback.
- The logical 12-byte auto-spill accounting stays conservative and unchanged.
- The final CB-BUCKET log reports compact records, count-overflow records and
  encoded record bytes so the real dataset validates its unit-count claim.

For the current full 320K producer, every observation has count one before
UMI collapse, so no count overflow is expected. This follows from the writer
and background-merge contracts and the full CBQ and BGZF wrappers both measured zero overflow records.

| RAM word field | Bits |
| --- | ---: |
| Bucket-relative CB index | 12 |
| Tag, including reserved zero | 5 |
| Gene | 15 |
| UMI | 24 |
| Probe-region state | 2 |
| Count or escape | 6 |

The 256 contiguous buckets over 737,280 CBs contain 2,880 CB indices each.
For other sizes, bucket starts use ceiling division and are checked in tests
on a non-divisible whitelist. No UMI, tag or gene bit is discarded.

Full-set arithmetic: 6,876,635,217 records * 4 bytes saved = 27,506,540,868
bytes (25.617462 GiB) of record payload before overflow/run metadata. The full CBQ and BGZF runs both measured exactly this payload saving;
RSS is a separate whole-process measurement below.

Tests are in `tests/test_compact_cb_run.cpp`; executable, build logs, expected
fixture-error attempts and final sanitizer pass are under
`/mnt/pikachu/star_suite_paper/analysis/flex_ram_20260912/r11/`.

Pre-aggregating equal keys before storage would change the unit-count
contract and require measuring overflow frequency. The present implementation
preserves individual observations until the existing UMI-collapse stage.

Measured prototype encoding: 94,253,684 compact records, zero overflow
entries, 754,029,472 encoded bytes. Scientific outputs are exact against
v1.9.3, with 22,157 calls. The run took 22.973217 s and 6.165997 GiB peak
RSS using the ordinary allocator.

L004 also matched all scientific outputs: 149.405691 s, 33.225533 GiB peak
RSS, 249,194 cells. All 1,716,952,145 stored records were compact and none
overflowed. The payload reduction is measured, but L004 RSS is higher than
R8's 31.595810 GiB; actual full-set peak RSS must determine promotion.

Full CBQ acceptance passed: 463.982182 s / 75.284855 GiB peak RSS, exact
scientific outputs and read-classification totals, 333,439 cells. The actual
encoding is 6,876,635,217 compact records, zero overflow entries, and
55,013,081,736 payload bytes. This default-allocator build is below the
120 decimal GB target. Full BGZF also passed: 723.714585 s /
74.402870 GiB (79.889474 decimal GB), all 95
scientific digests and read-classification totals exact, 333,439 cells.
Both formats measured the same compact-record/overflow/payload totals.

# Provider-neutral partition manifest

Public scatter/gather recipes consume a partition set; they do not prescribe or
implement how that set is created. The producer may be a site tool, workflow
engine, object store, preprocessing job, or another authorized provider.

The `biodepot.partition_manifest/v1` JSON contract requires:

- one portable dataset and partition-set identity;
- paired inputs with stable, contiguous ordinals;
- half-open logical read-pair ranges beginning at zero;
- complete, non-overlapping logical coverage;
- optional SHA-256 digests for the original logical input and each mate;
- no provider, storage implementation, credential, or transport metadata.

Example:

```json
{
  "schema": "biodepot.partition_manifest/v1",
  "dataset_id": "example-pe",
  "partition_set_id": "example-pe-parts-v1",
  "source_digest": "0000000000000000000000000000000000000000000000000000000000000000",
  "partitions": [
    {
      "id": "part-000000",
      "ordinal": 0,
      "mate1": "inputs/part-000000_R1.fastq.gz",
      "mate2": "inputs/part-000000_R2.fastq.gz",
      "logical_range": {"unit": "read_pairs", "start": 0, "end": 100000}
    },
    {
      "id": "part-000001",
      "ordinal": 1,
      "mate1": "inputs/part-000001_R1.fastq.gz",
      "mate2": "inputs/part-000001_R2.fastq.gz",
      "logical_range": {"unit": "read_pairs", "start": 100000, "end": 200000}
    }
  ]
}
```

The validator rejects gaps, overlaps, reordered ordinals, unexpected fields,
duplicate identities, malformed digests, and missing local files when
`--check-files` is selected. A scheduler must gather only after each required
ordinal has terminal success evidence. Failed or retried tasks retain the same
ordinal and locked manifest.

The contract deliberately says nothing about a partition provider's internal
representation. Provider publication and licensing are independent of this
open consumer interface.

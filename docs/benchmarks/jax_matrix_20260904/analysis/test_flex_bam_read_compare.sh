#!/usr/bin/env bash
set -euo pipefail

tool="${1:-./flex_bam_read_compare}"
tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT

python3 - "${tmp}" <<'PY'
import struct
import sys
from pathlib import Path

root = Path(sys.argv[1])
tags = {"BC001": "ACTTTAGG", "BC002": "AACGGGAA"}
(root / "tags.tsv").write_text("".join(f"{key}\t{value}\n" for key, value in tags.items()))

# The same CB16 occurs under both tags. They are deliberately different cells.
cb = "AAAAAAAAAAAAAAAA"
(root / "star_bc1.tsv").write_text(cb + "\n")
(root / "star_bc2.tsv").write_text(cb + "\n")
(root / "cr.tsv").write_text(
    cb + tags["BC001"] + "-1\n" +
    "CCCCCCCCCCCCCCCC" + tags["BC001"] + "-1\n" +
    "GGGGGGGGGGGGGGGG" + tags["BC002"] + "-1\n"
)
(root / "genes.txt").write_text("ENSG000001\nENSG000002\n")

def pack(seq):
    value = 0
    for base in seq:
        value = (value << 2) | "ACGT".index(base)
    return value & ((1 << 64) - 1), value >> 64

seq_h0 = "A" * 50
seq_h1 = "C" * 50
seq_deny = "G" * 50
# Cache v3 stores probe region in resolved_gene bits 30-31.
spliced = 1 << 30
unspliced = 2 << 30
records = [
    (*pack(seq_h0), spliced | 1, 0, 0, 1),
    # Identical probe sequence resolves differently for the second sample.
    (*pack(seq_h0), unspliced | 2, 0, 0, 2),
    (*pack(seq_h1), unspliced | 2, 1, 0, 0),
    (*pack(seq_deny), 0, 2, 1, 0),
]
records.sort(key=lambda row: (row[1], row[0], row[5]))
with (root / "cache.bin").open("wb") as out:
    out.write(struct.pack("<8sHHIQ", b"FH01SEQ1", 3, 50, 24, len(records)))
    for lo, hi, gene, cls, negative, sample in records:
        out.write(struct.pack("<QQIBBH", lo, hi, gene, cls, negative, sample))

qual = "I" * 50
def sam(qname, flag, cb24, seq, gx, xf, cigar="50M", mapq=255,
        ub="ACGTACGTACGT", ur=None):
    if ur is None:
        ur = ub
    tags = [f"CB:Z:{cb24}-1", "CR:Z:" + cb24[:16],
            "CY:Z:" + "I" * 16, "UY:Z:" + "I" * 12,
            f"fx:Z:{gx}", f"GX:Z:{gx}", f"GN:Z:fixture_{gx}", f"xf:i:{xf}",
            "pr:Z:fixture_probe", "RE:A:E"]
    if ub is not None:
        tags.append(f"UB:Z:{ub}")
    if ur is not None:
        tags.append(f"UR:Z:{ur}")
    return "\t".join(map(str, [qname, flag, "ENSG000001|probe", 1, mapq, cigar,
                                "*", 0, 0, seq, qual] + tags))

rows = [
    # shared: H0 and same gene
    sam("read_shared_h0", 0, cb + tags["BC001"], seq_h0, "ENSG000001", 8),
    # STAR-only: H1, but CR reports another gene
    sam("read_star_only_conflict", 0, cb + tags["BC002"], seq_h1, "ENSG000001", 8,
        "48M2S", 30, "CCCCCCCCCCCC"),
    # CR-only: cache DENY
    sam("read_cr_only_deny", 0, "CCCCCCCCCCCCCCCC" + tags["BC001"], seq_deny,
        "ENSG000002", 8, "25M1I24M", 20, "GGGGGGGGGGGG"),
    # CR-only under the second tag: cache miss
    sam("read_cr_only_miss", 0, "GGGGGGGGGGGGGGGG" + tags["BC002"], "T" * 50,
        "ENSG000002", 0, "*", 0, "TTTTTTTTTTTT"),
    # SAM stores reverse-strand SEQ reference-oriented; restoring original R2
    # must recover the all-A H0 window.
    sam("read_reverse_h0", 0x10, cb + tags["BC001"], "T" * 50, "ENSG000001", 8),
    # The first 50 bases miss; the physical +1 window is the H1 probe.
    sam("read_offset_h1", 0, cb + tags["BC002"], "A" + seq_h1,
        "ENSG000002", 8, "51M", 60, "AAAAACCCCCCC"),
    # Same cache key as BC001, but sample-aware lookup selects BC002's gene.
    sam("read_sample_h0", 0, cb + tags["BC002"], seq_h0, "ENSG000002", 8,
        "TTTTTCCCCCCC"),
    # The representative alone misses the cache, but another primary read for
    # the exact CB24+UB+fx molecule has a same-gene H0 hit.
    sam("read_rescue_representative", 0, cb + tags["BC001"], "T" * 50,
        "ENSG000001", 8, "50M", 255, "CATCATCATCAT"),
    sam("read_rescue_duplicate", 0x400, cb + tags["BC001"], seq_h0,
        "ENSG000001", 0, "50M", 255, "CATCATCATCAT"),
    # UR is accepted only as an explicit, measured fallback when UB is absent.
    sam("read_ur_fallback", 0, "GGGGGGGGGGGGGGGG" + tags["BC002"], "T" * 50,
        "ENSG000002", 8, "50M", 20, None, "GATTACAGATTA"),
    # Supplementary must be excluded by default.
    sam("read_secondary", 0x800, cb + tags["BC001"], seq_h0, "ENSG000001", 8),
]
# Exercise the flat molecule table through at least one rehash without changing
# the counted-molecule expectations.
for value in range(1, 801):
    number = value
    bases = []
    for _ in range(16):
        bases.append("ACGT"[number & 3])
        number >>= 2
    extra_cb = "".join(reversed(bases))
    rows.append(sam(f"rehash_{value}", 0, extra_cb + tags["BC001"], "T" * 50,
                    "ENSG000002", 0, "50M", 20, "ACACACACACAC"))
(root / "fixture.sam").write_text("@HD\tVN:1.6\tSO:unsorted\n" + "\n".join(rows) + "\n")
PY

"${tool}" \
  --input "${tmp}/fixture.sam" --input-kind sam \
  --out-prefix "${tmp}/out" --tag-map "${tmp}/tags.tsv" \
  --star-cells "BC001=${tmp}/star_bc1.tsv" \
  --star-cells "BC002=${tmp}/star_bc2.tsv" \
  --cr-cells "${tmp}/cr.tsv" \
  --hash-cache "${tmp}/cache.bin" --gene-list "${tmp}/genes.txt" \
  --details all --detail-max 0

python3 - "${tmp}" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
metrics = dict(csv.reader((root / "out.metrics.tsv").open(), delimiter="\t"))
assert metrics["star_cells"] == "2"
assert metrics["cr_cells"] == "3"
assert metrics["shared_cells"] == "1"
assert metrics["star_only_cells"] == "1"
assert metrics["cr_only_cells"] == "2"
assert metrics["included_primary_records"] == "810"
assert metrics["skipped_secondary_or_supplementary"] == "1"
assert metrics["cr_xf8_primary_records"] == "8"
assert metrics["unique_counted_molecule_keys"] == "7"
assert metrics["counted_molecule_keys_multiple_xf8"] == "1"
assert metrics["counted_molecule_keys_ur_fallback"] == "1"
assert metrics["molecule_missing_ub_records"] == "1"
assert metrics["molecule_ur_fallback_records"] == "1"
assert metrics["molecule_outcome_rescued_same_gene"] == "1"
assert metrics["molecule_keys_all_valid_records"] == "808"

rows = {row["qname"]: row for row in csv.DictReader((root / "out.details.tsv").open(), delimiter="\t")}
assert rows["read_shared_h0"]["cell_key"] == "AAAAAAAAAAAAAAAA|BC001"
assert rows["read_shared_h0"]["cell_class"] == "shared"
assert rows["read_shared_h0"]["hash_verdict"] == "KEEP_H0"
assert rows["read_shared_h0"]["gene_relation"] == "same"
assert rows["read_star_only_conflict"]["cell_key"] == "AAAAAAAAAAAAAAAA|BC002"
assert rows["read_star_only_conflict"]["cell_class"] == "star_only"
assert rows["read_star_only_conflict"]["hash_verdict"] == "KEEP_H1"
assert rows["read_star_only_conflict"]["gene_relation"] == "conflict"
assert rows["read_cr_only_deny"]["cell_class"] == "cr_only"
assert rows["read_cr_only_deny"]["hash_verdict"] == "DENY_CACHE"
assert rows["read_cr_only_miss"]["cell_key"] == "GGGGGGGGGGGGGGGG|BC002"
assert rows["read_cr_only_miss"]["hash_verdict"] == "MISS"
assert rows["read_reverse_h0"]["hash_verdict"] == "KEEP_H0"
assert rows["read_reverse_h0"]["probe_region"] == "spliced"
assert rows["read_offset_h1"]["hash_verdict"] == "KEEP_H1"
assert rows["read_offset_h1"]["hash_offset"] == "1"
assert rows["read_offset_h1"]["probe_region"] == "unspliced"
assert rows["read_sample_h0"]["hash_verdict"] == "KEEP_H0"
assert rows["read_sample_h0"]["hash_gene"] == "ENSG000002"
assert rows["read_sample_h0"]["sample_index"] == "2"
assert rows["read_rescue_representative"]["hash_verdict"] == "MISS"
assert rows["read_rescue_duplicate"]["hash_verdict"] == "KEEP_H0"
assert rows["read_rescue_representative"]["mapq_bin"] == "255/unavailable"

regions = list(csv.DictReader((root / "out.regions.tsv").open(), delimiter="\t"))
assert any(row["cr_region"] == "E" and row["probe_region"] == "spliced"
           and row["hash_gene"] == "ENSG000001" for row in regions)
assert (root / "out.counted_records.tsv").exists()
assert not (root / "out.counted_molecules.tsv").exists()

molecules = list(csv.DictReader((root / "out.molecules.tsv").open(), delimiter="\t"))
rescued = [row for row in molecules if row["outcome"] == "rescued_same_gene"]
assert len(rescued) == 1
assert rescued[0]["quant_gene_ids"] == "ENSG000001"
assert rescued[0]["hash_verdict"] == "KEEP_H0"
assert rescued[0]["hash_gene"] == "ENSG000001"
assert rescued[0]["molecules"] == "1", rescued
assert rescued[0]["input_records"] == "2", rescued
assert rescued[0]["xf8_records"] == "1", rescued

cigars = list(csv.DictReader((root / "out.cigars.tsv").open(), delimiter="\t"))
assert any(row["mapq_bin"] == "255/unavailable" for row in cigars)

# Regression against concordance_vs_cr.py's corrected composite identity helper.
sys.path.insert(0, str(Path.cwd().parent))
import concordance_vs_cr as cc
tag_map = cc.read_tag_map(root / "tags.tsv")
assert cc.cell_key("AAAAAAAAAAAAAAAAACTTTAGG-1", tag_map=tag_map) == rows["read_shared_h0"]["cell_key"]
assert cc.cell_key("AAAAAAAAAAAAAAAAAACGGGAA-1", tag_map=tag_map) == rows["read_star_only_conflict"]["cell_key"]
assert rows["read_shared_h0"]["cell_key"] != rows["read_star_only_conflict"]["cell_key"]
PY

echo "PASS: native Flex BAM evidence analyzer fixture"

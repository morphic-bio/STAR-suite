#!/usr/bin/env bash
# Materialize long-format HIV state count table for YW8 from Wei-et-al-2023 bc_* files.
set -euo pipefail

OUT_DIR="${1:-/tmp/hiv_dogma_gse239916/star_four_arm_downsample}"
PARTICIPANT="${2:-YW8}"
DNA_BC="${HIV_DNA_BC:-/tmp/hiv_dogma_gse239916/refs/bc_ATAC_HIV.txt}"
RNA_BC="${HIV_RNA_BC:-/tmp/hiv_dogma_gse239916/refs/bc_RNA_HIV.txt}"

mkdir -p "$OUT_DIR"
TABLE="$OUT_DIR/hiv_state_counts.tsv"

python3 - "$OUT_DIR" "$PARTICIPANT" "$DNA_BC" "$RNA_BC" "$TABLE" <<'PY'
import sys
from pathlib import Path

out_dir, participant, dna_path, rna_path, table_path = sys.argv[1:6]

def parse_counts(path: Path, participant: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            return counts
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            cell = parts[0]
            if not cell.startswith(f"{participant}_"):
                continue
            barcode = cell.split("_", 1)[1]
            counts[barcode] = int(parts[2])
    return counts

dna = parse_counts(Path(dna_path), participant)
rna = parse_counts(Path(rna_path), participant)
rows = ["barcode\tfeature_id\tcount"]
for barcode, count in sorted(dna.items()):
    rows.append(f"{barcode}\tHIV_DNA\t{count}")
for barcode, count in sorted(rna.items()):
    rows.append(f"{barcode}\tHIV_RNA\t{count}")
Path(table_path).write_text("\n".join(rows) + "\n", encoding="utf-8")
print(f"Wrote {len(rows)-1} rows to {table_path}")
PY

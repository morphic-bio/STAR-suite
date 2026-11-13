#!/usr/bin/env python3
"""Write the Visium HD barcode inputs that STARgrid needs from a slide-design file.

The slide-design file is distributed with Space Ranger (for example
lib/python/cellranger/barcodes/visium_hd_v1.slide); it is not redistributed
here. The file is read with the public protobuf wire format only: the script
collects the repeated fields whose values are all DNA sequences and expects
exactly two lists of equal length, the BC1 sequences (one per column) and the
BC2 sequences (one per row). BC1 is the list with the longer sequences: its
shortest and its longest sequences are each longer than BC2's.

Outputs, all in --out-dir:
  bc1_full_oligos.txt  BC1 sequences in column order (--soloSpatialBc1Oligos)
  bc2_full_oligos.txt  BC2 sequences in row order (--soloSpatialBc2Oligos)
  barcode_coords.tsv   every capture location in row-major order
                       (--soloSpatialBarcodeContract is the directory)
  summary.json         sizes and checksums
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path


GRID_SIDE = 3350
DNA = re.compile(rb"^[ACGTacgt]+$")
# Checksums of the oligo lists from visium_hd_v1.slide (identical in Space
# Ranger 3.0.0, 4.0.1 and 4.1.0), used for the results in the STARgrid paper.
KNOWN_V1 = {
    "bc1_full_oligos.txt": "24dfc754a706239d990ef8037c78f39bfe9c745b4ad424aa64dfc2317c4857f8",
    "bc2_full_oligos.txt": "2971bbbba0abdd7457b51a5646845df6ec2e3d620248dab17182fa69c660a1bb",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write STARgrid's Visium HD barcode inputs from a Space Ranger .slide file."
    )
    parser.add_argument("--slide", required=True, type=Path, help="Space Ranger Visium HD .slide file")
    parser.add_argument("--out-dir", required=True, type=Path, help="Output directory")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing outputs")
    return parser.parse_args(argv)


def read_varint(data: bytes, offset: int) -> tuple[int, int]:
    value = 0
    shift = 0
    while offset < len(data) and shift < 64:
        byte = data[offset]
        offset += 1
        value |= (byte & 0x7F) << shift
        if byte < 0x80:
            return value, offset
        shift += 7
    raise ValueError("malformed varint")


def dna_lists(data: bytes, path: tuple[int, ...] = (), depth: int = 0) -> dict[tuple[int, ...], list[str]]:
    """Return the DNA strings of a protobuf message by field path; raise if it does not parse."""
    found: dict[tuple[int, ...], list[str]] = {}
    offset = 0
    while offset < len(data):
        key, offset = read_varint(data, offset)
        field, wire = key >> 3, key & 7
        if field == 0:
            raise ValueError("field number 0")
        if wire == 0:
            _, offset = read_varint(data, offset)
        elif wire == 1:
            offset += 8
        elif wire == 5:
            offset += 4
        elif wire == 2:
            size, offset = read_varint(data, offset)
            value = data[offset : offset + size]
            offset += size
            if offset > len(data):
                raise ValueError("field runs past the end")
            if value and DNA.match(value):
                found.setdefault(path + (field,), []).append(value.decode("ascii").upper())
            elif depth < 8 and len(value) > 1:
                try:
                    nested = dna_lists(value, path + (field,), depth + 1)
                except ValueError:
                    continue  # a string or bytes field, not a message
                for key_path, values in nested.items():
                    found.setdefault(key_path, []).extend(values)
        else:
            raise ValueError(f"unsupported wire type {wire}")
        if offset > len(data):
            raise ValueError("field runs past the end")
    return found


def read_oligos(path: Path) -> tuple[list[str], list[str]]:
    try:
        found = dna_lists(path.read_bytes())
    except ValueError as exc:
        raise SystemExit(f"ERROR: {path} is not a protobuf slide-design file ({exc})") from exc
    lists = [values for values in found.values() if len(values) > 1]
    if len(lists) != 2 or len(lists[0]) != len(lists[1]):
        sizes = sorted(len(values) for values in lists)
        raise SystemExit(f"ERROR: expected two lists of barcode sequences of equal length; found {sizes}")
    for values in lists:
        if len(set(values)) != len(values):
            raise SystemExit("ERROR: a barcode list repeats a sequence")
    lengths = [(min(map(len, values)), max(map(len, values))) for values in lists]
    if lengths[0][0] > lengths[1][0] and lengths[0][1] > lengths[1][1]:
        bc1, bc2 = lists
    elif lengths[1][0] > lengths[0][0] and lengths[1][1] > lengths[0][1]:
        bc2, bc1 = lists
    else:
        raise SystemExit(f"ERROR: cannot tell BC1 from BC2 by length: {lengths}")
    if len(bc1) != GRID_SIDE:
        raise SystemExit(f"ERROR: STARgrid expects a {GRID_SIDE} x {GRID_SIDE} grid; this design has {len(bc1)}")
    return bc1, bc2


def output_path(out_dir: Path, name: str, overwrite: bool) -> Path:
    path = out_dir / name
    if path.exists() and not overwrite:
        raise SystemExit(f"ERROR: output exists, use --overwrite: {path}")
    return path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.slide.is_file():
        raise SystemExit(f"ERROR: slide-design file does not exist: {args.slide}")
    bc1, bc2 = read_oligos(args.slide)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    names = ("bc1_full_oligos.txt", "bc2_full_oligos.txt", "barcode_coords.tsv", "summary.json")
    paths = {name: output_path(args.out_dir, name, args.overwrite) for name in names}

    for name, values in (("bc1_full_oligos.txt", bc1), ("bc2_full_oligos.txt", bc2)):
        paths[name].write_text("".join(f"{value}\n" for value in values), encoding="ascii")
    with paths["barcode_coords.tsv"].open("w", encoding="ascii", buffering=1 << 20) as coords:
        coords.write("barcode\trow2\tcol2\n")
        for row, row_oligo in enumerate(bc2):
            coords.write("".join(f"{col_oligo}_{row_oligo}\t{row}\t{col}\n" for col, col_oligo in enumerate(bc1)))

    checksums = {name: sha256(paths[name]) for name in KNOWN_V1}
    summary = {
        "slide": str(args.slide),
        "slide_sha256": sha256(args.slide),
        "grid_rows": len(bc2),
        "grid_cols": len(bc1),
        "bc1_lengths": sorted({len(value) for value in bc1}),
        "bc2_lengths": sorted({len(value) for value in bc2}),
        "oligo_sha256": checksums,
        "matches_visium_hd_v1": checksums == KNOWN_V1,
    }
    paths["summary.json"].write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    for path in paths.values():
        print(f"wrote {path}")
    if not summary["matches_visium_hd_v1"]:
        print("WARNING: the oligo lists differ from those of visium_hd_v1.slide", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

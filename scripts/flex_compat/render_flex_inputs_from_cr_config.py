#!/usr/bin/env python3
import argparse
import csv
import shlex
import sys
from collections import OrderedDict
from pathlib import Path


def resolve_path(raw_value: str, base_dir: Path) -> str:
    value = (raw_value or "").strip()
    if not value:
        return ""
    path = Path(value)
    if not path.is_absolute():
        path = base_dir / path
    return str(path.resolve())


def parse_cr_config(path: Path) -> dict:
    sections = {}
    current = None
    with path.open(newline="") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                current = line[1:-1].strip().lower()
                sections[current] = []
                continue
            if current is None:
                continue
            sections[current].append(raw)

    config_dir = path.resolve().parent
    parsed = {"gene-expression": {}, "libraries": [], "samples": []}
    if "gene-expression" in sections:
        reader = csv.reader(sections["gene-expression"])
        for row in reader:
            if len(row) >= 2:
                key = row[0].strip().lower()
                value = row[1].strip()
                if key in {"reference", "probe-set", "probeset"}:
                    value = resolve_path(value, config_dir)
                parsed["gene-expression"][key] = value

    if "libraries" not in sections:
        raise SystemExit(f"Missing [libraries] section in {path}")
    reader = csv.DictReader(sections["libraries"])
    for row in reader:
        cleaned = {k.strip().lower(): (v.strip() if v is not None else "") for k, v in row.items()}
        if cleaned.get("fastqs"):
            cleaned["fastqs"] = resolve_path(cleaned["fastqs"], config_dir)
            parsed["libraries"].append(cleaned)

    if "samples" not in sections:
        raise SystemExit(f"Missing [samples] section in {path}")
    reader = csv.DictReader(sections["samples"])
    for row in reader:
        cleaned = {k.strip().lower(): (v.strip() if v is not None else "") for k, v in row.items()}
        if cleaned.get("sample_id"):
            parsed["samples"].append(cleaned)

    if not parsed["libraries"]:
        raise SystemExit(f"No library rows found in {path}")
    if not parsed["samples"]:
        raise SystemExit(f"No sample rows found in {path}")
    return parsed


def normalize_feature_type(value: str) -> str:
    return "".join(ch.lower() for ch in value if ch.isalnum())


def is_gex_library(row: dict) -> bool:
    feature_type = row.get("feature_types") or row.get("feature_type") or row.get("library_type") or ""
    norm = normalize_feature_type(feature_type)
    return "geneexpression" in norm or norm == "gex"


def enumerate_pairs(fastq_dir: Path) -> tuple[list[str], list[str]]:
    r1s = sorted(str(p.resolve()) for p in fastq_dir.glob("*_R1_001.fastq.gz"))
    if not r1s:
        raise SystemExit(f"No R1 FASTQs found under {fastq_dir}")
    r2s = []
    for r1 in r1s:
        r2 = r1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
        if not Path(r2).is_file():
            raise SystemExit(f"Missing R2 pair for {r1}")
        r2s.append(r2)
    return r1s, r2s


def split_probe_ids(value: str) -> list[str]:
    raw = value.replace(";", ",").replace("|", ",").replace(" ", ",")
    return [item for item in raw.split(",") if item]


def load_probe_catalog(path: Path) -> dict[str, list[tuple[str, str]]]:
    mapping: dict[str, list[tuple[str, str]]] = {}
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise SystemExit(f"Invalid sample probe catalog row in {path}: {raw.rstrip()}" )
            variant, canonical, sample_id = parts[0], parts[1], parts[2]
            mapping.setdefault(sample_id, []).append((variant, canonical))
    return mapping


def write_probe_list(probe_set_csv: Path, out_path: Path) -> None:
    genes = OrderedDict()
    header = None
    with probe_set_csv.open(newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            if header is None:
                header = [cell.strip().lower() for cell in row]
                continue
            if header is None:
                continue
            rec = {header[i]: row[i].strip() if i < len(row) else "" for i in range(len(header))}
            included = rec.get("included", "TRUE").strip().upper()
            if included not in {"", "TRUE", "YES", "1"}:
                continue
            gene_id = rec.get("gene_id", "")
            if gene_id:
                genes.setdefault(gene_id, None)
    if not genes:
        raise SystemExit(f"No included gene_id values found in probe set {probe_set_csv}")
    with out_path.open("w", newline="") as handle:
        for gene_id in genes.keys():
            handle.write(f"{gene_id}\n")


def emit_env(values: dict[str, str]) -> None:
    for key, value in values.items():
        print(f"{key}={shlex.quote(value)}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Render STAR FLEX inputs from a Cell Ranger config.csv")
    parser.add_argument("--config", required=True, help="Cell Ranger Flex config.csv")
    parser.add_argument("--sample-whitelist-out", required=True, help="Output sample_whitelist.tsv path")
    parser.add_argument("--sample-probes-out", required=True, help="Output sample probes TSV path")
    parser.add_argument("--probe-list-out", required=True, help="Output probe_list.txt path")
    parser.add_argument("--probe-catalog", required=True, help="Catalog of sample probe barcodes (variant, canonical, sample_id)")
    parser.add_argument("--emit-env", action="store_true", help="Emit shell-compatible KEY=VALUE lines")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    parsed = parse_cr_config(config_path)
    gex_libs = [row for row in parsed["libraries"] if is_gex_library(row)]
    if not gex_libs:
        raise SystemExit(f"No Gene Expression library rows found in {config_path}")

    probe_set = parsed.get("gene-expression", {}).get("probe-set") or parsed.get("gene-expression", {}).get("probeset")
    if not probe_set:
        raise SystemExit(f"No [gene-expression] probe-set found in {config_path}")
    probe_set_path = Path(probe_set)
    if not probe_set_path.is_file():
        raise SystemExit(f"Probe set file does not exist: {probe_set_path}")

    sample_whitelist_out = Path(args.sample_whitelist_out)
    sample_probes_out = Path(args.sample_probes_out)
    probe_list_out = Path(args.probe_list_out)
    for out_path in (sample_whitelist_out, sample_probes_out, probe_list_out):
        out_path.parent.mkdir(parents=True, exist_ok=True)

    gex_dirs = []
    gex_r1 = []
    gex_r2 = []
    seen_gex_dirs = set()
    for row in gex_libs:
        fastq_dir = row["fastqs"]
        if fastq_dir in seen_gex_dirs:
            continue
        seen_gex_dirs.add(fastq_dir)
        gex_dirs.append(fastq_dir)
        r1s, r2s = enumerate_pairs(Path(fastq_dir))
        gex_r1.extend(r1s)
        gex_r2.extend(r2s)

    catalog = load_probe_catalog(Path(args.probe_catalog).resolve())
    whitelist_rows: list[tuple[str, str]] = []
    probe_rows: list[tuple[str, str, str]] = []
    for sample in parsed["samples"]:
        label = sample["sample_id"]
        probe_ids = split_probe_ids(sample.get("probe_barcode_ids", ""))
        if not probe_ids:
            raise SystemExit(f"Sample {label} has no probe_barcode_ids in {config_path}")
        for probe_id in probe_ids:
            entries = catalog.get(probe_id)
            if not entries:
                raise SystemExit(f"probe_barcode_id {probe_id} for sample {label} not found in {args.probe_catalog}")
            canonical_values = {canonical for _, canonical in entries}
            if len(canonical_values) != 1:
                raise SystemExit(f"probe_barcode_id {probe_id} maps to multiple canonical tags in {args.probe_catalog}")
            canonical = next(iter(canonical_values))
            whitelist_rows.append((label, canonical))
            for variant, canonical in entries:
                probe_rows.append((variant, canonical, label))

    whitelist_rows = list(OrderedDict(((label, canonical), None) for label, canonical in whitelist_rows).keys())
    probe_rows = list(OrderedDict(((variant, canonical, label), None) for variant, canonical, label in probe_rows).keys())

    with sample_whitelist_out.open("w", newline="") as handle:
        for label, canonical in whitelist_rows:
            handle.write(f"{label}\t{canonical}\n")

    with sample_probes_out.open("w", newline="") as handle:
        for variant, canonical, label in probe_rows:
            handle.write(f"{variant}\t{canonical}\t{label}\n")

    write_probe_list(probe_set_path, probe_list_out)

    if args.emit_env:
        emit_env(
            {
                "CR_CONFIG": str(config_path),
                "CR_GENE_EXPRESSION_REFERENCE": parsed.get("gene-expression", {}).get("reference", ""),
                "CR_GENE_EXPRESSION_PROBE_SET": str(probe_set_path),
                "GEX_FASTQ_DIRS": ",".join(gex_dirs),
                "GEX_R1": ",".join(gex_r1),
                "GEX_R2": ",".join(gex_r2),
                "FLEX_SAMPLE_WHITELIST": str(sample_whitelist_out.resolve()),
                "FLEX_ALLOWED_TAGS": str(sample_whitelist_out.resolve()),
                "FLEX_SAMPLE_PROBES": str(sample_probes_out.resolve()),
                "FLEX_PROBE_LIST": str(probe_list_out.resolve()),
                "FLEX_SAMPLE_IDS": ",".join(sample["sample_id"] for sample in parsed["samples"]),
            }
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())

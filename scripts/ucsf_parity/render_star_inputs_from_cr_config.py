#!/usr/bin/env python3
import argparse
import csv
import shlex
import sys
from pathlib import Path
from typing import Iterable


def normalize_feature_type(value: str) -> str:
    return "".join(ch.lower() for ch in value if ch.isalnum())


def resolve_path(raw_value: str, base_dir: Path) -> str:
    value = (raw_value or "").strip()
    if not value:
        return ""
    path = Path(value)
    if not path.is_absolute():
        path = base_dir / path
    return str(path.resolve())


def normalize_row(row: dict, config_dir: Path) -> dict:
    cleaned = {k.strip().lower(): (v.strip() if v is not None else "") for k, v in row.items()}
    if cleaned.get("fastqs"):
        cleaned["fastqs"] = resolve_path(cleaned["fastqs"], config_dir)
    if cleaned.get("star_feature_ref"):
        cleaned["star_feature_ref"] = resolve_path(cleaned["star_feature_ref"], config_dir)
    if cleaned.get("starfeatureref"):
        cleaned["starfeatureref"] = resolve_path(cleaned["starfeatureref"], config_dir)
    return cleaned


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

    parsed = {"gene-expression": {}, "feature": {}, "libraries": []}
    config_dir = path.resolve().parent

    for section in ("gene-expression", "feature", "reference"):
        if section in sections:
            reader = csv.reader(sections[section])
            for row in reader:
                if len(row) >= 2:
                    key = row[0].strip().lower()
                    value = row[1].strip()
                    if key in {"reference", "ref", "fastqs", "probe-set", "probeset"}:
                        value = resolve_path(value, config_dir)
                    parsed.setdefault(section, {})[key] = value

    if "libraries" not in sections:
        raise SystemExit(f"Missing [libraries] section in {path}")

    reader = csv.DictReader(sections["libraries"])
    for row in reader:
        cleaned = normalize_row(row, config_dir)
        if cleaned.get("fastqs"):
            parsed["libraries"].append(cleaned)

    if not parsed["libraries"]:
        raise SystemExit(f"No library rows found in {path}")

    return parsed


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


def row_value(row: dict, *keys: str) -> str:
    for key in keys:
        value = row.get(key, "")
        if value:
            return value
    return ""


def is_gex_library(row: dict) -> bool:
    feature_type = row_value(row, "feature_types", "feature_type", "library_type", "librarytype")
    norm = normalize_feature_type(feature_type)
    return "geneexpression" in norm or norm == "gex"


def is_crispr_library(row: dict) -> bool:
    feature_type = row_value(row, "feature_types", "feature_type", "library_type", "librarytype")
    norm = normalize_feature_type(feature_type)
    return "crisprguidecapture" in norm or "guidecapture" in norm or norm == "crispr"


def classify_libraries(libraries: list[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    gex = []
    features = []
    guides = []
    for row in libraries:
        if is_gex_library(row):
            gex.append(row)
        else:
            features.append(row)
            if is_crispr_library(row):
                guides.append(row)
    return gex, features, guides


def determine_sample_id(libraries: list[dict], explicit: str | None) -> str:
    if explicit:
        return explicit
    ids = {
        row.get("sample") or row.get("fastq_id") or row.get("library_id") or ""
        for row in libraries
    }
    ids.discard("")
    if len(ids) == 1:
        return next(iter(ids))
    raise SystemExit(
        "Could not infer a single sample ID from Cell Ranger config; pass --sample-id explicitly"
    )


def sanitize_csv_field(value: str) -> str:
    return value.replace("\n", " ").replace("\r", " ")


def pf_row_for_library(row: dict, sample_id: str) -> list[str]:
    feature_type = row_value(row, "feature_types", "feature_type")
    library_type = row_value(row, "library_type", "librarytype") or feature_type or (
        "Gene Expression" if is_gex_library(row) else feature_type
    )
    return [
        row["fastqs"],
        sample_id,
        library_type,
        feature_type or library_type,
        row_value(row, "star_chemistry", "starchemistry"),
        row_value(row, "star_feature_ref", "starfeatureref"),
        row_value(row, "star_library_id", "starlibraryid"),
    ]


def write_pf_multi_config(path: Path, sample_id: str, libraries: list[dict], global_feature_ref: str) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["[libraries]"])
        writer.writerow([
            "fastqs",
            "sample",
            "library_type",
            "feature_types",
            "star_chemistry",
            "star_feature_ref",
            "star_library_id",
        ])
        for row in libraries:
            writer.writerow(pf_row_for_library(row, sample_id))
        if global_feature_ref:
            writer.writerow([])
            writer.writerow(["[feature]"])
            writer.writerow(["ref", global_feature_ref])


def split_csv_values(value: str) -> list[str]:
    return [item for item in value.split(",") if item]


def apply_fastq_rewrites(libraries: list[dict], rewrite_map: dict[str, str]) -> list[dict]:
    rewritten = []
    for row in libraries:
        new_row = dict(row)
        fastq_dir = row.get("fastqs", "")
        if fastq_dir in rewrite_map:
            new_row["fastqs"] = rewrite_map[fastq_dir]
        rewritten.append(new_row)
    return rewritten


def collect_pairs_for_libraries(libraries: Iterable[dict]) -> tuple[list[str], list[str], list[str]]:
    seen_dirs = set()
    dirs = []
    r1_all: list[str] = []
    r2_all: list[str] = []
    for row in libraries:
        fastq_dir = row["fastqs"]
        if fastq_dir in seen_dirs:
            continue
        seen_dirs.add(fastq_dir)
        dirs.append(fastq_dir)
        r1s, r2s = enumerate_pairs(Path(fastq_dir))
        r1_all.extend(r1s)
        r2_all.extend(r2s)
    return dirs, r1_all, r2_all


def emit_env(values: dict[str, str]) -> None:
    for key, value in values.items():
        print(f"{key}={shlex.quote(value)}")


def parse_rewrite_args(entries: list[str]) -> dict[str, str]:
    rewrites = {}
    for entry in entries:
        if "=" not in entry:
            raise SystemExit(f"Invalid --rewrite-fastq-dir value '{entry}' (expected OLD=NEW)")
        old, new = entry.split("=", 1)
        rewrites[old] = new
    return rewrites


def main() -> int:
    parser = argparse.ArgumentParser(description="Render STAR perturb inputs from a Cell Ranger config.csv")
    parser.add_argument("--config", required=True, help="Cell Ranger config.csv")
    parser.add_argument("--pf-multi-out", required=True, help="Output pf_multi_config.csv path")
    parser.add_argument("--sample-id", help="Override sample ID for generated pf_multi_config.csv")
    parser.add_argument(
        "--rewrite-fastq-dir",
        action="append",
        default=[],
        help="Rewrite a library fastq directory when rendering outputs (OLD=NEW); may be repeated",
    )
    parser.add_argument("--emit-env", action="store_true", help="Emit shell-compatible KEY=VALUE lines")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    parsed = parse_cr_config(config_path)
    gex_libs, feature_libs, guide_libs = classify_libraries(parsed["libraries"])

    if not gex_libs:
        raise SystemExit(f"No Gene Expression library rows found in {config_path}")
    if not feature_libs:
        raise SystemExit(f"No non-GEX feature library rows found in {config_path}")

    sample_id = determine_sample_id(parsed["libraries"], args.sample_id)
    global_feature_ref = parsed.get("feature", {}).get("reference") or parsed.get("feature", {}).get("ref", "")
    for row in feature_libs:
        if not row_value(row, "star_feature_ref", "starfeatureref") and not global_feature_ref:
            raise SystemExit(
                f"Feature library '{row.get('fastqs', '')}' is missing star_feature_ref and [feature] ref is absent in {config_path}"
            )

    rewrites = parse_rewrite_args(args.rewrite_fastq_dir)
    all_libs = apply_fastq_rewrites(parsed["libraries"], rewrites)
    gex_libs, feature_libs, guide_libs = classify_libraries(all_libs)

    pf_multi_out = Path(args.pf_multi_out)
    pf_multi_out.parent.mkdir(parents=True, exist_ok=True)
    write_pf_multi_config(pf_multi_out, sample_id, gex_libs + feature_libs, global_feature_ref)

    if args.emit_env:
        gex_dirs, gex_r1, gex_r2 = collect_pairs_for_libraries(gex_libs)
        feature_dirs, feature_r1, feature_r2 = collect_pairs_for_libraries(feature_libs)
        guide_dirs, guide_r1, guide_r2 = collect_pairs_for_libraries(guide_libs)
        emit_env(
            {
                "CR_CONFIG": str(config_path),
                "SAMPLE_ID": sample_id,
                "CR_FEATURE_REF": global_feature_ref,
                "CR_HAS_GLOBAL_FEATURE_REF": "1" if global_feature_ref else "0",
                "CR_GENE_EXPRESSION_REFERENCE": parsed.get("gene-expression", {}).get("reference", ""),
                "CR_GENE_EXPRESSION_CHEMISTRY": parsed.get("gene-expression", {}).get("chemistry", ""),
                "PF_MULTI_CONFIG": str(pf_multi_out),
                "GEX_FASTQ_DIRS": ",".join(gex_dirs),
                "FEATURE_FASTQ_DIRS": ",".join(feature_dirs),
                "GUIDE_FASTQ_DIRS": ",".join(guide_dirs),
                "GEX_R1": ",".join(gex_r1),
                "GEX_R2": ",".join(gex_r2),
                "FEATURE_R1": ",".join(feature_r1),
                "FEATURE_R2": ",".join(feature_r2),
                "GUIDE_R1": ",".join(guide_r1),
                "GUIDE_R2": ",".join(guide_r2),
            }
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())

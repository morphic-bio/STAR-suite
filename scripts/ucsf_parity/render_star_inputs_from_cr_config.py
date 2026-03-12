#!/usr/bin/env python3
import argparse
import csv
import shlex
import sys
from pathlib import Path


def normalize_feature_type(value: str) -> str:
    return "".join(ch.lower() for ch in value if ch.isalnum())


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

    for section in ("gene-expression", "feature", "reference"):
        if section in sections:
            reader = csv.reader(sections[section])
            for row in reader:
                if len(row) >= 2:
                    parsed.setdefault(section, {})[row[0].strip().lower()] = row[1].strip()

    if "libraries" not in sections:
        raise SystemExit(f"Missing [libraries] section in {path}")

    reader = csv.DictReader(sections["libraries"])
    for row in reader:
        cleaned = {k.strip().lower(): (v.strip() if v is not None else "") for k, v in row.items()}
        if cleaned.get("fastqs"):
            parsed["libraries"].append(cleaned)

    if not parsed["libraries"]:
        raise SystemExit(f"No library rows found in {path}")

    return parsed


def enumerate_pairs(fastq_dir: Path) -> tuple[list[str], list[str]]:
    r1s = sorted(str(p) for p in fastq_dir.glob("*_R1_001.fastq.gz"))
    if not r1s:
        raise SystemExit(f"No R1 FASTQs found under {fastq_dir}")
    r2s = []
    for r1 in r1s:
        r2 = r1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
        if not Path(r2).is_file():
            raise SystemExit(f"Missing R2 pair for {r1}")
        r2s.append(r2)
    return r1s, r2s


def classify_libraries(libraries: list[dict]) -> tuple[list[dict], list[dict]]:
    gex = []
    guides = []
    for row in libraries:
        feature_type = row.get("feature_types") or row.get("feature_type") or row.get("library_type") or ""
        norm = normalize_feature_type(feature_type)
        if "geneexpression" in norm or norm == "gex":
            gex.append(row)
        elif "crisprguidecapture" in norm or "guidecapture" in norm or "crispr" in norm:
            guides.append(row)
    return gex, guides


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


def write_pf_multi_config(path: Path, sample_id: str, gex_libs: list[dict], guide_libs: list[dict], feature_ref: str) -> None:
    with path.open("w", newline="") as handle:
        handle.write("[libraries]\n")
        handle.write("fastqs,sample,library_type,feature_types\n")
        for row in gex_libs:
            handle.write(f"{row['fastqs']},{sample_id},Gene Expression,{row.get('feature_types', 'Gene Expression')}\n")
        for row in guide_libs:
            handle.write(
                f"{row['fastqs']},{sample_id},CRISPR Guide Capture,{row.get('feature_types', 'CRISPR Guide Capture')}\n"
            )
        handle.write("\n[feature]\n")
        handle.write(f"ref,{feature_ref}\n")


def emit_env(values: dict[str, str]) -> None:
    for key, value in values.items():
        print(f"{key}={shlex.quote(value)}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Render STAR UCSF inputs from a Cell Ranger config.csv")
    parser.add_argument("--config", required=True, help="Cell Ranger config.csv")
    parser.add_argument("--pf-multi-out", required=True, help="Output pf_multi_config.csv path")
    parser.add_argument("--sample-id", help="Override sample ID for generated pf_multi_config.csv")
    parser.add_argument("--emit-env", action="store_true", help="Emit shell-compatible KEY=VALUE lines")
    args = parser.parse_args()

    config_path = Path(args.config)
    parsed = parse_cr_config(config_path)
    gex_libs, guide_libs = classify_libraries(parsed["libraries"])

    if not gex_libs:
        raise SystemExit(f"No Gene Expression library rows found in {config_path}")
    if not guide_libs:
        raise SystemExit(f"No CRISPR Guide Capture library rows found in {config_path}")

    sample_id = determine_sample_id(parsed["libraries"], args.sample_id)
    feature_ref = parsed.get("feature", {}).get("reference") or parsed.get("feature", {}).get("ref")
    if not feature_ref:
        raise SystemExit(f"No feature reference found in [feature] section of {config_path}")

    gex_r1, gex_r2, guide_r1, guide_r2 = [], [], [], []
    for row in gex_libs:
        r1s, r2s = enumerate_pairs(Path(row["fastqs"]))
        gex_r1.extend(r1s)
        gex_r2.extend(r2s)
    for row in guide_libs:
        r1s, r2s = enumerate_pairs(Path(row["fastqs"]))
        guide_r1.extend(r1s)
        guide_r2.extend(r2s)

    pf_multi_out = Path(args.pf_multi_out)
    pf_multi_out.parent.mkdir(parents=True, exist_ok=True)
    write_pf_multi_config(pf_multi_out, sample_id, gex_libs, guide_libs, feature_ref)

    if args.emit_env:
        emit_env(
            {
                "CR_CONFIG": str(config_path),
                "SAMPLE_ID": sample_id,
                "CR_FEATURE_REF": feature_ref,
                "CR_GENE_EXPRESSION_REFERENCE": parsed.get("gene-expression", {}).get("reference", ""),
                "CR_GENE_EXPRESSION_CHEMISTRY": parsed.get("gene-expression", {}).get("chemistry", ""),
                "PF_MULTI_CONFIG": str(pf_multi_out),
                "GEX_R1": ",".join(gex_r1),
                "GEX_R2": ",".join(gex_r2),
                "GUIDE_R1": ",".join(guide_r1),
                "GUIDE_R2": ",".join(guide_r2),
            }
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())

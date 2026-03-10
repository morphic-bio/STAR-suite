#!/usr/bin/env python3
import argparse
import importlib.util
import re
import sys
from pathlib import Path


def load_compare_module(script_path: Path):
    spec = importlib.util.spec_from_file_location("compare_salmon_star", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module from {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_log(log_path: Path) -> dict[str, object]:
    text = log_path.read_text()
    patterns = {
        "votes": r"Library format votes:\s*(.+)",
        "detected_format": r"Detected library format:\s+([A-Z0-9_]+)",
        "dropped_incompat": r"dropped_incompat:\s+(\d+)",
        "dropped_missing_mate_fields": r"dropped_missing_mate_fields:\s+(\d+)",
        "dropped_unknown_obs_fmt": r"dropped_unknown_obs_fmt:\s+(\d+)",
    }
    parsed: dict[str, object] = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if not match:
            raise RuntimeError(f"Could not parse {key} from {log_path}")
        parsed[key] = match.group(1)
    for key in ("dropped_incompat", "dropped_missing_mate_fields", "dropped_unknown_obs_fmt"):
        parsed[key] = int(parsed[key])  # type: ignore[index]
    return parsed


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate PE TranscriptVB auto-detect output")
    parser.add_argument("--outdir", required=True, help="STAR output directory")
    parser.add_argument("--expected-format", default="", help="Exact detected format expected")
    parser.add_argument("--forbid-format", action="append", default=["UNKNOWN"], help="Detected formats that must not appear")
    parser.add_argument("--max-dropped-incompat", type=int, default=0)
    parser.add_argument("--max-dropped-missing-mate-fields", type=int, default=0)
    parser.add_argument("--max-dropped-unknown-obs-fmt", type=int, default=0)
    parser.add_argument("--salmon-quant", default="", help="Optional Salmon quant.sf for parity checks")
    parser.add_argument("--star-quant", default="", help="Optional STAR quant.sf override")
    parser.add_argument("--min-spearman-all", type=float, default=0.0)
    parser.add_argument("--min-spearman-expressed", type=float, default=0.0)
    parser.add_argument("--min-pearson-expressed", type=float, default=0.0)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    log_path = outdir / "Log.out"
    if not log_path.exists():
        print(f"ERROR: missing {log_path}", file=sys.stderr)
        return 2

    parsed = parse_log(log_path)
    detected = str(parsed["detected_format"])
    print(f"Votes: {parsed['votes']}")
    print(f"Detected format: {detected}")
    print(f"dropped_incompat: {parsed['dropped_incompat']}")
    print(f"dropped_missing_mate_fields: {parsed['dropped_missing_mate_fields']}")
    print(f"dropped_unknown_obs_fmt: {parsed['dropped_unknown_obs_fmt']}")

    ok = True

    if args.expected_format and detected != args.expected_format:
        print(f"FAIL: expected detected format {args.expected_format}, saw {detected}", file=sys.stderr)
        ok = False
    if detected in set(args.forbid_format):
        print(f"FAIL: forbidden detected format {detected}", file=sys.stderr)
        ok = False
    if args.max_dropped_incompat >= 0 and parsed["dropped_incompat"] > args.max_dropped_incompat:
        print(
            f"FAIL: dropped_incompat {parsed['dropped_incompat']} exceeds {args.max_dropped_incompat}",
            file=sys.stderr,
        )
        ok = False
    if (
        args.max_dropped_missing_mate_fields >= 0
        and parsed["dropped_missing_mate_fields"] > args.max_dropped_missing_mate_fields
    ):
        print(
            "FAIL: dropped_missing_mate_fields "
            f"{parsed['dropped_missing_mate_fields']} exceeds {args.max_dropped_missing_mate_fields}",
            file=sys.stderr,
        )
        ok = False
    if args.max_dropped_unknown_obs_fmt >= 0 and parsed["dropped_unknown_obs_fmt"] > args.max_dropped_unknown_obs_fmt:
        print(
            "FAIL: dropped_unknown_obs_fmt "
            f"{parsed['dropped_unknown_obs_fmt']} exceeds {args.max_dropped_unknown_obs_fmt}",
            file=sys.stderr,
        )
        ok = False

    if args.salmon_quant:
        compare_script = Path(__file__).with_name("compare_salmon_star.py")
        compare_module = load_compare_module(compare_script)
        salmon_quant = Path(args.salmon_quant)
        star_quant = Path(args.star_quant) if args.star_quant else outdir / "quant.sf"
        if not salmon_quant.exists():
            print(f"ERROR: missing Salmon quant file {salmon_quant}", file=sys.stderr)
            return 2
        if not star_quant.exists():
            print(f"ERROR: missing STAR quant file {star_quant}", file=sys.stderr)
            return 2
        results = compare_module.compare_quantifications(
            compare_module.load_quant_sf(str(salmon_quant)),
            compare_module.load_quant_sf(str(star_quant)),
            verbose=False,
        )
        compare_module.print_results(results, verbose=False)

        metric_thresholds = {
            "spearman_all": args.min_spearman_all,
            "spearman_expressed": args.min_spearman_expressed,
            "pearson_expressed": args.min_pearson_expressed,
        }
        for metric, threshold in metric_thresholds.items():
            if threshold <= 0:
                continue
            value = results[metric]
            if value < threshold:
                print(f"FAIL: {metric} {value:.4f} < {threshold:.4f}", file=sys.stderr)
                ok = False

    if ok:
        print("PASS: TranscriptVB PE auto-detect validation succeeded")
        return 0
    print("FAIL: TranscriptVB PE auto-detect validation failed", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())

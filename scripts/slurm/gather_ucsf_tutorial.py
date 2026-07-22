#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Gather finalized UCSF sample manifests.")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sample-manifest", action="append", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    samples = []
    for raw_path in args.sample_manifest:
        path = Path(raw_path)
        payload = json.loads(path.read_text())
        if payload.get("schema") != "star_suite.ucsf_downstream_finalize/v1":
            raise SystemExit(f"unexpected sample manifest schema: {path}")
        if payload.get("status") != "complete":
            raise SystemExit(f"sample is not complete: {path}")
        sample = path.parent.parent.name
        samples.append(
            {
                "sample": sample,
                "manifest": str(path),
                "shape": payload.get("shape"),
                "final_counts_h5ad": payload["artifacts"]["final_counts.h5ad"]["path"],
                "cellbender_h5": payload["cellbender_h5"],
            }
        )

    samples.sort(key=lambda item: item["sample"])
    payload = {
        "schema": "star_suite.ucsf_tutorial_set/v1",
        "status": "complete",
        "run_id": args.run_id,
        "completed_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "sample_count": len(samples),
        "samples": samples,
    }
    manifest = output_dir / "tutorial-set.json"
    manifest.write_text(json.dumps(payload, indent=2) + "\n")

    lines = ["sample\tn_obs\tn_vars\tfinal_counts_h5ad\tcellbender_h5"]
    for sample in samples:
        n_obs, n_vars = sample["shape"]
        lines.append(
            f"{sample['sample']}\t{n_obs}\t{n_vars}\t"
            f"{sample['final_counts_h5ad']}\t{sample['cellbender_h5']}"
        )
    (output_dir / "samples.tsv").write_text("\n".join(lines) + "\n")
    print(f"Wrote {manifest}")


if __name__ == "__main__":
    main()

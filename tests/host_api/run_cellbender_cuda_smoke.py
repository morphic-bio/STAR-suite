#!/usr/bin/env python3
"""CUDA/layer-integration smoke on a raw MEX with real cells and empty droplets.

Five epochs test execution, not biological convergence. Never replace production
deliverables with these smoke outputs. Run under the shared host benchmark lock.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import os


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-mex", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--recipes-root", type=Path, required=True)
    ap.add_argument("--expected-cells", type=int, required=True)
    ap.add_argument("--total-droplets", type=int, default=5000)
    ap.add_argument("--retained-droplets", type=int, default=20000)
    ap.add_argument("--epochs", type=int, default=5)
    ap.add_argument("--image", default="biodepot/cellbender:0.3.2")
    ap.add_argument("--cellbender-gpu", action="store_true", required=True)
    args = ap.parse_args()
    if not 0 < args.expected_cells < args.total_droplets < args.retained_droplets:
        ap.error("require 0 < expected cells < total droplets < retained droplets")
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    import scanpy as sc
    import numpy as np
    import anndata as ad

    raw = sc.read_10x_mtx(args.raw_mex, var_names="gene_ids", gex_only=True)
    totals = np.asarray(raw.X.sum(axis=1)).ravel()
    order = np.argsort(-totals, kind="stable")[:args.retained_droplets]
    raw = raw[order].copy()
    raw.X = raw.X.astype(np.int32)
    raw.write_h5ad(out / "raw_counts.h5ad")
    assert raw.n_obs == args.retained_droplets
    assert totals[order[args.total_droplets]] > 5, "fixture lacks ambient droplets"

    image_id = subprocess.check_output(
        ["docker", "image", "inspect", args.image, "--format", "{{.Id}}"], text=True).strip()
    docker = ["docker", "run", "--rm", "--user", f"{os.getuid()}:{os.getgid()}",
              "-v", f"{out}:{out}", "-w", str(out), "-e", f"HOME={out}"]
    command = [*docker, "--gpus", "all", image_id, "cellbender", "remove-background",
               "--cuda", "--input", str(out / "raw_counts.h5ad"),
               "--output", str(out / "cellbender_counts.h5"),
               "--expected-cells", str(args.expected_cells),
               "--total-droplets-included", str(args.total_droplets),
               "--epochs", str(args.epochs), "--cpu-threads", "4",
               "--posterior-batch-size", "32", "--num-training-tries", "1"]
    helper = args.recipes_root.resolve() / "scripts/add_cellbender_layer_from_h5.py"
    (out / "manifest.json").write_text(json.dumps({
        "purpose": "CUDA functional smoke, not converged production denoising",
        "raw_mex": str(args.raw_mex.resolve()), "input_shape": list(raw.shape),
        "rank_umis": {str(i): int(totals[order[i-1]]) for i in
                      (1, args.expected_cells, args.total_droplets, args.retained_droplets)},
        "image_id": image_id, "command": command,
        "layer_helper": str(helper),
        "layer_helper_sha256": hashlib.sha256(helper.read_bytes()).hexdigest(),
    }, indent=2) + "\n")
    (out / "command.sh").write_text(shlex.join(command) + "\n")
    with (out / "cellbender.log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    assert (out / "cellbender_counts.h5").stat().st_size > 0
    with (out / "layer.log").open("w") as log:
        subprocess.run([*docker, "-v", f"{helper}:{helper}:ro", image_id, "python", str(helper),
                        "--cellbender-h5", str(out / "cellbender_counts.h5"),
                        "--input-h5ad", str(out / "raw_counts.h5ad"),
                        "--output-h5ad", str(out / "denoised_counts.h5ad")],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    result = ad.read_h5ad(out / "denoised_counts.h5ad")
    assert result.obs_names.equals(raw.obs_names) and result.var_names.equals(raw.var_names)
    assert (result.X != raw.X).nnz == 0, "layer integration changed raw counts"
    layer = result.layers["denoised"]
    assert layer.shape == raw.shape and layer.nnz > 0
    assert np.isfinite(layer.data).all() and (layer.data >= 0).all()
    assert not ((layer - raw.X).data > 0).any(), "denoised counts exceed raw counts"
    (out / "PASS.json").write_text(json.dumps({
        "shape": list(result.shape), "raw_umis": int(raw.X.sum()),
        "denoised_umis": int(layer.sum()), "denoised_nnz": layer.nnz,
        "epochs": args.epochs, "production_convergence_validated": False,
    }, indent=2) + "\n")
    print(f"PASS: CUDA inference and H5AD layer integration; evidence: {out}")


if __name__ == "__main__":
    main()

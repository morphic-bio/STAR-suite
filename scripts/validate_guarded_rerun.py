#!/usr/bin/env python3
"""Validate guarded EmptyDrops rerun outputs against old production outputs.

Checks per sample:
1. Barcode consistency across all surfaces (GeneFull, Velocyto, outs MEX, h5ads)
2. Jaccard overlap old vs new filtered barcodes
3. CellBender 'denoised' layer present in h5ads
4. Feature library (CRISPR) annotations carried over
5. Summary stat comparison (median genes, UMIs, mito%)
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np


SAMPLES = [
    "EBs1_1", "EBs1_2", "EBs1_3", "EBs1_4", "EBs1_5",
    "EBs2_1", "EBs2_2", "EBs2_3", "EBs2_4", "EBs2_5",
    "iPSC1_1", "iPSC1_2", "iPSC1_3",
    "iPSC2_1", "iPSC2_2", "iPSC2_3",
]


def read_barcodes(path: Path) -> list[str]:
    if not path.exists():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return [line.strip() for line in f if line.strip()]


def strip_suffix(bc: str) -> str:
    return bc[:-2] if bc.endswith("-1") else bc


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    union = a | b
    return len(a & b) / len(union) if union else 0.0


def check_barcode_consistency(sample_dir: Path, new_ds: Path) -> dict:
    """Check barcode counts agree across all surfaces."""
    results = {}

    # GeneFull filtered
    gf = read_barcodes(sample_dir / "run/Solo.out/GeneFull/filtered/barcodes.tsv")
    results["genefull_filtered"] = len(gf)

    # Velocyto filtered
    vf = read_barcodes(sample_dir / "run/Solo.out/Velocyto/filtered/barcodes.tsv")
    results["velocyto_filtered"] = len(vf)

    # outs filtered MEX
    ofm = read_barcodes(sample_dir / "run/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
    results["outs_filtered_mex"] = len(ofm)

    # outs filtered velocyto MEX
    ovfm = read_barcodes(sample_dir / "run/outs/filtered_velocyto_feature_bc_matrix/barcodes.tsv.gz")
    results["outs_filtered_velo_mex"] = len(ovfm)

    # h5ads
    for name in ("counts.h5ad", "filtered_counts.h5ad", "default_singlet_filtered_counts.h5ad"):
        h5path = new_ds / name
        if h5path.exists():
            try:
                adata = ad.read_h5ad(h5path, backed="r")
                results[f"h5ad_{name}"] = adata.n_obs
                adata.file.close()
            except Exception as e:
                results[f"h5ad_{name}"] = f"ERROR: {e}"
        else:
            results[f"h5ad_{name}"] = "MISSING"

    # Normalize to stripped barcodes for comparison
    gf_set = {strip_suffix(bc) for bc in gf}
    ofm_set = {strip_suffix(bc) for bc in ofm}
    vf_set = {strip_suffix(bc) for bc in vf}
    ovfm_set = {strip_suffix(bc) for bc in ovfm}

    results["genefull_vs_outs_match"] = gf_set == ofm_set
    results["genefull_vs_velocyto_match"] = gf_set == vf_set
    results["genefull_vs_velo_outs_match"] = gf_set == ovfm_set

    return results


def check_jaccard(sample_dir: Path) -> dict:
    """Jaccard overlap between old and new filtered barcodes."""
    old_path = sample_dir / "run/Solo.out/GeneFull/filtered_pre_guarded_backup/barcodes.tsv"
    new_path = sample_dir / "run/Solo.out/GeneFull/filtered/barcodes.tsv"

    old_bcs = {strip_suffix(bc) for bc in read_barcodes(old_path)}
    new_bcs = {strip_suffix(bc) for bc in read_barcodes(new_path)}

    return {
        "old_count": len(old_bcs),
        "new_count": len(new_bcs),
        "intersection": len(old_bcs & new_bcs),
        "old_only": len(old_bcs - new_bcs),
        "new_only": len(new_bcs - old_bcs),
        "jaccard": jaccard(old_bcs, new_bcs),
    }


def check_cellbender_layer(new_ds: Path) -> dict:
    """Check that 'denoised' layer is present in h5ads."""
    results = {}
    for name in ("counts.h5ad", "filtered_counts.h5ad", "default_singlet_filtered_counts.h5ad"):
        h5path = new_ds / name
        if not h5path.exists():
            results[name] = "MISSING"
            continue
        try:
            adata = ad.read_h5ad(h5path, backed="r")
            has_layer = "denoised" in adata.layers
            results[name] = "OK" if has_layer else "NO denoised layer"
            adata.file.close()
        except Exception as e:
            results[name] = f"ERROR: {e}"
    return results


def check_feature_libraries(new_ds: Path) -> dict:
    """Check CRISPR feature library integration."""
    results = {}
    feat_dir = new_ds / "feature_libraries"
    results["feature_libraries_dir"] = feat_dir.exists()

    for name in ("counts.h5ad", "filtered_counts.h5ad"):
        h5path = new_ds / name
        if not h5path.exists():
            results[f"{name}_obsm"] = "MISSING"
            continue
        try:
            adata = ad.read_h5ad(h5path, backed="r")
            crispr_keys = [k for k in adata.obsm.keys() if "crispr" in k.lower() or "guide" in k.lower()]
            obs_crispr = [c for c in adata.obs.columns if "crispr" in c.lower() or "guide" in c.lower()]
            results[f"{name}_crispr_obsm"] = crispr_keys if crispr_keys else "none"
            results[f"{name}_crispr_obs"] = obs_crispr if obs_crispr else "none"
            adata.file.close()
        except Exception as e:
            results[f"{name}_obsm"] = f"ERROR: {e}"
    return results


def compare_summaries(old_ds: Path, new_ds: Path) -> dict:
    """Compare key stats from old vs new summary.txt."""
    results = {}
    for label, ds in [("old", old_ds), ("new", new_ds)]:
        summary_path = ds / "summary.txt"
        if not summary_path.exists():
            results[label] = "MISSING"
            continue
        text = summary_path.read_text()
        stats = {}
        for line in text.splitlines():
            for key in ("n_obs", "n_vars", "median_genes", "median_counts", "median_pct_mito"):
                if key in line.lower() or key.replace("_", " ") in line.lower():
                    parts = line.split(":")
                    if len(parts) >= 2:
                        try:
                            stats[key] = parts[-1].strip()
                        except Exception:
                            pass
        results[label] = stats if stats else text[:200]
    return results


def validate_sample(prod_root: Path, sample: str) -> dict:
    sample_dir = prod_root / "samples" / sample
    old_ds = sample_dir / "downstream_genefull_velocyto_cellbender"
    new_ds = sample_dir / "downstream_genefull_velocyto_cellbender_guarded_rerun"

    report = {"sample": sample}

    # 1. Barcode consistency
    report["barcode_consistency"] = check_barcode_consistency(sample_dir, new_ds)

    # 2. Jaccard
    report["jaccard"] = check_jaccard(sample_dir)

    # 3. CellBender layer
    report["cellbender_layer"] = check_cellbender_layer(new_ds)

    # 4. Feature libraries
    report["feature_libraries"] = check_feature_libraries(new_ds)

    # 5. Summary comparison
    report["summary_comparison"] = compare_summaries(old_ds, new_ds)

    return report


def print_summary_table(reports: list[dict]) -> None:
    print("\n" + "=" * 120)
    print("VALIDATION SUMMARY")
    print("=" * 120)

    # Barcode consistency table
    print("\n--- Barcode Consistency ---")
    print(f"{'Sample':<10} {'GF filt':>8} {'Velo filt':>10} {'Outs MEX':>9} {'Velo MEX':>9} "
          f"{'GF=Outs':>8} {'GF=Velo':>8} {'GF=VeloO':>9}")
    for r in reports:
        bc = r["barcode_consistency"]
        print(f"{r['sample']:<10} {bc['genefull_filtered']:>8} {bc['velocyto_filtered']:>10} "
              f"{bc['outs_filtered_mex']:>9} {bc['outs_filtered_velo_mex']:>9} "
              f"{'OK' if bc['genefull_vs_outs_match'] else 'FAIL':>8} "
              f"{'OK' if bc['genefull_vs_velocyto_match'] else 'FAIL':>8} "
              f"{'OK' if bc['genefull_vs_velo_outs_match'] else 'FAIL':>9}")

    # Jaccard table
    print("\n--- Old vs New Jaccard ---")
    print(f"{'Sample':<10} {'Old':>8} {'New':>8} {'Inter':>8} {'Old-only':>9} {'New-only':>9} {'Jaccard':>9}")
    for r in reports:
        j = r["jaccard"]
        print(f"{r['sample']:<10} {j['old_count']:>8} {j['new_count']:>8} {j['intersection']:>8} "
              f"{j['old_only']:>9} {j['new_only']:>9} {j['jaccard']:>9.6f}")

    # CellBender layer table
    print("\n--- CellBender 'denoised' Layer ---")
    print(f"{'Sample':<10} {'counts.h5ad':>15} {'filtered_counts':>18} {'default_singlet':>18}")
    for r in reports:
        cb = r["cellbender_layer"]
        print(f"{r['sample']:<10} {cb.get('counts.h5ad', '?'):>15} "
              f"{cb.get('filtered_counts.h5ad', '?'):>18} "
              f"{cb.get('default_singlet_filtered_counts.h5ad', '?'):>18}")

    # Feature libraries
    print("\n--- Feature Libraries ---")
    print(f"{'Sample':<10} {'Dir exists':>11} {'counts crispr_obs':>20} {'filtered crispr_obs':>22}")
    for r in reports:
        fl = r["feature_libraries"]
        cobs_c = str(fl.get("counts.h5ad_crispr_obs", "?"))
        cobs_f = str(fl.get("filtered_counts.h5ad_crispr_obs", "?"))
        print(f"{r['sample']:<10} {'OK' if fl.get('feature_libraries_dir') else 'MISSING':>11} "
              f"{cobs_c[:20]:>20} {cobs_f[:22]:>22}")

    # h5ad obs counts
    print("\n--- h5ad Cell Counts ---")
    print(f"{'Sample':<10} {'counts':>10} {'filtered':>10} {'def_singlet':>12}")
    for r in reports:
        bc = r["barcode_consistency"]
        c = bc.get("h5ad_counts.h5ad", "?")
        f = bc.get("h5ad_filtered_counts.h5ad", "?")
        d = bc.get("h5ad_default_singlet_filtered_counts.h5ad", "?")
        print(f"{r['sample']:<10} {str(c):>10} {str(f):>10} {str(d):>12}")

    # Overall pass/fail
    print("\n--- Overall ---")
    all_pass = True
    for r in reports:
        bc = r["barcode_consistency"]
        cb = r["cellbender_layer"]
        issues = []
        if not bc["genefull_vs_outs_match"]:
            issues.append("GF!=Outs")
        if not bc["genefull_vs_velocyto_match"]:
            issues.append("GF!=Velo")
        if not bc["genefull_vs_velo_outs_match"]:
            issues.append("GF!=VeloOuts")
        for name in ("counts.h5ad", "filtered_counts.h5ad", "default_singlet_filtered_counts.h5ad"):
            if cb.get(name) != "OK":
                issues.append(f"CB:{name}")
        status = "PASS" if not issues else f"FAIL: {', '.join(issues)}"
        if issues:
            all_pass = False
        print(f"  {r['sample']:<10} {status}")

    print(f"\n{'ALL SAMPLES PASS' if all_pass else 'SOME SAMPLES HAVE ISSUES'}")
    print("=" * 120)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prod-root", required=True, help="Production run root directory")
    parser.add_argument("--samples", nargs="*", default=SAMPLES, help="Sample names (default: all 16)")
    parser.add_argument("--json-output", help="Optional JSON output path")
    args = parser.parse_args()

    prod_root = Path(args.prod_root).resolve()

    reports = []
    for sample in args.samples:
        print(f"Validating {sample}...", file=sys.stderr)
        report = validate_sample(prod_root, sample)
        reports.append(report)

    print_summary_table(reports)

    if args.json_output:
        out_path = Path(args.json_output)
        # Convert sets and non-serializable types
        def make_serializable(obj):
            if isinstance(obj, set):
                return sorted(obj)
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            return obj

        import json as json_mod
        with open(out_path, "w") as f:
            json_mod.dump(reports, f, indent=2, default=make_serializable)
        print(f"\nJSON report written to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()

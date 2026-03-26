#!/usr/bin/env python3
"""Exact (or soft) comparison of Velocyto outputs between two STAR run directories.

Modes:
  solo      — Solo.out/Velocyto/raw matrices + mandatory Gene/raw axis tables (barcodes + features)
  packaged  — gz MEX: per-run manifest vs counts, matrix axes vs tables, raw/filtered feature identity
              and filtered⊆raw barcodes, layer-sum consistency, then cross-run diff
  all       — both (required for full Stage 1 acceptance)
  genes     — Solo.out/Gene and Solo.out/GeneFull raw+filtered only (Phase 6: Velocyto path must not alter GEX)

Velocyto/raw usually has only *.mtx; axes are taken from Solo.out/Gene/raw (same as STAR output).
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

from scipy.io import mmread
from scipy.sparse import issparse


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_single_column(path: Path) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with open_text(path) as fh:
        return [line.split("\t")[0].rstrip("\n") for line in fh if line.strip()]


def read_tsv_rows(path: Path) -> list[list[str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with open_text(path) as fh:
        return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


def read_feature_rows_gz(path: Path) -> list[list[str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with gzip.open(path, "rt") as fh:
        return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


def read_gz_lines(path: Path) -> list[str]:
    """Non-empty logical lines from a gz text file (matches perturb velocyto MEX smoke)."""
    if not path.is_file():
        raise FileNotFoundError(path)
    with gzip.open(path, "rt") as fh:
        return [line.rstrip("\n") for line in fh if line.strip()]


def load_sparse_mtx(path: Path):
    if not path.is_file():
        raise FileNotFoundError(path)
    opener = gzip.open if path.suffix == ".gz" else open
    mode = "rb" if path.suffix == ".gz" else "r"
    with opener(path, mode) as fh:
        m = mmread(fh)
    if not issparse(m):
        raise TypeError(f"Expected sparse Matrix Market input: {path}")
    return m.tocsr()


def matrices_equal(a, b, atol: float = 0.0) -> bool:
    if a.shape != b.shape:
        return False
    if atol == 0.0:
        return (a != b).nnz == 0
    return (a - b).nnz == 0 or abs((a - b).data).max() <= atol


def compare_pair_mtx(
    label: str,
    a: Path,
    b: Path,
    exact: bool,
) -> list[str]:
    errors: list[str] = []
    ma = load_sparse_mtx(a)
    mb = load_sparse_mtx(b)
    if ma.shape != mb.shape:
        errors.append(f"{label}: shape mismatch {ma.shape} vs {mb.shape}")
        return errors
    if exact and not matrices_equal(ma, mb, 0.0):
        errors.append(f"{label}: matrix data differs (nnz {ma.nnz} vs {mb.nnz})")
    elif not exact and not matrices_equal(ma, mb, 1e-9):
        errors.append(f"{label}: matrix data differs beyond tolerance (nnz {ma.nnz} vs {mb.nnz})")
    return errors


def compare_solo_velocyto_raw(
    run_a: Path,
    run_b: Path,
    velocyto_subdir: Path,
    gene_raw_subdir: Path,
    exact: bool,
) -> list[str]:
    base_a = run_a / velocyto_subdir
    base_b = run_b / velocyto_subdir
    all_err: list[str] = []

    # Velocyto matrices use the same feature / raw-barcode axes as Gene (STAR writes mtx in that space).
    ga = run_a / gene_raw_subdir
    gb = run_b / gene_raw_subdir
    for rel, reader, label in (
        ("barcodes.tsv", read_single_column, "solo/Gene/raw barcodes.tsv"),
        ("features.tsv", read_tsv_rows, "solo/Gene/raw features.tsv"),
    ):
        pa = ga / rel
        pb = gb / rel
        if not pa.is_file():
            all_err.append(f"{label}: missing in run A ({pa})")
            continue
        if not pb.is_file():
            all_err.append(f"{label}: missing in run B ({pb})")
            continue
        ca = reader(pa)
        cb = reader(pb)
        if ca != cb:
            all_err.append(f"{label}: axis mismatch ({len(ca)} vs {len(cb)} rows or content differs)")

    if all_err:
        return all_err

    n_feat_a = len(read_tsv_rows(ga / "features.tsv"))
    n_bc_a = len(read_single_column(ga / "barcodes.tsv"))
    n_feat_b = len(read_tsv_rows(gb / "features.tsv"))
    n_bc_b = len(read_single_column(gb / "barcodes.tsv"))

    layers = ("spliced.mtx", "unspliced.mtx", "ambiguous.mtx")
    for layer in layers:
        all_err.extend(compare_pair_mtx(f"solo/{layer}", base_a / layer, base_b / layer, exact))

    sa = load_sparse_mtx(base_a / "spliced.mtx")
    ua = load_sparse_mtx(base_a / "unspliced.mtx")
    aa = load_sparse_mtx(base_a / "ambiguous.mtx")
    sb = load_sparse_mtx(base_b / "spliced.mtx")
    ub = load_sparse_mtx(base_b / "unspliced.mtx")
    ab = load_sparse_mtx(base_b / "ambiguous.mtx")

    for name, m in (("spliced", sa), ("unspliced", ua), ("ambiguous", aa)):
        if m.shape != (n_feat_a, n_bc_a):
            all_err.append(
                f"solo/{name}.mtx run A: shape {m.shape} != (features={n_feat_a}, barcodes={n_bc_a})"
            )
    for name, m in (("spliced", sb), ("unspliced", ub), ("ambiguous", ab)):
        if m.shape != (n_feat_b, n_bc_b):
            all_err.append(
                f"solo/{name}.mtx run B: shape {m.shape} != (features={n_feat_b}, barcodes={n_bc_b})"
            )

    total_a = sa + ua + aa
    total_b = sb + ub + ab
    if exact and not matrices_equal(total_a, total_b, 0.0):
        all_err.append("solo: reconstructed total matrix mismatch")
    elif not exact and not matrices_equal(total_a, total_b, 1e-9):
        all_err.append("solo: reconstructed total matrix mismatch (soft)")

    return all_err


def packaged_matrix_sum_consistency(directory: Path, tag: str) -> list[str]:
    """matrix.mtx.gz must equal spliced + unspliced + ambiguous (per packaged tree)."""
    names = ("spliced.mtx.gz", "unspliced.mtx.gz", "ambiguous.mtx.gz", "matrix.mtx.gz")
    paths = [directory / n for n in names]
    if not all(p.is_file() for p in paths):
        return []
    errs: list[str] = []
    try:
        sp = load_sparse_mtx(paths[0])
        up = load_sparse_mtx(paths[1])
        ap = load_sparse_mtx(paths[2])
        tot = load_sparse_mtx(paths[3])
        if not matrices_equal(tot, sp + up + ap, 0.0):
            errs.append(f"{tag}: matrix.mtx.gz != spliced+unspliced+ambiguous")
    except (OSError, FileNotFoundError, TypeError) as exc:
        errs.append(f"{tag}: packaged layer/total consistency check failed: {exc}")
    return errs


def compare_one_packaged_tree(a_dir: Path, b_dir: Path, label: str, exact: bool) -> list[str]:
    errs: list[str] = []
    for fname in (
        "barcodes.tsv.gz",
        "features.tsv.gz",
        "matrix.mtx.gz",
        "spliced.mtx.gz",
        "unspliced.mtx.gz",
        "ambiguous.mtx.gz",
    ):
        pa = a_dir / fname
        pb = b_dir / fname
        if not pa.is_file():
            errs.append(f"{label}/{fname}: missing in run A")
            continue
        if not pb.is_file():
            errs.append(f"{label}/{fname}: missing in run B")
            continue
        if fname.endswith(".tsv.gz"):
            if fname.startswith("barcodes"):
                ca = read_single_column(pa)
                cb = read_single_column(pb)
                if ca != cb:
                    errs.append(f"{label}/{fname}: barcode list mismatch")
            else:
                ra = read_feature_rows_gz(pa)
                rb = read_feature_rows_gz(pb)
                if ra != rb:
                    errs.append(f"{label}/{fname}: feature table mismatch")
        else:
            errs.extend(compare_pair_mtx(f"{label}/{fname}", pa, pb, exact))

    errs.extend(packaged_matrix_sum_consistency(a_dir, f"{label}/run_A"))
    errs.extend(packaged_matrix_sum_consistency(b_dir, f"{label}/run_B"))

    return errs


def check_packaged_tree_matrix_axes(pkg_dir: Path, label: str) -> list[str]:
    """Every layer must be (n_features, n_barcodes) for that tree's gz axis tables."""
    bc_path = pkg_dir / "barcodes.tsv.gz"
    ft_path = pkg_dir / "features.tsv.gz"
    errs: list[str] = []
    try:
        bc_lines = read_gz_lines(bc_path)
        ft_lines = read_gz_lines(ft_path)
    except FileNotFoundError as exc:
        return [f"{label}: missing axis file ({exc})"]
    n_bc = len(bc_lines)
    n_ft = len(ft_lines)
    expected = (n_ft, n_bc)
    for fname in ("spliced.mtx.gz", "unspliced.mtx.gz", "ambiguous.mtx.gz", "matrix.mtx.gz"):
        p = pkg_dir / fname
        if not p.is_file():
            errs.append(f"{label}: missing {fname}")
            continue
        try:
            m = load_sparse_mtx(p)
        except (OSError, TypeError) as exc:
            errs.append(f"{label}/{fname}: load failed: {exc}")
            continue
        if m.shape != expected:
            errs.append(
                f"{label}/{fname}: shape {m.shape} != (features={n_ft}, barcodes={n_bc}) from axis tables"
            )
    return errs


def validate_packaged_velocyto_cross_split(run_root: Path, tag: str) -> list[str]:
    """Per-run invariants: axis shapes, raw/filtered feature identity, filtered barcodes subset of raw."""
    raw_dir = run_root / "outs" / "raw_velocyto_feature_bc_matrix"
    filt_dir = run_root / "outs" / "filtered_velocyto_feature_bc_matrix"
    errs: list[str] = []

    for d, lbl in ((raw_dir, f"{tag} raw"), (filt_dir, f"{tag} filtered")):
        if not d.is_dir():
            errs.append(f"{lbl}: missing directory {d}")
            continue
        errs.extend(check_packaged_tree_matrix_axes(d, lbl))

    if not raw_dir.is_dir() or not filt_dir.is_dir():
        return errs

    try:
        raw_feat = read_gz_lines(raw_dir / "features.tsv.gz")
        raw_bc = read_gz_lines(raw_dir / "barcodes.tsv.gz")
        filt_feat = read_gz_lines(filt_dir / "features.tsv.gz")
        filt_bc = read_gz_lines(filt_dir / "barcodes.tsv.gz")
    except FileNotFoundError as exc:
        errs.append(f"{tag}: cannot read packaged axis gz ({exc})")
        return errs

    if raw_feat != filt_feat:
        errs.append(f"{tag}: raw and filtered velocyto feature tables differ (line-wise)")
    if not set(filt_bc).issubset(set(raw_bc)):
        errs.append(f"{tag}: filtered barcodes are not a subset of raw barcodes")

    return errs


def measure_packaged_tree(pkg_dir: Path) -> dict[str, int] | None:
    """Return manifest-shaped counts from on-disk packaged MEX, or None if incomplete."""
    names = (
        "barcodes.tsv.gz",
        "features.tsv.gz",
        "spliced.mtx.gz",
        "unspliced.mtx.gz",
        "ambiguous.mtx.gz",
        "matrix.mtx.gz",
    )
    paths = [pkg_dir / n for n in names]
    if not all(p.is_file() for p in paths):
        return None
    try:
        bc = read_single_column(paths[0])
        ft = read_feature_rows_gz(paths[1])
        sp = load_sparse_mtx(paths[2])
        up = load_sparse_mtx(paths[3])
        ap = load_sparse_mtx(paths[4])
        tot = load_sparse_mtx(paths[5])
    except (OSError, FileNotFoundError, TypeError):
        return None
    return {
        "features": len(ft),
        "barcodes": len(bc),
        "nnz_total": int(tot.nnz),
        "nnz_spliced": int(sp.nnz),
        "nnz_unspliced": int(up.nnz),
        "nnz_ambiguous": int(ap.nnz),
    }


def validate_velocyto_manifest_vs_packaged(run_root: Path, tag: str) -> list[str]:
    """Each run's manifest must match that run's packaged dirs (not only cross-run equality)."""
    mpath = run_root / "outs" / "velocyto_feature_bc_matrix_manifest.json"
    if not mpath.is_file():
        return [f"{tag}: missing {mpath}"]
    try:
        data = json.loads(mpath.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        return [f"{tag}: invalid JSON in manifest: {exc}"]

    errs: list[str] = []
    for split, dirname in (
        ("raw", "raw_velocyto_feature_bc_matrix"),
        ("filtered", "filtered_velocyto_feature_bc_matrix"),
    ):
        section = data.get(split)
        if not isinstance(section, dict):
            errs.append(f"{tag}: manifest missing or invalid {split!r} section")
            continue
        pkg = run_root / "outs" / dirname
        measured = measure_packaged_tree(pkg)
        if measured is None:
            errs.append(f"{tag}: incomplete packaged tree {pkg}")
            continue
        for k in ("features", "barcodes", "nnz_total", "nnz_spliced", "nnz_unspliced", "nnz_ambiguous"):
            mv = section.get(k)
            gv = measured.get(k)
            if mv != gv:
                errs.append(f"{tag}: manifest {split}.{k}={mv!r} != on-disk {gv!r} ({pkg})")
    return errs


def compare_packaged_velocyto(
    run_a: Path,
    run_b: Path,
    exact: bool,
) -> list[str]:
    errs: list[str] = []
    errs.extend(validate_velocyto_manifest_vs_packaged(run_a, "packaged/manifest run_A"))
    errs.extend(validate_velocyto_manifest_vs_packaged(run_b, "packaged/manifest run_B"))
    errs.extend(validate_packaged_velocyto_cross_split(run_a, "packaged/invariants run_A"))
    errs.extend(validate_packaged_velocyto_cross_split(run_b, "packaged/invariants run_B"))

    for split in ("raw_velocyto_feature_bc_matrix", "filtered_velocyto_feature_bc_matrix"):
        da = run_a / "outs" / split
        db = run_b / "outs" / split
        errs.extend(compare_one_packaged_tree(da, db, f"packaged/{split}", exact))

    ma = run_a / "outs" / "velocyto_feature_bc_matrix_manifest.json"
    mb = run_b / "outs" / "velocyto_feature_bc_matrix_manifest.json"
    if ma.is_file() and mb.is_file():
        ja = json.loads(ma.read_text(encoding="utf-8"))
        jb = json.loads(mb.read_text(encoding="utf-8"))
        for key in ("raw", "filtered"):
            if ja.get(key) != jb.get(key):
                errs.append(f"packaged/manifest.json: {key} summary differs between runs")
    elif ma.is_file() != mb.is_file():
        errs.append("packaged/manifest.json: present in only one run")

    return errs


def compare_solo_gene_genefull(run_a: Path, run_b: Path, exact: bool) -> list[str]:
    """Exact diff of Solo Gene / GeneFull MEX trees (raw + filtered). Runbook Phase 6 gate."""
    errs: list[str] = []
    for feat in ("Gene", "GeneFull"):
        for split in ("raw", "filtered"):
            da = run_a / "Solo.out" / feat / split
            db = run_b / "Solo.out" / feat / split
            if not da.is_dir() and not db.is_dir():
                continue
            if not da.is_dir():
                errs.append(f"genes: missing directory run A {da}")
                continue
            if not db.is_dir():
                errs.append(f"genes: missing directory run B {db}")
                continue
            base = f"genes/{feat}/{split}"
            for rel, reader, label in (
                ("barcodes.tsv", read_single_column, f"{base}/barcodes.tsv"),
                ("features.tsv", read_tsv_rows, f"{base}/features.tsv"),
            ):
                pa = da / rel
                pb = db / rel
                if not pa.is_file():
                    errs.append(f"{label}: missing in run A ({pa})")
                    continue
                if not pb.is_file():
                    errs.append(f"{label}: missing in run B ({pb})")
                    continue
                try:
                    ca = reader(pa)
                    cb = reader(pb)
                except (OSError, FileNotFoundError) as exc:
                    errs.append(f"{label}: read error {exc}")
                    continue
                if ca != cb:
                    errs.append(f"{label}: content mismatch ({len(ca)} vs {len(cb)} rows or bytes differ)")
            ma = da / "matrix.mtx"
            mb = db / "matrix.mtx"
            if not ma.is_file():
                errs.append(f"{base}/matrix.mtx: missing in run A")
            elif not mb.is_file():
                errs.append(f"{base}/matrix.mtx: missing in run B")
            else:
                errs.extend(compare_pair_mtx(f"{base}/matrix.mtx", ma, mb, exact))
    return errs


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run_a", type=Path, help="STAR --outFileNamePrefix directory")
    p.add_argument("run_b", type=Path, help="Second run directory")
    p.add_argument(
        "--mode",
        choices=("solo", "packaged", "all", "genes"),
        default="solo",
        help="solo|packaged|all|genes (genes=Gene+GeneFull Solo.out only)",
    )
    p.add_argument(
        "--soft",
        action="store_true",
        help="Allow tiny floating tolerance",
    )
    p.add_argument(
        "--velocyto-subdir",
        default="Solo.out/Velocyto/raw",
        help="Relative path for solo mode layer matrices",
    )
    p.add_argument(
        "--gene-raw-subdir",
        default="Solo.out/Gene/raw",
        help="Relative path for mandatory barcodes.tsv / features.tsv (solo mode)",
    )
    args = p.parse_args()
    exact = not args.soft

    all_err: list[str] = []

    if args.mode in ("solo", "all"):
        all_err.extend(
            compare_solo_velocyto_raw(
                args.run_a,
                args.run_b,
                Path(args.velocyto_subdir),
                Path(args.gene_raw_subdir),
                exact,
            )
        )

    if args.mode in ("packaged", "all"):
        all_err.extend(compare_packaged_velocyto(args.run_a, args.run_b, exact))

    if args.mode == "genes":
        all_err.extend(compare_solo_gene_genefull(args.run_a, args.run_b, exact))

    if all_err:
        print("compare_velocyto_mex: FAIL", file=sys.stderr)
        for e in all_err:
            print(f"  {e}", file=sys.stderr)
        return 1

    print("compare_velocyto_mex: PASS", f"mode={args.mode}", f"run_a={args.run_a}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

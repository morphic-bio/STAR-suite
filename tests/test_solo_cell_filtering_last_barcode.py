#!/usr/bin/env python3
"""--runMode soloCellFiltering must consider every barcode of the raw matrix,
and --soloCellFilter TopCells N must report N barcodes.

The loader used the cell counter as a 0-based index while filling its arrays
and returned it without restoring the count, so the last barcode that has
counts was never offered to the cell caller. This test builds a small raw
matrix in which that barcode has the largest total and checks that it is
called.

TopCells took the count at 0-based position N of the sorted totals as its
threshold and so reported N + 1 barcodes. It now takes the count of the N-th
barcode: N barcodes, plus any that tie with the N-th.

    python3 tests/test_solo_cell_filtering_last_barcode.py --star core/legacy/source/STAR --out /tmp/out

No genome and no external data are needed.
"""
import argparse
import subprocess
import sys
from pathlib import Path

GENES = ["G1", "G2", "G3"]
# barcode -> counts for G1, G2, G3. The last barcode with counts has the
# largest total; "TTTT" is listed in barcodes.tsv and has no counts.
CELLS = [
    ("AAAA", (200, 150, 50)),   # 400
    ("AACC", (150, 100, 50)),   # 300
    ("ACAC", (3, 2, 0)),        # 5
    ("CAAC", (0, 2, 1)),        # 3
    ("CCAA", (1, 0, 1)),        # 2
    ("GGGG", (250, 150, 100)),  # 500, last barcode with counts
    ("TTTT", (0, 0, 0)),
]


def write_raw(raw, cell_order):
    raw.mkdir(parents=True)
    (raw / "features.tsv").write_text("".join("%s\t%s\tGene Expression\n" % (g, g) for g in GENES))
    (raw / "barcodes.tsv").write_text("".join(b + "\n" for b, _ in CELLS))
    entries = []
    for ci in cell_order:
        for gi, n in enumerate(CELLS[ci][1]):
            if n:
                entries.append((gi + 1, ci + 1, n))
    with (raw / "matrix.mtx").open("w") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        f.write("%d %d %d\n" % (len(GENES), len(CELLS), len(entries)))
        for e in entries:
            f.write("%d %d %d\n" % e)


def call(star, raw, out, cell_filter):
    out.mkdir(parents=True)
    prefix = str(out) + "/"
    argv = [str(star), "--runMode", "soloCellFiltering", str(raw), prefix,
            "--outFileNamePrefix", prefix, "--soloCellFilter"] + cell_filter
    with (out / "stdout.log").open("w") as log:
        proc = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT)
    assert proc.returncode == 0, "STAR exited %d: %s" % (proc.returncode, " ".join(argv))
    barcodes = (out / "barcodes.tsv").read_text().split()
    totals = {}
    lines = [l for l in (out / "matrix.mtx").read_text().splitlines() if not l.startswith("%")]
    for line in lines[1:]:
        _, c, n = line.split()
        b = barcodes[int(c) - 1]
        totals[b] = totals.get(b, 0) + int(float(n))
    return barcodes, totals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--star", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    star, out = args.star.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    total = {b: sum(c) for b, c in CELLS}
    failures = []

    def check(label, got, want):
        ok = got == want
        print("%s %s: got %s, expected %s" % ("ok  " if ok else "FAIL", label, got, want))
        if not ok:
            failures.append(label)

    # the matrix entries in barcode order, and with the last barcode's entries first:
    # the result must not depend on the order of the entries
    for name, order in (("sorted", list(range(len(CELLS)))), ("last_first", [5, 0, 1, 2, 3, 4, 6])):
        raw = out / name / "raw"
        write_raw(raw, order)
        b, t = call(star, raw, out / name / "top1", ["TopCells", "1"])
        check(name + " TopCells 1 barcodes", sorted(b), ["GGGG"])
        b, t = call(star, raw, out / name / "top2", ["TopCells", "2"])
        check(name + " TopCells 2 barcodes", sorted(b), ["AAAA", "GGGG"])
        check(name + " TopCells 2 totals", t, {"AAAA": total["AAAA"], "GGGG": total["GGGG"]})
        b, t = call(star, raw, out / name / "top6", ["TopCells", "6"])
        check(name + " TopCells 6 barcodes", sorted(b), ["AAAA", "AACC", "ACAC", "CAAC", "CCAA", "GGGG"])
        check(name + " TopCells 6 total counts", sum(t.values()), sum(total.values()))
        # three expected cells, 99th percentile of them is the largest total (500); cells with at least a tenth of it
        b, t = call(star, raw, out / name / "cr22", ["CellRanger2.2", "3", "0.99", "10"])
        check(name + " CellRanger2.2 barcodes", sorted(b), ["AAAA", "AACC", "GGGG"])

    # TopCells N reports N barcodes, and any that tie with the N-th
    raw = out / "ties" / "raw"
    raw.mkdir(parents=True)
    (raw / "features.tsv").write_text("G1\tG1\tGene Expression\n")
    (raw / "barcodes.tsv").write_text("AAAA\nCCCC\nGGGG\nTTTT\n")
    (raw / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n%\n1 4 4\n1 1 10\n1 2 7\n1 3 7\n1 4 1\n")
    b, _ = call(star, raw, out / "ties" / "top1", ["TopCells", "1"])
    check("ties TopCells 1 barcodes", sorted(b), ["AAAA"])
    b, _ = call(star, raw, out / "ties" / "top2", ["TopCells", "2"])
    check("ties TopCells 2 barcodes (two tied at the second place)", sorted(b), ["AAAA", "CCCC", "GGGG"])
    b, _ = call(star, raw, out / "ties" / "top9", ["TopCells", "9"])
    check("ties TopCells 9 barcodes (more than there are)", sorted(b), ["AAAA", "CCCC", "GGGG", "TTTT"])
    # the number of cells has to be a whole decimal integer >= 1 that fits the
    # counter: nothing may be read from a numeric prefix (1e3 is not 1)
    for i, arg in enumerate(["0", "-1", "+2", "2junk", "1.5", "1e3", "0x2", "", "4294967296", "99999999999999999999999"]):
        bad = out / "ties" / ("rejected_%d" % i)
        bad.mkdir()
        proc = subprocess.run([str(star), "--runMode", "soloCellFiltering", str(raw), str(bad) + "/", "--outFileNamePrefix", str(bad) + "/",
                               "--soloCellFilter", "TopCells", arg], capture_output=True, text=True)
        log = bad / "Log.out"
        said = proc.stdout + proc.stderr + (log.read_text() if log.exists() else "")
        check("TopCells '%s' is rejected" % arg,
              proc.returncode != 0 and ("integer >= 1" in said or "number of cells not specified" in said)
              and not (bad / "barcodes.tsv").exists(), True)
    # accepted: leading zeros, and the largest value of the counter
    b, _ = call(star, raw, out / "ties" / "top002", ["TopCells", "002"])
    check("ties TopCells 002 barcodes", sorted(b), ["AAAA", "CCCC", "GGGG"])
    b, _ = call(star, raw, out / "ties" / "topmax", ["TopCells", "4294967295"])
    check("ties TopCells 4294967295 barcodes", sorted(b), ["AAAA", "CCCC", "GGGG", "TTTT"])

    # a matrix in which one barcode has counts: before the loader fix the
    # number of cells was 0 and no filtering was done
    raw = out / "single" / "raw"
    raw.mkdir(parents=True)
    (raw / "features.tsv").write_text("G1\tG1\tGene Expression\nG2\tG2\tGene Expression\n")
    (raw / "barcodes.tsv").write_text("AAAA\nCCCC\nGGGG\n")
    (raw / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n%\n2 3 2\n1 2 40\n2 2 2\n")
    for label, cell_filter in (("top1", ["TopCells", "1"]), ("top5", ["TopCells", "5"]), ("cr22", ["CellRanger2.2", "3", "0.99", "10"])):
        b, t = call(star, raw, out / "single" / label, cell_filter)
        check("single barcode %s barcodes" % " ".join(cell_filter), b, ["CCCC"])
        check("single barcode %s total" % " ".join(cell_filter), t, {"CCCC": 42})

    if failures:
        print("FAILED: %d checks" % len(failures))
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()

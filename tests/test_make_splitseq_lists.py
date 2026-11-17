#!/usr/bin/env python3
"""Checks scripts/make_splitseq_lists.py on small generated lists and on the
bundled v3 lists.

  python3 tests/test_make_splitseq_lists.py --out tests/splitseq_lists_output

--tool PATH runs the checks on another copy of the tool.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "scripts" / "make_splitseq_lists.py"
BUNDLE = ROOT / "share" / "star-suite" / "barcode_lists" / "splitseq_v3"

# round 3, round 2, round 1 columns; round 1 has two wells with an oligo-dT (T)
# and a random-hexamer (R) entry each, one further entry (X) mapped to well 1
# and one entry with no replacement line
R3 = ["AAAA", "CCCC", "GGGG"]
R2 = ["ACAC", "CACA"]
R1 = ["ATAT", "TATA", "AGAG", "GAGA", "TTAA", "GGCC"]
REPLACE = [("TATA", "ATAT"), ("GAGA", "AGAG"), ("TTAA", "ATAT")]


def run(argv, ok=True):
    proc = subprocess.run([sys.executable, str(run.tool)] + [str(x) for x in argv], capture_output=True, text=True)
    assert (proc.returncode == 0) == ok, proc.stdout + proc.stderr
    return proc


def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--tool", type=Path, default=TOOL)
    args = ap.parse_args()
    run.tool = args.tool.resolve()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    onlist, replace = out / "onlist.txt", out / "replace.txt"
    n = max(len(R3), len(R2), len(R1))
    pad = lambda c: c + ["-"] * (n - len(c))
    onlist.write_text("".join("%s %s %s\n" % t for t in zip(pad(R3), pad(R2), pad(R1))))
    replace.write_text("".join("%s *%s\n" % t for t in REPLACE))

    # all entries
    d = out / "all"
    run(["--onlist", onlist, "--replace", replace, "--out-dir", d])
    assert rows(d / "round1.tsv") == [["ATAT", "ATAT"], ["TATA", "ATAT"], ["AGAG", "AGAG"], ["GAGA", "AGAG"],
                                      ["TTAA", "ATAT"], ["GGCC", "GGCC"]]
    assert (d / "round2.tsv").read_text() == "ACAC\nCACA\n"
    assert (d / "round3.tsv").read_text() == "AAAA\nCCCC\nGGGG\n"
    cross = {(r[0], r[1]): r[2:] for r in rows(d / "crosswalk.tsv")[1:]}
    assert cross[("round1", "ATAT")] == ["T", "A1", "1", "ATAT"]
    assert cross[("round1", "TATA")] == ["R", "A1", "1", "ATAT"]
    assert cross[("round1", "GAGA")] == ["R", "A2", "2", "AGAG"]
    assert cross[("round1", "TTAA")] == ["X", "A1", "1", "ATAT"]
    assert cross[("round1", "GGCC")] == ["X", "", "0", "GGCC"]
    assert cross[("round3", "GGGG")] == ["L", "A3", "3", "GGGG"]
    summary = json.loads((d / "summary.json").read_text())
    assert summary["round1_kept_primers"] == {"T": 2, "R": 2, "X": 2}
    assert summary["outputs"]["round1.tsv"] == hashlib.sha256((d / "round1.tsv").read_bytes()).hexdigest()

    # the first well only; the entry with no replacement line is left out
    d = out / "one_well"
    run(["--onlist", onlist, "--replace", replace, "--wells", "1", "--out-dir", d])
    assert rows(d / "round1.tsv") == [["ATAT", "ATAT"], ["TATA", "ATAT"], ["TTAA", "ATAT"]]

    # wells seen in a list of observed sequences
    seen = out / "seen.txt"
    seen.write_text("GAGA\n")
    d = out / "seen"
    run(["--onlist", onlist, "--replace", replace, "--restrict-round1", seen, "--out-dir", d])
    assert rows(d / "round1.tsv") == [["AGAG", "AGAG"], ["GAGA", "AGAG"]]

    # well labels of cell barcodes (canonical round 1 + round 2 + round 3)
    barcodes = out / "barcodes.tsv"
    barcodes.write_text("ATATACACAAAA\nAGAGCACAGGGG\nGGCCACACCCCC\n")
    d = out / "labels"
    run(["--onlist", onlist, "--replace", replace, "--label-barcodes", barcodes, "--out-dir", d])
    assert rows(d / "barcode_wells.tsv")[1:] == [["ATATACACAAAA", "A1", "A1", "A1", "01_01_01"],
                                                 ["AGAGCACAGGGG", "A2", "A2", "A3", "02_02_03"],
                                                 ["GGCCACACCCCC", "", "A1", "A2", ""]]

    # other column order
    swapped = out / "onlist_swapped.txt"
    swapped.write_text("".join("%s %s %s\n" % t for t in zip(pad(R1), pad(R2), pad(R3))))
    d = out / "swapped"
    run(["--onlist", swapped, "--replace", replace, "--rounds", "1,2,3", "--out-dir", d])
    assert (d / "round1.tsv").read_text() == (out / "all" / "round1.tsv").read_text()

    # errors
    bad = out / "bad_replace.txt"
    bad.write_text("TATA *ATAT\nTATA *AGAG\n")
    assert "source listed twice" in run(["--onlist", onlist, "--replace", bad, "--out-dir", out / "e1"], ok=False).stderr
    bad.write_text("NNNN *ATAT\n")
    assert "not in the round-1 column" in run(["--onlist", onlist, "--replace", bad, "--out-dir", out / "e2"], ok=False).stderr
    bad.write_text("TATA *ATAT\nATAT *AGAG\n")
    assert "both source and target" in run(["--onlist", onlist, "--replace", bad, "--out-dir", out / "e3"], ok=False).stderr
    dup = out / "dup_onlist.txt"
    dup.write_text("AAAA ACAC ATAT\nAAAA CACA TATA\n")
    assert "duplicate entry" in run(["--onlist", dup, "--replace", replace, "--out-dir", out / "e4"], ok=False).stderr
    assert "not both" in run(["--onlist", onlist, "--replace", replace, "--wells", "1", "--restrict-round1", seen,
                              "--out-dir", out / "e5"], ok=False).stderr
    for name in ("e1", "e2", "e3", "e4", "e5"):
        assert not (out / name).exists(), name

    # an output directory that is used again holds one set of files that its summary.json describes
    def snapshot(d):
        return {p.name: (p.read_bytes() if p.is_file() else "<directory>") for p in sorted(d.iterdir())}

    def consistent(d, other=()):
        summary = json.loads((d / "summary.json").read_text())
        assert set(snapshot(d)) == set(summary["outputs"]) | {"summary.json"} | set(other), sorted(snapshot(d))
        for name, digest in summary["outputs"].items():
            assert hashlib.sha256((d / name).read_bytes()).hexdigest() == digest, name
        return summary

    d = out / "reuse"
    run(["--onlist", onlist, "--replace", replace, "--wells", "2", "--label-barcodes", barcodes, "--out-dir", d])
    assert consistent(d)["labelled"] == {"barcodes": 3, "without_well": 1}
    assert (d / "barcode_wells.tsv").exists()
    # one well and no labels into the same directory: the labels of the earlier run are removed
    (d / "notes.txt").write_text("not written by the tool\n")
    run(["--onlist", onlist, "--replace", replace, "--wells", "1", "--out-dir", d])
    summary = consistent(d, other=["notes.txt"])
    assert summary["labelled"] is None and summary["wells"] == 1
    assert not (d / "barcode_wells.tsv").exists()
    assert rows(d / "round1.tsv") == [["ATAT", "ATAT"], ["TATA", "ATAT"], ["TTAA", "ATAT"]]
    # a run that fails while its files are made (a cell barcode of the wrong length) leaves the directory as it was
    before = snapshot(d)
    bad_cells = out / "bad_cells.tsv"
    bad_cells.write_text("ATATACACAAAA\nA\n")
    proc = run(["--onlist", onlist, "--replace", replace, "--wells", "2", "--label-barcodes", bad_cells, "--out-dir", d], ok=False)
    assert "unexpected length" in proc.stderr
    assert snapshot(d) == before
    consistent(d, other=["notes.txt"])
    # the same failing run into a new directory leaves no directory
    run(["--onlist", onlist, "--replace", replace, "--wells", "2", "--label-barcodes", bad_cells, "--out-dir", out / "e6"], ok=False)
    assert not (out / "e6").exists()
    # a label file that cannot be read, and an output path that is a file
    run(["--onlist", onlist, "--replace", replace, "--label-barcodes", out / "no_such_file.tsv", "--out-dir", d], ok=False)
    assert snapshot(d) == before
    assert "not a directory" in run(["--onlist", onlist, "--replace", replace, "--out-dir", onlist], ok=False).stderr
    assert onlist.read_text().startswith("AAAA ACAC ATAT")
    # labels again: the set is complete again
    run(["--onlist", onlist, "--replace", replace, "--wells", "2", "--label-barcodes", barcodes, "--out-dir", d])
    assert consistent(d, other=["notes.txt"])["labelled"] == {"barcodes": 3, "without_well": 1}
    assert (d / "notes.txt").read_text() == "not written by the tool\n"
    # a move into place that fails part way (round3.tsv cannot be replaced: it is a directory with a file in it):
    # the run fails, says so, and leaves no summary.json; the files are a mix of the two runs
    before = snapshot(d)
    (d / "round3.tsv").unlink()
    (d / "round3.tsv").mkdir()
    (d / "round3.tsv" / "in_the_way").write_text("x\n")
    proc = run(["--onlist", onlist, "--replace", replace, "--wells", "1", "--out-dir", d], ok=False)
    assert "no summary.json" in proc.stderr and "Traceback" not in proc.stderr, proc.stderr
    after = snapshot(d)
    assert "summary.json" not in after
    assert not [n for n in after if n.startswith(".make_splitseq_lists.")], sorted(after)
    assert after["round1.tsv"] != before["round1.tsv"] and after["round3.tsv"] == "<directory>"
    assert after["notes.txt"] == before["notes.txt"]
    # with the cause removed the next run completes the set
    (d / "round3.tsv" / "in_the_way").unlink()
    (d / "round3.tsv").rmdir()
    run(["--onlist", onlist, "--replace", replace, "--wells", "1", "--out-dir", d])
    assert consistent(d, other=["notes.txt"])["wells"] == 1

    # the bundled lists are what the tool writes from the bundled source files
    for wells in (96, 48):
        d = out / ("bundle_%d" % wells)
        run(["--onlist", BUNDLE / "source" / "splitseqv3_barcodes.txt", "--replace", BUNDLE / "source" / "splitseqv3_replace.txt",
             "--wells", wells, "--out-dir", d])
        for name in ("round1.tsv", "round2.tsv", "round3.tsv", "crosswalk.tsv", "summary.json"):
            assert (d / name).read_bytes() == (BUNDLE / ("%d_wells" % wells) / name).read_bytes(), (wells, name)
    print("PASS make_splitseq_lists: lists, wells, labels, column order, errors, reuse of an output directory, bundled files")


if __name__ == "__main__":
    main()

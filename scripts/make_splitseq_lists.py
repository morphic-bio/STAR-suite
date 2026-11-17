#!/usr/bin/env python3
"""Allowlists for --soloSplitseqWhitelist from a kallisto-bustools on-list and
replacement list, and well labels for the resulting cell barcodes.

The kallisto-bustools documentation (pachterlab/kb_docs, BSD-2-Clause)
publishes, per SPLiT-seq chemistry, an on-list and a replacement list:

  on-list           whitespace-separated columns, one per barcode segment in
                    read order, shorter columns padded with "-"
  replacement list  "SOURCE *TARGET" per line: a source round-1 barcode is
                    counted as its target (random-hexamer barcodes into the
                    oligo-dT barcode of the same well)

  make_splitseq_lists.py --onlist barcodes.txt --replace replace.txt --out-dir DIR
      [--rounds 3,2,1] [--wells N] [--restrict-round1 observed.txt]
      [--label-barcodes barcodes.tsv]

--rounds names the round of each on-list column (default: column 1 is round 3,
column 2 round 2, column 3 round 1, the order of the segments in the barcode
read of the v3 chemistry).

Written to DIR:
  round1.tsv      "sequence<TAB>canonical" for every round-1 entry; canonical
                  is the replacement target, or the sequence itself
  round2.tsv, round3.tsv   the on-list columns, one sequence per line
  crosswalk.tsv   round, sequence, primer, well, well_index, canonical
  summary.json    inputs with SHA256, counts, SHA256 of the outputs

Pass round1.tsv round2.tsv round3.tsv to --soloSplitseqWhitelist in this
order, with --soloSplitseqCBposition in the same round order.

DIR may hold the files of an earlier run. The new files are written to a
temporary directory inside DIR and moved into place only when all of them are
complete, so a run that fails while it makes or checks them leaves DIR as it
was. The move removes summary.json of the earlier run first and writes the new
one last; barcode_wells.tsv of an earlier run is removed when --label-barcodes
is not given. The move is file by file, not one step: if it fails or is
interrupted, DIR can hold old and new files side by side and has no
summary.json. Use DIR only when summary.json is there; otherwise run the tool
again. Other files in DIR are not touched.

Wells. For rounds 2 and 3 the well index is the line number in the column. For
round 1 it is the order in which targets first appear in the replacement list;
a source has its target's index. primer is read off the lists: T a target, R
the first source of a target, X any further source or an entry with no
replacement line (well index 0). The files do not state that these orders are
plate orders (row-major: A1 = 1, A12 = 12, B1 = 13); docs/SPLITSEQ_RNA.md
records the check against a public library.

--wells N keeps the round-1 entries of wells 1..N (a 48-well kit: 48); entries
with no replacement line are then left out. --restrict-round1 FILE keeps the
wells that have at least one sequence in FILE instead.

--label-barcodes FILE reads cell barcodes as the mode writes them (canonical
round 1, round 2 and round 3 sequences concatenated) and writes
barcode_wells.tsv: barcode, the three well names and the label
"%02d_%02d_%02d" of the three well indices.
"""
import argparse, gzip, hashlib, json, os, shutil, sys, tempfile

# every file this tool writes into --out-dir; summary.json is the record of a complete set
LISTS = ["round1.tsv", "round2.tsv", "round3.tsv", "crosswalk.tsv"]
LABELS = "barcode_wells.tsv"
SUMMARY = "summary.json"


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def well_name(i):
    return "" if i <= 0 else "ABCDEFGH"[(i - 1) // 12] + str((i - 1) % 12 + 1)


def fail(msg):
    sys.exit("make_splitseq_lists.py: " + msg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--onlist", required=True)
    ap.add_argument("--replace", required=True)
    ap.add_argument("--rounds", default="3,2,1")
    ap.add_argument("--wells", type=int, default=0)
    ap.add_argument("--restrict-round1")
    ap.add_argument("--label-barcodes")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    if a.wells and a.restrict_round1:
        fail("give --wells or --restrict-round1, not both")
    rounds = ["round" + r for r in a.rounds.split(",")]
    if sorted(rounds) != ["round1", "round2", "round3"]:
        fail("--rounds must be a permutation of 1,2,3")
    cols = {r: [] for r in rounds}
    for n, line in enumerate(open(a.onlist), 1):
        p = line.split()
        if not p:
            continue
        if len(p) != len(rounds):
            fail("%s line %d has %d columns, expected %d" % (a.onlist, n, len(p), len(rounds)))
        for r, w in zip(rounds, p):
            if w != "-":
                if not w or set(w) - set("ACGT"):
                    fail("%s line %d: not an A/C/G/T sequence: %s" % (a.onlist, n, w))
                cols[r].append(w)
    for r in rounds:
        if not cols[r]:
            fail("empty column for " + r)
        if len(set(cols[r])) != len(cols[r]):
            fail("duplicate entry in the %s column" % r)
        if len({len(w) for w in cols[r]}) != 1:
            fail("entries of different lengths in the %s column" % r)
    r1 = cols["round1"]
    r1set = set(r1)
    target, index, primer = {}, {}, {}
    for n, line in enumerate(open(a.replace), 1):
        p = line.split()
        if not p:
            continue
        if len(p) != 2:
            fail("%s line %d: expected SOURCE *TARGET" % (a.replace, n))
        s, t = p[0], p[1].lstrip("*")
        if s not in r1set or t not in r1set:
            fail("%s line %d: entry not in the round-1 column" % (a.replace, n))
        if s in target:
            fail("%s line %d: source listed twice: %s" % (a.replace, n, s))
        if t not in index:
            index[t] = len(set(target.values())) + 1
            primer[t] = "T"
            primer[s] = "R"
        else:
            primer[s] = "X"
        target[s] = t
        index[s] = index[t]
    both = set(target) & set(target.values())
    if both:
        fail("a barcode is both source and target: " + sorted(both)[0])
    keep = list(r1)
    if a.wells:
        if a.wells < 1 or a.wells > len(set(target.values())):
            fail("--wells must be between 1 and %d" % len(set(target.values())))
        keep = [s for s in r1 if 1 <= index.get(s, 0) <= a.wells]
    elif a.restrict_round1:
        seen = {l.strip() for l in open(a.restrict_round1) if l.strip()}
        wells = {index[s] for s in r1 if s in seen and s in index}
        if not wells:
            fail("no sequence of --restrict-round1 is a round-1 entry with a well")
        keep = [s for s in r1 if index.get(s, 0) in wells]
    if os.path.exists(a.out_dir) and not os.path.isdir(a.out_dir):
        fail("--out-dir is not a directory: " + a.out_dir)
    created = not os.path.isdir(a.out_dir)
    os.makedirs(a.out_dir, exist_ok=True)
    stage = tempfile.mkdtemp(prefix=".make_splitseq_lists.", dir=a.out_dir)
    done = False
    try:
        summary = write_set(a, stage, rounds, cols, r1, keep, target, index, primer)
        try:
            publish(stage, a.out_dir, list(summary["outputs"]))
        except OSError as e:
            fail("could not move the new files into %s (%s). The directory has no summary.json and may hold "
                 "old and new files side by side: do not use it, correct the cause and run again" % (a.out_dir, e))
        done = True
    finally:
        shutil.rmtree(stage, ignore_errors=True)
        if created and not done:
            try:
                os.rmdir(a.out_dir)
            except OSError:
                pass
    print(json.dumps(summary, indent=1))


def publish(stage, out_dir, outputs):
    """Replace the files of an earlier run in out_dir by the complete set in stage.

    summary.json goes first and comes back last, so a publish that fails or is
    interrupted leaves no summary; files of an earlier run that this run does
    not write are removed. There is no rollback: after a failure out_dir can
    hold old and new files side by side.
    """
    for f in [SUMMARY] + [f for f in LISTS + [LABELS] if f not in outputs]:
        try:
            os.remove(os.path.join(out_dir, f))
        except FileNotFoundError:
            pass
    for f in outputs + [SUMMARY]:
        os.replace(os.path.join(stage, f), os.path.join(out_dir, f))


def write_set(a, out_dir, rounds, cols, r1, keep, target, index, primer):
    """Write the lists, the optional labels and summary.json into out_dir; return the summary."""
    rows = []
    with open(os.path.join(out_dir, "round1.tsv"), "w") as m:
        for s in keep:
            c = target.get(s, s)
            m.write("%s\t%s\n" % (s, c))
            i = index.get(s, 0)
            rows.append(("round1", s, primer.get(s, "X"), well_name(i), i, c))
    for r in ("round2", "round3"):
        with open(os.path.join(out_dir, r + ".tsv"), "w") as o:
            for i, s in enumerate(cols[r], 1):
                o.write(s + "\n")
                rows.append((r, s, "L", well_name(i), i, s))
    with open(os.path.join(out_dir, "crosswalk.tsv"), "w") as o:
        o.write("round\tsequence\tprimer\twell\twell_index\tcanonical\n")
        for row in rows:
            o.write("\t".join(str(x) for x in row) + "\n")
    outputs = list(LISTS)
    labelled = None
    if a.label_barcodes:
        widx = {r: {} for r in ("round1", "round2", "round3")}
        for r, s, _, _, i, c in rows:
            if r != "round1" or s == c:
                widx[r][s] = i
        ln = {r: len(cols[r][0]) for r in widx}
        op = gzip.open if a.label_barcodes.endswith(".gz") else open
        labelled = {"barcodes": 0, "without_well": 0}
        with op(a.label_barcodes, "rt") as f, open(os.path.join(out_dir, LABELS), "w") as o:
            o.write("barcode\tround1_well\tround2_well\tround3_well\tlabel\n")
            for line in f:
                b = line.strip()
                if not b:
                    continue
                if len(b) != ln["round1"] + ln["round2"] + ln["round3"]:
                    fail("barcode of unexpected length: " + b)
                w = (b[:ln["round1"]], b[ln["round1"]:ln["round1"] + ln["round2"]], b[ln["round1"] + ln["round2"]:])
                idx = [widx[r].get(x, 0) for r, x in zip(("round1", "round2", "round3"), w)]
                labelled["barcodes"] += 1
                labelled["without_well"] += int(0 in idx)
                o.write("%s\t%s\t%s\t%s\t%s\n" % (b, well_name(idx[0]), well_name(idx[1]), well_name(idx[2]),
                                                    "%02d_%02d_%02d" % tuple(idx) if 0 not in idx else ""))
        outputs.append(LABELS)
    kept = [r for r in rows if r[0] == "round1"]
    summary = {"onlist": {"file": os.path.basename(a.onlist), "sha256": sha(a.onlist)},
               "replace": {"file": os.path.basename(a.replace), "sha256": sha(a.replace)},
               "column_rounds": rounds, "column_sizes": {r: len(cols[r]) for r in rounds},
               "replace_lines": len(target), "replace_targets": len(set(target.values())),
               "round1_without_replacement_line": len([s for s in r1 if s not in index]),
               "wells": a.wells or None, "restricted_to": os.path.basename(a.restrict_round1) if a.restrict_round1 else None,
               "round1_kept": len(kept),
               "round1_kept_primers": {t: sum(r[2] == t for r in kept) for t in ("T", "R", "X")},
               "round1_kept_wells": len({r[4] for r in kept if r[4]}),
               "labelled": labelled,
               "outputs": {f: sha(os.path.join(out_dir, f)) for f in outputs}}
    with open(os.path.join(out_dir, SUMMARY), "w") as o:
        json.dump(summary, o, indent=1)
        o.write("\n")
    return summary


if __name__ == "__main__":
    main()

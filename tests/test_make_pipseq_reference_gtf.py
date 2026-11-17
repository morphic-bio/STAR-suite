#!/usr/bin/env python3
"""Checks scripts/make_pipseq_reference_gtf.py on a small generated GTF.

  python3 tests/test_make_pipseq_reference_gtf.py --out tests/pipseq_reference_gtf_output
"""
import argparse
import gzip
import json
from pathlib import Path
import subprocess
import sys

TOOL = Path(__file__).resolve().parent.parent / "scripts" / "make_pipseq_reference_gtf.py"

RT = "readthrough_transcript"
# gene id, gene type, class it is dropped in (None = kept), tags of its transcripts
GENES = [
    ("ENSG01.1", "protein_coding", None, [None]),
    ("ENSG02.3", "lncRNA", None, [None]),
    ("ENSG03.1", "TEC", None, [None]),
    ("ENSG04.2", "Mt_tRNA", None, [None]),
    ("ENSG05.1", "IG_V_gene", None, [None]),
    ("ENSG06.1", "ribozyme", None, [None]),
    ("ENSG07.1", "processed_pseudogene", "pseudogene", [None]),
    ("ENSG08.4", "transcribed_unprocessed_pseudogene", "pseudogene", [None]),
    ("ENSG09.1", "IG_V_pseudogene", "pseudogene", [None]),
    ("ENSG10.1", "rRNA_pseudogene", "pseudogene", [None]),
    ("ENSG11.1", "misc_RNA", "shortRNA", [None]),
    ("ENSG12.1", "snRNA", "shortRNA", [None]),
    ("ENSG13.1", "miRNA", "shortRNA", [None]),
    ("ENSG14.1", "snoRNA", "shortRNA", [None]),
    ("ENSG15.1", "scaRNA", "shortRNA", [None]),
    ("ENSG16.1", "sRNA", "shortRNA", [None]),
    ("ENSG17.1", "scRNA", "shortRNA", [None]),
    ("ENSG18.1", "rRNA", "rRNA", [None]),
    ("ENSG19.1", "Mt_rRNA", "rRNA", [None]),
    # every transcript is read-through: the gene goes
    ("ENSG20.1", "protein_coding", "readthrough", [RT, RT]),
    # one of two transcripts is read-through: the gene stays, that transcript goes
    ("ENSG21.2", "lncRNA", None, [None, RT]),
    # a different tag with a similar name is not a read-through transcript
    ("ENSG22.1", "protein_coding", None, ["stop_codon_readthrough"]),
    # a pseudogene with read-through transcripts is counted as a pseudogene
    ("ENSG23.1", "unprocessed_pseudogene", "pseudogene", [RT]),
]
BIOTYPES = {"pseudogene", "shortRNA", "rRNA"}


def gtf_text(type_attribute):
    lines = ["##description: generated for a test\n", "##provider: test\n"]
    for i, (gid, gtype, _, tags) in enumerate(GENES):
        start = 1000 * (i + 1)
        attrs = 'gene_id "%s"; %s "%s"; gene_name "G%d";' % (gid, type_attribute, gtype, i)
        lines.append("chr1\tTEST\tgene\t%d\t%d\t.\t+\t.\t%s\n" % (start, start + 500, attrs))
        for j, tag in enumerate(tags):
            tattrs = attrs + ' transcript_id "T%d_%d.1";' % (i, j)
            if tag:
                tattrs += ' tag "basic"; tag "%s";' % tag
            lines.append("chr1\tTEST\ttranscript\t%d\t%d\t.\t+\t.\t%s\n" % (start, start + 500, tattrs))
            lines.append("chr1\tTEST\texon\t%d\t%d\t.\t+\t.\t%s\n" % (start, start + 200, tattrs))
            lines.append("chr1\tTEST\texon\t%d\t%d\t.\t+\t.\t%s\n" % (start + 300, start + 500, tattrs))
    return "".join(lines)


def expected(drop):
    """Kept gene ids and the number of records written for a set of dropped classes."""
    ids, n = [], 0
    for gid, _, cls, tags in GENES:
        if cls in drop:
            continue
        ids.append(gid)
        kept_transcripts = [t for t in tags if not ("readthrough" in drop and t == RT)]
        n += 1 + 3 * len(kept_transcripts)
    return ids, n


def run(argv, expect):
    result = subprocess.run([sys.executable, str(TOOL)] + [str(v) for v in argv],
                            capture_output=True, text=True)
    if result.returncode != expect:
        raise AssertionError("exit %d, expected %d\n%s\n%s" % (result.returncode, expect, result.stdout, result.stderr))
    return result.stdout


def genes_in(path):
    ids, records = [], 0
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            continue
        records += 1
        fields = line.split("\t")
        if fields[2] == "gene":
            ids.append(fields[8].split('"')[1])
    return ids, records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)

    plain = args.out / "gencode_style.gtf"
    plain.write_text(gtf_text("gene_type"))
    packed = args.out / "ensembl_style.gtf.gz"
    with gzip.open(packed, "wt") as f:
        f.write(gtf_text("gene_biotype"))

    checks = {}

    # default: the three biotype classes and read-through
    out = args.out / "default.gtf"
    run(["--gtf", plain, "--out-gtf", out, "--report", args.out / "default.json"], 0)
    assert genes_in(out) == expected(BIOTYPES | {"readthrough"}), genes_in(out)
    text = out.read_text()
    assert text.startswith("##description: generated for a test\n")
    assert 'transcript_id "T20_0.1"' in text and 'transcript_id "T20_1.1"' not in text   # ENSG21: one of two
    assert 'gene_id "ENSG20.1"' not in text and 'gene_id "ENSG22.1"' in text
    report = json.loads((args.out / "default.json").read_text())
    assert report["genes_in_gtf"] == len(GENES)
    assert report["genes_dropped_by_class"] == {"pseudogene": 5, "shortRNA": 7, "rRNA": 2, "readthrough": 1}, report
    assert report["readthrough_transcripts_dropped_in_kept_genes"] == 1, report
    assert report["by_gene_type"]["processed_pseudogene"] == {"dropped": 1}
    assert report["by_gene_type"]["Mt_tRNA"] == {"kept": 1}
    assert report["by_gene_type"]["protein_coding"] == {"kept": 2, "dropped": 1}
    checks["default"] = "ok"

    # the documented list: read-through genes and transcripts stay
    out = args.out / "documented.gtf"
    run(["--gtf", plain, "--out-gtf", out, "--drop", "pseudogene,shortRNA,rRNA", "--report", args.out / "documented.json"], 0)
    assert genes_in(out) == expected(BIOTYPES), genes_in(out)
    report = json.loads((args.out / "documented.json").read_text())
    assert report["genes_dropped_by_class"] == {"pseudogene": 5, "shortRNA": 7, "rRNA": 2}, report
    assert report["readthrough_transcripts_dropped_in_kept_genes"] == 0
    checks["documented_list"] = "ok"

    # Ensembl attribute name and gzip input give the same list
    out = args.out / "ensembl.gtf"
    run(["--gtf", packed, "--out-gtf", out], 0)
    assert genes_in(out) == expected(BIOTYPES | {"readthrough"})
    checks["gene_biotype_gzip"] = "ok"

    # one class only
    out = args.out / "rrna_only.gtf"
    run(["--gtf", plain, "--out-gtf", out, "--drop", "rRNA"], 0)
    assert genes_in(out) == expected({"rRNA"})
    out = args.out / "readthrough_only.gtf"
    run(["--gtf", plain, "--out-gtf", out, "--drop", "readthrough"], 0)
    ids, _ = genes_in(out)
    # ENSG23 is a pseudogene made of read-through transcripts: without the pseudogene class it goes as read-through
    assert ids == [g for g, _, cls, tags in GENES if not all(t == RT for t in tags)], ids
    checks["drop_one_class"] = "ok"

    # nothing dropped: records equal the input
    out = args.out / "none.gtf"
    run(["--gtf", plain, "--out-gtf", out, "--drop", "none"], 0)
    assert out.read_text() == plain.read_text()
    checks["drop_none"] = "ok"

    # unknown class is an error
    result = subprocess.run([sys.executable, str(TOOL), "--gtf", str(plain), "--out-gtf",
                             str(args.out / "bad.gtf"), "--drop", "tRNA"], capture_output=True, text=True)
    assert result.returncode != 0 and "unknown class" in result.stderr, result.stderr
    checks["unknown_class"] = "ok"

    # features list: identical, then one gene missing
    kept_default = expected(BIOTYPES | {"readthrough"})[0]
    features = args.out / "features.tsv.gz"
    with gzip.open(features, "wt") as f:
        for i, gid in enumerate(kept_default):
            f.write("%s\tG%d\tGene Expression\n" % (gid, i))
    stdout = run(["--gtf", plain, "--out-gtf", args.out / "checked.gtf", "--check-features", features,
                  "--report", args.out / "checked.json"], 0)
    assert json.loads(stdout)["identical_to_features_list"] is True
    short = args.out / "features_short.tsv"
    short.write_text("".join("%s\tx\tGene Expression\n" % g for g in kept_default[1:]))
    run(["--gtf", plain, "--out-gtf", args.out / "checked_short.gtf", "--check-features", short,
         "--report", args.out / "checked_short.json"], 1)
    report = json.loads((args.out / "checked_short.json").read_text())
    assert report["check_features"]["kept_not_listed"] == 1
    assert report["check_features"]["examples_kept_not_listed"] == [kept_default[0]]
    checks["check_features"] = "ok"

    (args.out / "summary.json").write_text(json.dumps(checks, indent=1) + "\n")
    print(json.dumps(checks))


if __name__ == "__main__":
    main()

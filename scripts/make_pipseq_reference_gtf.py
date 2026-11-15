#!/usr/bin/env python3
"""Annotation for --soloPipseqMode runs: a GENCODE GTF without pseudogene,
short RNA, rRNA and read-through genes.

Illumina's documentation for this chemistry names GENCODE v44 (comprehensive
annotation on the reference chromosomes, gencode.v44.annotation.gtf.gz) as the
default annotation and says the pipeline ignores the biotypes "pseudogene",
"shortRNA" and "rRNA". It does not list the gene types in each class. The
assignment used here was read from the feature list of a public deposit
(GEO GSE295097, sample GSM8942285), where it reproduces the list exactly:

  pseudogene   every gene type whose name contains "pseudogene"
  shortRNA     misc_RNA, snRNA, miRNA, snoRNA, scaRNA, sRNA, scRNA
  rRNA         rRNA, Mt_rRNA

A fourth class is dropped by default because of how STAR counts:

  readthrough  transcripts tagged "readthrough_transcript", and genes left
               without a transcript

A read-through gene joins neighbouring genes and shares their exons. STAR
leaves out a read that fits two genes, so with such a gene in the annotation
its neighbours lose most of their reads. The documented list keeps these
genes; they hold about 0.2% of the counts of the deposit named above.

Everything else is kept, including lncRNA, TEC, Mt_tRNA and the IG and TR
genes. Records of kept genes and kept transcripts are written unchanged:
nothing is renamed, and a gene record keeps its coordinates when one of its
transcripts is dropped. Reads that align inside a dropped gene are counted for
a kept gene that contains them, or not at all.

  make_pipseq_reference_gtf.py --gtf gencode.v44.annotation.gtf.gz --out-gtf genes.gtf
      [--drop pseudogene,shortRNA,rRNA,readthrough] [--report report.json]
      [--check-features features.tsv.gz]

--drop pseudogene,shortRNA,rRNA writes the documented list itself.
--check-features compares the kept gene identifiers with the first column of a
features file (plain or gzip) and makes the exit status 1 if the lists differ.
Build the STAR index from the genome FASTA and the written GTF:

  STAR --runMode genomeGenerate --genomeDir index --genomeFastaFiles genome.fa \
       --sjdbGTFfile genes.gtf --sjdbOverhang 100

Standard library only. GENCODE (gene_type) and Ensembl (gene_biotype)
attribute names are both read; the read-through tag is GENCODE's.
"""
import argparse
import collections
import gzip
import hashlib
import json
import re
import sys

SHORT_RNA = {"misc_RNA", "snRNA", "miRNA", "snoRNA", "scaRNA", "sRNA", "scRNA"}
RRNA = {"rRNA", "Mt_rRNA"}
BIOTYPE_CLASSES = ("pseudogene", "shortRNA", "rRNA")
CLASSES = BIOTYPE_CLASSES + ("readthrough",)
READTHROUGH_TAG = ' tag "readthrough_transcript"'

GENE_ID = re.compile(r'gene_id "([^"]+)"')
TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
GENE_TYPE = re.compile(r'gene_(?:type|biotype) "([^"]+)"')


def gene_class(gene_type):
    """The documented biotype class a gene type falls in, or None."""
    if "pseudogene" in gene_type:
        return "pseudogene"
    if gene_type in SHORT_RNA:
        return "shortRNA"
    if gene_type in RRNA:
        return "rRNA"
    return None


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def records(path):
    """Yield (line, fields) for every non-comment line; comments as (line, None)."""
    with open_text(path) as src:
        for line in src:
            if line.startswith("#"):
                yield line, None
                continue
            fields = line.split("\t", 9)
            if len(fields) < 9:
                sys.exit("not a GTF record: %r" % line[:80])
            yield line, fields


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--drop", default=",".join(CLASSES),
                    help="classes to drop, comma-separated; 'none' keeps every gene")
    ap.add_argument("--report")
    ap.add_argument("--check-features")
    a = ap.parse_args()

    drop = set() if a.drop.lower() == "none" else set(a.drop.split(","))
    unknown = drop - set(CLASSES)
    if unknown:
        sys.exit("unknown class in --drop: %s (choose from %s)" % (",".join(sorted(unknown)), ",".join(CLASSES)))

    # First pass: which transcripts are read-through, and which genes consist of them only.
    transcripts = collections.Counter()
    readthrough = collections.Counter()
    dropped_transcripts = set()
    if "readthrough" in drop:
        for line, fields in records(a.gtf):
            if fields is None or fields[2] != "transcript":
                continue
            gid = GENE_ID.search(fields[8]).group(1)
            transcripts[gid] += 1
            if READTHROUGH_TAG in fields[8]:
                readthrough[gid] += 1
                dropped_transcripts.add(TRANSCRIPT_ID.search(fields[8]).group(1))
    readthrough_genes = {g for g, n in readthrough.items() if n == transcripts[g]}

    kept = set()
    dropped = collections.Counter()
    by_type = collections.defaultdict(collections.Counter)
    written = collections.Counter()
    transcripts_dropped_in_kept_genes = 0
    untyped = 0
    with open(a.out_gtf, "w") as out:
        for line, fields in records(a.gtf):
            if fields is None:
                out.write(line)
                continue
            gid = GENE_ID.search(fields[8])
            if gid is None:
                sys.exit("record without gene_id: %r" % line[:80])
            gid = gid.group(1)
            gtype = GENE_TYPE.search(fields[8])
            if gtype is None:
                untyped += 1
            cls = gene_class(gtype.group(1)) if gtype else None
            if cls not in drop:
                cls = "readthrough" if gid in readthrough_genes else None
            if fields[2] == "gene":
                by_type[gtype.group(1) if gtype else "(no gene type)"]["dropped" if cls else "kept"] += 1
                if cls:
                    dropped[cls] += 1
                else:
                    kept.add(gid)
            if cls is not None:
                continue
            if dropped_transcripts and fields[2] != "gene":
                tid = TRANSCRIPT_ID.search(fields[8])
                if tid is not None and tid.group(1) in dropped_transcripts:
                    if fields[2] == "transcript":
                        transcripts_dropped_in_kept_genes += 1
                    continue
            out.write(line)
            written[fields[2]] += 1

    report = {
        "gtf": a.gtf, "gtf_sha256": sha256(a.gtf),
        "out_gtf": a.out_gtf, "out_gtf_sha256": sha256(a.out_gtf),
        "classes_dropped": sorted(drop),
        "genes_in_gtf": sum(v["kept"] + v["dropped"] for v in by_type.values()),
        "genes_kept": len(kept),
        "genes_dropped_by_class": {c: dropped[c] for c in CLASSES if c in drop},
        "readthrough_transcripts_dropped_in_kept_genes": transcripts_dropped_in_kept_genes,
        "records_written": dict(written),
        "records_without_gene_type": untyped,
        "by_gene_type": {t: dict(v) for t, v in sorted(by_type.items(), key=lambda kv: -sum(kv[1].values()))},
    }
    same = True
    if a.check_features:
        with open_text(a.check_features) as f:
            listed = set(line.rstrip("\n").split("\t")[0] for line in f if line.strip())
        same = listed == kept
        report["check_features"] = {
            "file": a.check_features, "genes_listed": len(listed), "identical": same,
            "kept_not_listed": len(kept - listed), "listed_not_kept": len(listed - kept),
            "examples_kept_not_listed": sorted(kept - listed)[:20],
            "examples_listed_not_kept": sorted(listed - kept)[:20],
        }
    if a.report:
        with open(a.report, "w") as f:
            json.dump(report, f, indent=1)
            f.write("\n")
    summary = {k: report[k] for k in ("genes_in_gtf", "genes_kept", "genes_dropped_by_class",
                                      "readthrough_transcripts_dropped_in_kept_genes")}
    if a.check_features:
        summary["identical_to_features_list"] = same
    print(json.dumps(summary))
    sys.exit(0 if same else 1)


if __name__ == "__main__":
    main()

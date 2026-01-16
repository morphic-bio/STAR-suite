#!/usr/bin/env python3
import argparse
import os
import sys


DEBUG_MARKER = "SLAM_DEBUG_GENES"


def find_default_root():
    env = os.environ.get("GEDI_ROOT")
    if env:
        return env
    candidate = os.path.join(os.getcwd(), "test", "tmp_slam_fixture", "gedi_install", "gedi")
    if os.path.isdir(candidate):
        return candidate
    return ""


def read_file(path):
    with open(path, "r", encoding="utf-8") as handle:
        return handle.read()


def write_file(path, text):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def insert_after(text, anchor, addition, label):
    if addition.strip() in text:
        return text
    idx = text.find(anchor)
    if idx == -1:
        raise RuntimeError(f"Anchor not found for {label}")
    return text[: idx + len(anchor)] + addition + text[idx + len(anchor) :]


def insert_before(text, anchor, addition, label):
    if addition.strip() in text:
        return text
    idx = text.find(anchor)
    if idx == -1:
        raise RuntimeError(f"Anchor not found for {label}")
    return text[:idx] + addition + text[idx:]


def apply_patch(text):
    if DEBUG_MARKER in text:
        return text, False

    text = insert_after(
        text,
        "import gedi.util.StringUtils;\n",
        "import gedi.util.io.text.LineOrientedFile;\nimport gedi.util.io.text.LineWriter;\n",
        "imports",
    )

    text = insert_after(
        text,
        "\tprivate int verbose = 0;\n",
        "\tprivate HashSet<String> debugGenes = null;\n"
        "\tprivate LineWriter debugWriter = null;\n"
        "\tprivate int debugMaxReads = 0;\n"
        "\tprivate AtomicInteger debugReadCount = new AtomicInteger(0);\n",
        "debug fields",
    )

    text = insert_after(
        text,
        "\t\tthis.highmem = highmem;\n",
        "\t\tinitDebug();\n",
        "constructor init",
    )

    debug_methods = (
        "\n"
        "\tprivate void initDebug() {\n"
        "\t\tdebugGenes = parseDebugGenes();\n"
        "\t\tif (debugGenes == null || debugGenes.isEmpty()) return;\n"
        "\t\tdebugMaxReads = parseDebugMaxReads();\n"
        "\t\tString out = readPropEnv(\"slam.debug.out\", \"SLAM_DEBUG_OUT\");\n"
        "\t\tif (out == null || out.isEmpty()) out = \"slam_debug.tsv\";\n"
        "\t\tString verboseStr = readPropEnv(\"slam.debug.verbose\", \"SLAM_DEBUG_VERBOSE\");\n"
        "\t\tif (verboseStr != null && verboseStr.length() > 0) {\n"
        "\t\t\ttry {\n"
        "\t\t\t\tsetVerbose(Integer.parseInt(verboseStr));\n"
        "\t\t\t} catch (NumberFormatException e) {\n"
        "\t\t\t\tthrow new RuntimeException(\"Invalid SLAM debug verbose value: \" + verboseStr, e);\n"
        "\t\t\t}\n"
        "\t\t}\n"
        "\t\ttry {\n"
        "\t\t\tdebugWriter = new LineOrientedFile(out).write();\n"
        "\t\t\tdebugWriter.write(\"Gene\\tRead\\tStrand\\tOppositeStrand\\tGeneConsistent\\tTranscriptCount\\tConsistentTranscripts\\tOverlapGeneCount\\tOverlapGenes\\tReadLen\\tReadLen1\\tReadLen2\\tDistinctIndex\\tWeight\\n\");\n"
        "\t\t\tfinal LineWriter writerRef = debugWriter;\n"
        "\t\t\tRuntime.getRuntime().addShutdownHook(new Thread(() -> {\n"
        "\t\t\t\ttry { writerRef.close(); } catch (Exception e) {}\n"
        "\t\t\t}));\n"
        "\t\t} catch (Exception e) {\n"
        "\t\t\tthrow new RuntimeException(\"Cannot open SLAM debug output: \" + out, e);\n"
        "\t\t}\n"
        "\t}\n"
        "\n"
        "\tprivate static String readPropEnv(String prop, String env) {\n"
        "\t\tString val = System.getProperty(prop);\n"
        "\t\tif (val == null || val.length() == 0) val = System.getenv(env);\n"
        "\t\treturn val;\n"
        "\t}\n"
        "\n"
        "\tprivate static HashSet<String> parseDebugGenes() {\n"
        "\t\tString raw = readPropEnv(\"slam.debug.genes\", \"SLAM_DEBUG_GENES\");\n"
        "\t\tif (raw == null) return null;\n"
        "\t\traw = raw.trim();\n"
        "\t\tif (raw.isEmpty()) return null;\n"
        "\t\tHashSet<String> set = new HashSet<>();\n"
        "\t\tfor (String token : raw.split(\"[,\\\\s]+\")) {\n"
        "\t\t\tif (!token.isEmpty()) set.add(token);\n"
        "\t\t}\n"
        "\t\treturn set.isEmpty() ? null : set;\n"
        "\t}\n"
        "\n"
        "\tprivate static int parseDebugMaxReads() {\n"
        "\t\tString raw = readPropEnv(\"slam.debug.max\", \"SLAM_DEBUG_MAX_READS\");\n"
        "\t\tif (raw == null || raw.isEmpty()) return 0;\n"
        "\t\ttry {\n"
        "\t\t\treturn Integer.parseInt(raw);\n"
        "\t\t} catch (NumberFormatException e) {\n"
        "\t\t\tthrow new RuntimeException(\"Invalid SLAM debug max reads: \" + raw, e);\n"
        "\t\t}\n"
        "\t}\n"
        "\n"
        "\tprivate boolean isDebugGene(ImmutableReferenceGenomicRegion<String> gene) {\n"
        "\t\treturn debugWriter != null && debugGenes != null && debugGenes.contains(gene.getData());\n"
        "\t}\n"
        "\n"
        "\tprivate void writeDebugLine(String geneId,\n"
        "\t\t\tImmutableReferenceGenomicRegion<AlignedReadsData> read,\n"
        "\t\t\tReadInfo readinfo,\n"
        "\t\t\tint transcriptCount,\n"
        "\t\t\tint consistentTranscripts,\n"
        "\t\t\tHashSet<String> overlapGenes,\n"
        "\t\t\tint distinctIndex,\n"
        "\t\t\tdouble weight,\n"
        "\t\t\tint readLen1,\n"
        "\t\t\tint readLen2) {\n"
        "\t\tif (debugWriter == null) return;\n"
        "\t\tint idx = debugReadCount.incrementAndGet();\n"
        "\t\tif (debugMaxReads > 0 && idx > debugMaxReads) return;\n"
        "\t\tString overlapList = overlapGenes == null ? \"\" : StringUtils.concat(\",\", overlapGenes);\n"
        "\t\tint overlapCount = overlapGenes == null ? 0 : overlapGenes.size();\n"
        "\t\tint readLen = read.getRegion().getTotalLength();\n"
        "\t\tsynchronized (debugWriter) {\n"
        "\t\t\ttry {\n"
        "\t\t\t\tdebugWriter.writef(\"%s\\t%s\\t%s\\t%d\\t%d\\t%d\\t%d\\t%d\\t%s\\t%d\\t%d\\t%d\\t%d\\t%.6f\\n\",\n"
        "\t\t\t\t\tgeneId,\n"
        "\t\t\t\t\tread.toLocationString(),\n"
        "\t\t\t\t\tread.getReference().getStrand().toString(),\n"
        "\t\t\t\t\treadinfo.isOppositeStrand() ? 1 : 0,\n"
        "\t\t\t\t\treadinfo.genes.length > 0 ? 1 : 0,\n"
        "\t\t\t\t\ttranscriptCount,\n"
        "\t\t\t\t\tconsistentTranscripts,\n"
        "\t\t\t\t\toverlapCount,\n"
        "\t\t\t\t\toverlapList,\n"
        "\t\t\t\t\treadLen,\n"
        "\t\t\t\t\treadLen1,\n"
        "\t\t\t\t\treadLen2,\n"
        "\t\t\t\t\tdistinctIndex,\n"
        "\t\t\t\t\tweight);\n"
        "\t\t\t} catch (Exception e) {\n"
        "\t\t\t\tthrow new RuntimeException(\"Cannot write SLAM debug output!\", e);\n"
        "\t\t\t}\n"
        "\t\t}\n"
        "\t}\n"
    )

    text = insert_before(
        text,
        "\tpublic void countNumis",
        debug_methods,
        "debug methods",
    )

    text = insert_after(
        text,
        "\t\tReadInfo readinfo = new ReadInfo(gene,tr);\n",
        "\t\tboolean debugGene = isDebugGene(gene);\n",
        "debug gene flag",
    )

    text = insert_after(
        text,
        "\t\t\treadinfo.overlapgenes = EI.wrap(tr).map(t->t.getData().getGeneId()).unique(false).toArray(String.class);\n",
        "\t\t\tHashSet<String> debugOverlapGenes = null;\n"
        "\t\t\tint debugConsistentTranscripts = -1;\n"
        "\t\t\tif (debugGene) {\n"
        "\t\t\t\tdebugConsistentTranscripts = EI.wrap(tr).filter(t->isConsistent(t, read)).count();\n"
        "\t\t\t\tdebugOverlapGenes = genomic.getTranscripts().ei(read).filter(t->isConsistent(t, read)).map(t->t.getData().getGeneId()).set();\n"
        "\t\t\t}\n",
        "debug per-read summary",
    )

    text = insert_after(
        text,
        "\t\t\t\t\treadLen2=this.readLength2.accumulateAndGet(trl2, Math::max);\n",
        "\t\t\t\t\tif (debugGene) {\n"
        "\t\t\t\t\t\tdouble debugWeight = rd.getTotalCountForDistinct(d, mode);\n"
        "\t\t\t\t\t\twriteDebugLine(gene.getData(), read, readinfo, tr.size(), debugConsistentTranscripts, debugOverlapGenes, d, debugWeight, trl1, trl2);\n"
        "\t\t\t\t\t}\n",
        "debug per-distinct log",
    )

    return text, True


def main():
    parser = argparse.ArgumentParser(
        description="Patch GEDI SlamCollector to emit gene-targeted debug TSV output."
    )
    parser.add_argument(
        "--gedi-root",
        default=find_default_root(),
        help="Path to gedi repo root (contains Gedi/). Default: $GEDI_ROOT or test/tmp_slam_fixture/gedi_install/gedi",
    )
    parser.add_argument(
        "--revert",
        action="store_true",
        help="Restore SlamCollector.java from backup (.bak).",
    )
    args = parser.parse_args()

    if not args.gedi_root:
        parser.error("Could not determine GEDI root; pass --gedi-root.")

    slam_collector = os.path.join(
        args.gedi_root, "Gedi", "src", "gedi", "slam", "SlamCollector.java"
    )
    if not os.path.isfile(slam_collector):
        raise FileNotFoundError(f"Missing SlamCollector.java at {slam_collector}")

    backup_path = slam_collector + ".bak"

    if args.revert:
        if not os.path.isfile(backup_path):
            raise FileNotFoundError(f"Backup not found: {backup_path}")
        write_file(slam_collector, read_file(backup_path))
        print(f"Restored {slam_collector} from {backup_path}")
        return 0

    text = read_file(slam_collector)
    new_text, changed = apply_patch(text)
    if not changed:
        print("Patch already applied; no changes made.")
        return 0

    if not os.path.isfile(backup_path):
        write_file(backup_path, text)

    write_file(slam_collector, new_text)
    print(f"Patched {slam_collector}")
    print("Next: rebuild GEDI (e.g., mvn -f <gedi-root>/Gedi package).")
    print("Set debug env vars: SLAM_DEBUG_GENES, SLAM_DEBUG_OUT, SLAM_DEBUG_MAX_READS, SLAM_DEBUG_VERBOSE.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

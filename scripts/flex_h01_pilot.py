#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import os
import struct
import sys
from collections import defaultdict

# Must match core/legacy/source/FlexHashScreen.cpp
CACHE_MAGIC = b"FH01SEQ1"
CACHE_VERSION_SAMPLE_AWARE = 2
CACHE_KMER_LENGTH = 50
CACHE_RECORD_SIZE = 24
NEG_PROBE_AMBIG = 1
# Hash-screen cache class: H2 KEEP (2-substitution) — FlexHashScreen treats like H1 for non-exact path
CACHE_CLASS_H2_KEEP = 3


def read_probe_list_index(path):
    gene_to_idx = {}
    idx_to_gene = {}
    with open(path, "r", encoding="utf-8") as handle:
        for idx, line in enumerate(handle, start=1):
            gene_id = line.strip()
            if not gene_id:
                continue
            if gene_id not in gene_to_idx:
                gene_to_idx[gene_id] = idx
                idx_to_gene[idx] = gene_id
    return gene_to_idx, idx_to_gene


def normalize_qname(qname):
    return qname[1:] if qname.startswith("@") else qname


def read_probe_fasta(path):
    probes = []
    header = None
    seq_lines = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    probes.append(parse_probe_record(header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        probes.append(parse_probe_record(header, "".join(seq_lines)))
    return probes


def parse_probe_record(header, seq):
    parts = header.split("|")
    gene_id = parts[0]
    gene_name = parts[1] if len(parts) > 1 else gene_id
    probe_id = parts[2] if len(parts) > 2 else header
    return {
        "gene_id": gene_id,
        "gene_name": gene_name,
        "probe_id": probe_id,
        "seq": seq,
    }


def load_feature_rows(features_tsv):
    rows = []
    with open(features_tsv, "r", encoding="utf-8") as handle:
        for idx, line in enumerate(handle, start=1):
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            gene_id = parts[0]
            gene_name = parts[1] if len(parts) > 1 and parts[1] else gene_id
            rows.append((idx, gene_id, gene_name))
    return rows


def load_matrix_row_totals(matrix_mtx):
    totals = defaultdict(int)
    saw_dimensions = False
    with open(matrix_mtx, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) == 3 and not saw_dimensions:
                # Matrix dimensions header.
                saw_dimensions = True
                continue
            if len(parts) != 3:
                continue
            row = int(parts[0])
            value = int(parts[2])
            totals[row] += value
    return totals


def cmd_select_probes(args):
    gene_to_idx, _ = read_probe_list_index(args.probe_list)
    probes = read_probe_fasta(args.probes_fasta)
    feature_rows = load_feature_rows(args.features_tsv)
    row_totals = load_matrix_row_totals(args.matrix_mtx)

    selected = []
    seen_genes = set()
    probe_by_gene = {}
    for probe in probes:
        probe_by_gene.setdefault(probe["gene_id"], []).append(probe)

    for row_idx, gene_id, gene_name in feature_rows:
        total = row_totals.get(row_idx, 0)
        if total <= 0 or gene_id in seen_genes:
            continue
        gene_probes = probe_by_gene.get(gene_id, [])
        if not gene_probes:
            continue
        probe = gene_probes[0]
        selected.append({
            "gene_id": gene_id,
            "gene_name": gene_name,
            "gene_idx15": gene_to_idx.get(gene_id, 0),
            "probe_id": probe["probe_id"],
            "probe_seq": probe["seq"],
            "row_idx": row_idx,
            "matrix_total": total,
        })
        seen_genes.add(gene_id)
        if len(selected) >= args.count:
            break

    if len(selected) < args.count:
        raise SystemExit(
            f"requested {args.count} probes but only found {len(selected)} expressed genes with probes"
        )

    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "gene_id",
                "gene_name",
                "gene_idx15",
                "probe_id",
                "probe_seq",
                "row_idx",
                "matrix_total",
            ]
        )
        for row in selected:
            writer.writerow(
                [
                    row["gene_id"],
                    row["gene_name"],
                    row["gene_idx15"],
                    row["probe_id"],
                    row["probe_seq"],
                    row["row_idx"],
                    row["matrix_total"],
                ]
            )


def first_whitelist_barcode(path):
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line:
                return line.split("\t")[0]
    raise SystemExit(f"no barcode found in {path}")


def read_whitelist_barcodes(path, count):
    out = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            out.append(line.split("\t")[0])
            if len(out) >= count:
                return out
    raise SystemExit(f"requested {count} whitelist barcodes but only found {len(out)} in {path}")


def first_sample_tag(path):
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2 and parts[1]:
                return parts[1], parts[0]
    raise SystemExit(f"no sample tag found in {path}")


def variant_records(seq):
    bases = "ACGT"
    yield (0, -1, "=", seq)
    seq_chars = list(seq)
    for pos, ref in enumerate(seq_chars):
        for alt in bases:
            if alt == ref:
                continue
            var = seq_chars[:]
            var[pos] = alt
            yield (1, pos, alt, "".join(var))


def encode_window_50(seq):
    """Pack 50 bp into (seq_lo, seq_hi) matching FlexHashScreenCache::encodeWindow."""
    if len(seq) != CACHE_KMER_LENGTH:
        raise ValueError(f"expected {CACHE_KMER_LENGTH} bp, got {len(seq)}")
    mp = {"A": 0, "C": 1, "G": 2, "T": 3}
    lo = 0
    hi = 0
    for c in seq.upper():
        if c not in mp:
            raise ValueError(f"invalid base {c!r}")
        code = mp[c]
        hi = ((hi << 2) | (lo >> 62)) & 0xFFFFFFFFFFFFFFFF
        lo = ((lo << 2) | code) & 0xFFFFFFFFFFFFFFFF
    return lo, hi


def decode_window_50(seq_lo, seq_hi):
    """Unpack (seq_lo, seq_hi) to 50 bp string."""
    inv = "ACGT"
    combined = (seq_hi << 64) | seq_lo
    out = []
    for i in range(CACHE_KMER_LENGTH):
        shift = 98 - 2 * i
        code = (combined >> shift) & 3
        out.append(inv[code])
    return "".join(out)


def load_binary_hash_cache(path):
    """Read FH01SEQ1 cache; return list of dicts with seq_lo, seq_hi, sequence, gene_idx15, cache_class, negative_code, sample_idx."""
    with open(path, "rb") as handle:
        header = handle.read(24)
        if len(header) != 24:
            raise SystemExit(f"truncated cache header: {path}")
        magic, version, kmer_len, rec_size, nrec = struct.unpack("<8sHHIQ", header)
        if magic != CACHE_MAGIC:
            raise SystemExit(f"bad magic in {path}")
        if version not in (1, CACHE_VERSION_SAMPLE_AWARE) or kmer_len != CACHE_KMER_LENGTH or rec_size != CACHE_RECORD_SIZE:
            raise SystemExit(f"unsupported cache format version={version} kmer={kmer_len} rec={rec_size}")
        out = []
        for _ in range(int(nrec)):
            raw = handle.read(24)
            if len(raw) != 24:
                raise SystemExit(f"truncated cache body: {path}")
            seq_lo, seq_hi, gene15, cclass, negcode, res = struct.unpack("<QQIBBH", raw)
            sample_idx = res if version >= CACHE_VERSION_SAMPLE_AWARE else 0
            seq = decode_window_50(seq_lo, seq_hi)
            out.append(
                {
                    "seq_lo": seq_lo,
                    "seq_hi": seq_hi,
                    "sequence": seq,
                    "resolved_gene_idx15": gene15,
                    "cache_class": cclass,
                    "negative_code": negcode,
                    "sample_idx": sample_idx,
                }
            )
        return out


def write_binary_hash_cache(path, records):
    """Write FH01SEQ1 v2. records: iterable of (seq_lo, seq_hi, gene_idx15, cache_class, negative_code, sample_idx)."""
    recs = list(records)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "wb") as handle:
        handle.write(
            struct.pack(
                "<8sHHIQ",
                CACHE_MAGIC,
                CACHE_VERSION_SAMPLE_AWARE,
                CACHE_KMER_LENGTH,
                CACHE_RECORD_SIZE,
                len(recs),
            )
        )
        for seq_lo, seq_hi, gene15, cclass, negcode, sample_idx in recs:
            handle.write(struct.pack("<QQIBBH", seq_lo, seq_hi, gene15, cclass, negcode, sample_idx & 0xFFFF))


def variant_records_h2(seq):
    """All 2-substitution variants (Hamming distance 2). Yields (i, j, alt_i, alt_j, variant_seq)."""
    bases = "ACGT"
    chars = list(seq.upper())
    n = len(chars)
    if n != CACHE_KMER_LENGTH:
        raise ValueError(f"H2 variants need length {CACHE_KMER_LENGTH}")
    for i in range(n):
        ref_i = chars[i]
        for j in range(i + 1, n):
            ref_j = chars[j]
            for alt_i in bases:
                if alt_i == ref_i:
                    continue
                for alt_j in bases:
                    if alt_j == ref_j:
                        continue
                    var = chars[:]
                    var[i] = alt_i
                    var[j] = alt_j
                    yield (i, j, alt_i, alt_j, "".join(var))


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8")
    return open(path, "w", encoding="utf-8")


def read_selected_probes(path):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row["gene_idx15"] = int(row["gene_idx15"])
            row["row_idx"] = int(row["row_idx"])
            row["matrix_total"] = int(row["matrix_total"])
            rows.append(row)
    return rows


def load_selected_gene_ids(path):
    return {row["gene_id"] for row in read_selected_probes(path)}


def cmd_make_synth_fastq(args):
    probes = read_selected_probes(args.selected_probes)
    n_records = len(probes) * 151
    barcodes = read_whitelist_barcodes(args.cb_whitelist, n_records)
    umi = args.synthetic_umi
    if len(umi) != 12:
        raise SystemExit(f"expected 12 bp synthetic UMI, got {len(umi)}")
    sample_tag = args.sample_tag.upper() if args.sample_tag else ""
    sample_label = ""
    if args.sample_whitelist:
        file_tag, file_label = first_sample_tag(args.sample_whitelist)
        sample_label = file_label
        if not sample_tag:
            sample_tag = file_tag.upper()
    use_real_flex_layout = bool(sample_tag)
    if use_real_flex_layout and len(sample_tag) != 8:
        raise SystemExit(f"expected 8 bp sample tag, got {len(sample_tag)} for {sample_tag}")
    if use_real_flex_layout:
        if args.probe_offset < 0 or args.sample_offset < 0 or args.read_length <= 0:
            raise SystemExit("probe_offset, sample_offset, and read_length must be positive")
        if args.probe_offset + 50 > args.read_length:
            raise SystemExit("probe placement exceeds read length")
        if args.sample_offset + 8 > args.read_length:
            raise SystemExit("sample tag placement exceeds read length")
        fill_base = args.fill_base.upper()
        if fill_base not in {"A", "C", "G", "T"}:
            raise SystemExit(f"fill base must be A/C/G/T, got {fill_base}")
    else:
        left_flank = args.left_flank.upper()
        right_flank = args.right_flank.upper()

    os.makedirs(os.path.dirname(args.r1_fastq), exist_ok=True)
    os.makedirs(os.path.dirname(args.r2_fastq), exist_ok=True)
    manifest_rows = []
    record_id = 1
    with open_maybe_gzip(args.r1_fastq) as r1_handle, open_maybe_gzip(args.r2_fastq) as r2_handle:
        for probe in probes:
            for distance, pos, alt, seq in variant_records(probe["probe_seq"]):
                qname = f"H01|{record_id:06d}"
                note = "ref" if distance == 0 else f"p{pos:02d}{alt}"
                barcode = barcodes[record_id - 1]
                if len(barcode) != 16:
                    raise SystemExit(f"expected 16 bp barcode, got {len(barcode)} for {barcode}")
                r1_seq = barcode + umi
                if use_real_flex_layout:
                    read_chars = [fill_base] * args.read_length
                    for idx, base in enumerate(seq):
                        read_chars[args.probe_offset + idx] = base
                    for idx, base in enumerate(sample_tag):
                        read_chars[args.sample_offset + idx] = base
                    read_seq = "".join(read_chars)
                else:
                    read_seq = left_flank + seq + right_flank
                r1_qual = "I" * len(r1_seq)
                r2_qual = "I" * len(read_seq)
                r1_handle.write(f"@{qname}\n{r1_seq}\n+\n{r1_qual}\n")
                r2_handle.write(f"@{qname}\n{read_seq}\n+\n{r2_qual}\n")
                manifest_rows.append(
                    {
                        "qname": qname,
                        "core_sequence": seq,
                        "read_sequence": read_seq,
                        "cell_barcode": barcode,
                        "sample_tag": sample_tag,
                        "sample_label": sample_label,
                        "distance": distance,
                        "mut_pos": pos,
                        "mut_alt": alt,
                        "variant_note": note,
                        "gene_id": probe["gene_id"],
                        "gene_name": probe["gene_name"],
                        "gene_idx15": probe["gene_idx15"],
                        "probe_id": probe["probe_id"],
                    }
                )
                record_id += 1

    with open(args.manifest_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "qname",
                "core_sequence",
                "read_sequence",
                "cell_barcode",
                "sample_tag",
                "sample_label",
                "distance",
                "mut_pos",
                "mut_alt",
                "variant_note",
                "gene_id",
                "gene_name",
                "gene_idx15",
                "probe_id",
            ],
        )
        writer.writeheader()
        writer.writerows(manifest_rows)


def parse_reject_log(path):
    events = defaultdict(list)
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            qname = normalize_qname(row.get("qname", ""))
            if not qname:
                continue
            events[qname].append(row)
    return events


def classify_qname_event(event_rows):
    keep_gene = 0
    keep_reason = ""
    for row in event_rows:
        reason = row.get("reason", "")
        gene_idx15 = int(row.get("geneIdx15", "0") or 0)
        if reason == "KEEP_HASH" and gene_idx15 > 0:
            keep_gene = gene_idx15
            keep_reason = reason
    if keep_gene > 0:
        return ("KEEP", keep_gene, keep_reason)
    return ("DENY", 0, "")


def cmd_build_cache(args):
    manifest = {}
    with open(args.manifest_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row["distance"] = int(row["distance"])
            row["gene_idx15"] = int(row["gene_idx15"])
            row["mut_pos"] = int(row["mut_pos"])
            manifest[row["qname"]] = row

    reject_events = parse_reject_log(args.reject_log)

    by_sequence = defaultdict(list)
    qname_rows = []
    for qname, meta in manifest.items():
        state, keep_gene_idx15, keep_reason = classify_qname_event(reject_events.get(qname, []))
        qname_row = {
            "qname": qname,
            "sequence": meta["core_sequence"],
            "input_distance": meta["distance"],
            "input_gene_id": meta["gene_id"],
            "input_gene_idx15": meta["gene_idx15"],
            "probe_id": meta["probe_id"],
            "state": state,
            "resolved_gene_idx15": keep_gene_idx15,
            "resolved_gene_id": meta["gene_id"] if keep_gene_idx15 == meta["gene_idx15"] else "",
            "reason": keep_reason if keep_reason else "NON_KEEP",
        }
        qname_rows.append(qname_row)
        by_sequence[meta["core_sequence"]].append(qname_row)

    with open(args.qname_cache_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "qname",
                "sequence",
                "input_distance",
                "input_gene_id",
                "input_gene_idx15",
                "probe_id",
                "state",
                "resolved_gene_idx15",
                "resolved_gene_id",
                "reason",
            ],
        )
        writer.writeheader()
        writer.writerows(qname_rows)

    sequence_rows = []
    for sequence, rows in by_sequence.items():
        keep_genes = {row["resolved_gene_idx15"] for row in rows if row["state"] == "KEEP"}
        has_deny = any(row["state"] != "KEEP" for row in rows)
        if has_deny or len(keep_genes) != 1:
            cache_class = 2
            resolved_gene_idx15 = 0
        else:
            cache_class = min(row["input_distance"] for row in rows if row["state"] == "KEEP")
            resolved_gene_idx15 = next(iter(keep_genes))
        source_gene_ids = sorted({row["input_gene_id"] for row in rows})
        sequence_rows.append(
            {
                "sequence": sequence,
                "cache_class": cache_class,
                "resolved_gene_idx15": resolved_gene_idx15,
                "source_gene_ids": ",".join(source_gene_ids),
                "n_records": len(rows),
                "n_keep": sum(1 for row in rows if row["state"] == "KEEP"),
                "n_deny": sum(1 for row in rows if row["state"] != "KEEP"),
            }
        )

    sequence_rows.sort(key=lambda row: (row["cache_class"], row["resolved_gene_idx15"], row["sequence"]))
    with open(args.sequence_cache_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "sequence",
                "cache_class",
                "resolved_gene_idx15",
                "source_gene_ids",
                "n_records",
                "n_keep",
                "n_deny",
            ],
        )
        writer.writeheader()
        writer.writerows(sequence_rows)


def parse_fastq_handle(handle):
    while True:
        h = handle.readline()
        if not h:
            return
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not qual:
            raise SystemExit("truncated FASTQ record")
        yield h.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def open_fastq(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def open_fastq_writer(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if path.endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8")
    return open(path, "w", encoding="utf-8")


def load_sequence_cache(path):
    cache = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cache[row["sequence"]] = (
                int(row["cache_class"]),
                int(row["resolved_gene_idx15"]),
            )
    return cache


def load_barcode_order(path):
    out = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            barcode = raw.strip()
            if barcode:
                out.append(barcode)
    return out


def load_features_order(path):
    out = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            parts = raw.rstrip("\n").split("\t")
            if parts and parts[0]:
                out.append(parts[0])
    return out


def load_mtx_entries(path):
    entries = defaultdict(list)
    with open(path, "r", encoding="utf-8") as handle:
        saw_dimensions = False
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) == 3 and not saw_dimensions:
                saw_dimensions = True
                continue
            if len(parts) != 3:
                continue
            row_idx = int(parts[0])
            col_idx = int(parts[1])
            value = int(parts[2])
            entries[col_idx].append((row_idx, value))
    return entries


def mex_barcode_keys(row):
    keys = []
    cell_barcode = row.get("cell_barcode", "")
    sample_tag = row.get("sample_tag", "")
    if cell_barcode and sample_tag:
        keys.append(cell_barcode + sample_tag)
    if cell_barcode:
        keys.append(cell_barcode)
    return keys


def cmd_h2_make_synth_fastq(args):
    """Generate R1/R2 FASTQ for all H2 variants of unique H0 sequences from a binary hash cache."""
    _, idx_to_gene = read_probe_list_index(args.probe_list)
    cache_rows = load_binary_hash_cache(args.binary_cache)
    h0_by_seq = {}
    for row in cache_rows:
        if row["cache_class"] != 0:
            continue
        key = (row["seq_lo"], row["seq_hi"])
        if key not in h0_by_seq:
            gid = idx_to_gene.get(row["resolved_gene_idx15"], "")
            h0_by_seq[key] = {
                "sequence": row["sequence"],
                "gene_idx15": row["resolved_gene_idx15"],
                "gene_id": gid,
                "probe_id": f"gene{row['resolved_gene_idx15']}",
            }
    parents = list(h0_by_seq.values())
    if args.parent_limit is not None:
        parents = parents[: args.parent_limit]
    num_shards = max(1, args.num_shards)
    shard_id = args.shard_id
    if shard_id < 0 or shard_id >= num_shards:
        raise SystemExit(f"shard_id must be in [0, {num_shards})")

    selected = []
    for p in parents:
        h = int(hashlib.md5(p["sequence"].encode()).hexdigest(), 16) % num_shards
        if h == shard_id:
            selected.append(p)

    variant_tasks = []
    for probe in selected:
        for i, j, alt_i, alt_j, var_seq in variant_records_h2(probe["sequence"]):
            variant_tasks.append((probe, i, j, alt_i, alt_j, var_seq))

    n_records = len(variant_tasks)
    if n_records == 0:
        raise SystemExit("no H2 variants for this shard/parent set")
    max_cb = args.max_reads_per_shard
    if n_records > max_cb:
        raise SystemExit(
            f"shard has {n_records} variants but max_reads_per_shard={max_cb} (whitelist size). "
            f"Increase --num-shards (>= {((len(parents) * 11025 + max_cb - 1) // max_cb)}) or --parent-limit."
        )

    barcodes = read_whitelist_barcodes(args.cb_whitelist, n_records)
    umi = args.synthetic_umi
    if len(umi) != 12:
        raise SystemExit(f"expected 12 bp synthetic UMI, got {len(umi)}")
    sample_tag = args.sample_tag.upper() if args.sample_tag else ""
    sample_label = ""
    if args.sample_whitelist:
        file_tag, file_label = first_sample_tag(args.sample_whitelist)
        sample_label = file_label
        if not sample_tag:
            sample_tag = file_tag.upper()
    use_real_flex_layout = bool(sample_tag)
    if use_real_flex_layout and len(sample_tag) != 8:
        raise SystemExit(f"expected 8 bp sample tag, got {len(sample_tag)}")
    if use_real_flex_layout:
        if args.probe_offset + CACHE_KMER_LENGTH > args.read_length:
            raise SystemExit("probe placement exceeds read length")
        if args.sample_offset + 8 > args.read_length:
            raise SystemExit("sample tag placement exceeds read length")
        fill_base = args.fill_base.upper()
        if fill_base not in {"A", "C", "G", "T"}:
            raise SystemExit(f"fill base must be A/C/G/T, got {fill_base}")

    os.makedirs(os.path.dirname(args.r1_fastq), exist_ok=True)
    os.makedirs(os.path.dirname(args.r2_fastq), exist_ok=True)
    manifest_rows = []
    record_id = 1
    with open_maybe_gzip(args.r1_fastq) as r1_handle, open_maybe_gzip(args.r2_fastq) as r2_handle:
        for probe, i, j, alt_i, alt_j, var_seq in variant_tasks:
            qname = f"H02|{record_id:09d}"
            note = f"p{i:02d}{alt_i}_p{j:02d}{alt_j}"
            barcode = barcodes[record_id - 1]
            r1_seq = barcode + umi
            if use_real_flex_layout:
                read_chars = [fill_base] * args.read_length
                for idx, base in enumerate(var_seq):
                    read_chars[args.probe_offset + idx] = base
                for idx, base in enumerate(sample_tag):
                    read_chars[args.sample_offset + idx] = base
                read_seq = "".join(read_chars)
            else:
                read_seq = args.left_flank + var_seq + args.right_flank
            r1_qual = "I" * len(r1_seq)
            r2_qual = "I" * len(read_seq)
            r1_handle.write(f"@{qname}\n{r1_seq}\n+\n{r1_qual}\n")
            r2_handle.write(f"@{qname}\n{read_seq}\n+\n{r2_qual}\n")
            manifest_rows.append(
                {
                    "qname": qname,
                    "core_sequence": var_seq,
                    "read_sequence": read_seq,
                    "cell_barcode": barcode,
                    "sample_tag": sample_tag,
                    "sample_label": sample_label,
                    "distance": 2,
                    "mut_pos": i,
                    "mut_alt": alt_i,
                    "mut_pos2": j,
                    "mut_alt2": alt_j,
                    "variant_note": note,
                    "gene_id": probe["gene_id"],
                    "gene_name": probe["gene_id"],
                    "gene_idx15": probe["gene_idx15"],
                    "probe_id": probe["probe_id"],
                }
            )
            record_id += 1

    with open(args.manifest_tsv, "w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "qname",
            "core_sequence",
            "read_sequence",
            "cell_barcode",
            "sample_tag",
            "sample_label",
            "distance",
            "mut_pos",
            "mut_alt",
            "mut_pos2",
            "mut_alt2",
            "variant_note",
            "gene_id",
            "gene_name",
            "gene_idx15",
            "probe_id",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(
        f"h2-make-synth-fastq: shard {shard_id}/{num_shards} wrote {n_records} reads "
        f"from {len(selected)} parent probes",
        file=sys.stderr,
    )


def cmd_h2_build_cache_from_mex(args):
    """Same as build-cache-from-mex, but KEEP rows use cache_class=3 (H2) instead of min(distance)."""
    gene_to_idx, _ = read_probe_list_index(args.probe_list)
    manifest_rows = []
    by_qname = {}
    with open(args.manifest_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row["distance"] = int(row["distance"])
            row["gene_idx15"] = int(row["gene_idx15"])
            row["mut_pos"] = int(row["mut_pos"])
            manifest_rows.append(row)
            by_qname[row["qname"]] = row

    barcodes = load_barcode_order(args.barcodes_tsv)
    features = load_features_order(args.features_tsv)
    col_entries = load_mtx_entries(args.matrix_mtx)
    barcode_to_col = {barcode: idx for idx, barcode in enumerate(barcodes, start=1)}

    qname_rows = []
    by_sequence = defaultdict(list)
    for row in manifest_rows:
        entries = []
        for barcode_key in mex_barcode_keys(row):
            col_idx = barcode_to_col.get(barcode_key)
            if col_idx is not None:
                entries = col_entries.get(col_idx, [])
                break
        nonzero = [(features[r - 1], value) for r, value in entries if value > 0]
        if len(nonzero) == 1 and nonzero[0][1] == 1:
            resolved_gene_id = nonzero[0][0]
            resolved_gene_idx15 = gene_to_idx.get(resolved_gene_id, 0)
            state = (
                "KEEP"
                if resolved_gene_idx15 > 0 and resolved_gene_id == row["gene_id"]
                else "DENY"
            )
        else:
            resolved_gene_id = ""
            resolved_gene_idx15 = 0
            state = "DENY"
        qname_row = {
            "qname": row["qname"],
            "sequence": row["core_sequence"],
            "input_distance": row["distance"],
            "input_gene_id": row["gene_id"],
            "input_gene_idx15": row["gene_idx15"],
            "probe_id": row["probe_id"],
            "state": state,
            "resolved_gene_idx15": resolved_gene_idx15,
            "resolved_gene_id": resolved_gene_id,
            "reason": "SYNTH_MEX_H2",
        }
        qname_rows.append(qname_row)
        by_sequence[row["core_sequence"]].append(qname_row)

    with open(args.qname_cache_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "qname",
                "sequence",
                "input_distance",
                "input_gene_id",
                "input_gene_idx15",
                "probe_id",
                "state",
                "resolved_gene_idx15",
                "resolved_gene_id",
                "reason",
            ],
        )
        writer.writeheader()
        writer.writerows(qname_rows)

    sequence_rows = []
    for sequence, rows in by_sequence.items():
        keep_genes = {row["resolved_gene_idx15"] for row in rows if row["state"] == "KEEP"}
        has_deny = any(row["state"] != "KEEP" for row in rows)
        if has_deny or len(keep_genes) != 1:
            cache_class = 2
            resolved_gene_idx15 = 0
            neg_code = NEG_PROBE_AMBIG
        else:
            cache_class = CACHE_CLASS_H2_KEEP
            resolved_gene_idx15 = next(iter(keep_genes))
            neg_code = 0
        source_gene_ids = sorted({row["input_gene_id"] for row in rows})
        sequence_rows.append(
            {
                "sequence": sequence,
                "cache_class": cache_class,
                "negative_code": neg_code,
                "resolved_gene_idx15": resolved_gene_idx15,
                "source_gene_ids": ",".join(source_gene_ids),
                "n_records": len(rows),
                "n_keep": sum(1 for row in rows if row["state"] == "KEEP"),
                "n_deny": sum(1 for row in rows if row["state"] != "KEEP"),
            }
        )

    sequence_rows.sort(key=lambda row: (row["cache_class"], row["resolved_gene_idx15"], row["sequence"]))
    with open(args.sequence_cache_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "sequence",
                "cache_class",
                "negative_code",
                "resolved_gene_idx15",
                "source_gene_ids",
                "n_records",
                "n_keep",
                "n_deny",
            ],
        )
        writer.writeheader()
        writer.writerows(sequence_rows)


def cmd_h2_write_binary_cache(args):
    """Write FH01SEQ1 from sequence_cache TSV (optionally filter to KEEP rows only)."""
    rows = []
    with open(args.sequence_cache_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cclass = int(row["cache_class"])
            if args.min_class is not None and cclass < args.min_class:
                continue
            if args.max_class is not None and cclass > args.max_class:
                continue
            if args.keep_only and cclass != CACHE_CLASS_H2_KEEP:
                continue
            rows.append(row)

    recs = []
    sample_idx = args.sample_idx
    neg_default = 0
    for row in rows:
        seq = row["sequence"]
        if len(seq) != CACHE_KMER_LENGTH:
            raise SystemExit(f"bad sequence length {len(seq)}")
        lo, hi = encode_window_50(seq)
        gene15 = int(row["resolved_gene_idx15"])
        cclass = int(row["cache_class"])
        negcode = int(row.get("negative_code", neg_default))
        if cclass == 2 and gene15 == 0:
            negcode = NEG_PROBE_AMBIG
        recs.append((lo, hi, gene15, cclass, negcode, sample_idx))

    recs.sort(key=lambda r: (r[1], r[0], r[5]))
    write_binary_hash_cache(args.output_bin, recs)
    print(f"wrote {len(recs)} records to {args.output_bin}", file=sys.stderr)


def cmd_build_cache_from_mex(args):
    gene_to_idx, _ = read_probe_list_index(args.probe_list)
    manifest_rows = []
    by_qname = {}
    with open(args.manifest_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row["distance"] = int(row["distance"])
            row["gene_idx15"] = int(row["gene_idx15"])
            row["mut_pos"] = int(row["mut_pos"])
            manifest_rows.append(row)
            by_qname[row["qname"]] = row

    barcodes = load_barcode_order(args.barcodes_tsv)
    features = load_features_order(args.features_tsv)
    col_entries = load_mtx_entries(args.matrix_mtx)
    barcode_to_col = {barcode: idx for idx, barcode in enumerate(barcodes, start=1)}

    qname_rows = []
    by_sequence = defaultdict(list)
    for row in manifest_rows:
        entries = []
        for barcode_key in mex_barcode_keys(row):
            col_idx = barcode_to_col.get(barcode_key)
            if col_idx is not None:
                entries = col_entries.get(col_idx, [])
                break
        nonzero = [(features[r - 1], value) for r, value in entries if value > 0]
        if len(nonzero) == 1 and nonzero[0][1] == 1:
            resolved_gene_id = nonzero[0][0]
            resolved_gene_idx15 = gene_to_idx.get(resolved_gene_id, 0)
            state = (
                "KEEP"
                if resolved_gene_idx15 > 0 and resolved_gene_id == row["gene_id"]
                else "DENY"
            )
        else:
            resolved_gene_id = ""
            resolved_gene_idx15 = 0
            state = "DENY"
        qname_row = {
            "qname": row["qname"],
            "sequence": row["core_sequence"],
            "input_distance": row["distance"],
            "input_gene_id": row["gene_id"],
            "input_gene_idx15": row["gene_idx15"],
            "probe_id": row["probe_id"],
            "state": state,
            "resolved_gene_idx15": resolved_gene_idx15,
            "resolved_gene_id": resolved_gene_id,
            "reason": "SYNTH_MEX",
        }
        qname_rows.append(qname_row)
        by_sequence[row["core_sequence"]].append(qname_row)

    with open(args.qname_cache_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "qname",
                "sequence",
                "input_distance",
                "input_gene_id",
                "input_gene_idx15",
                "probe_id",
                "state",
                "resolved_gene_idx15",
                "resolved_gene_id",
                "reason",
            ],
        )
        writer.writeheader()
        writer.writerows(qname_rows)

    sequence_rows = []
    for sequence, rows in by_sequence.items():
        keep_genes = {row["resolved_gene_idx15"] for row in rows if row["state"] == "KEEP"}
        has_deny = any(row["state"] != "KEEP" for row in rows)
        if has_deny or len(keep_genes) != 1:
            cache_class = 2
            resolved_gene_idx15 = 0
        else:
            cache_class = min(row["input_distance"] for row in rows if row["state"] == "KEEP")
            resolved_gene_idx15 = next(iter(keep_genes))
        source_gene_ids = sorted({row["input_gene_id"] for row in rows})
        sequence_rows.append(
            {
                "sequence": sequence,
                "cache_class": cache_class,
                "resolved_gene_idx15": resolved_gene_idx15,
                "source_gene_ids": ",".join(source_gene_ids),
                "n_records": len(rows),
                "n_keep": sum(1 for row in rows if row["state"] == "KEEP"),
                "n_deny": sum(1 for row in rows if row["state"] != "KEEP"),
            }
        )

    sequence_rows.sort(key=lambda row: (row["cache_class"], row["resolved_gene_idx15"], row["sequence"]))
    with open(args.sequence_cache_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "sequence",
                "cache_class",
                "resolved_gene_idx15",
                "source_gene_ids",
                "n_records",
                "n_keep",
                "n_deny",
            ],
        )
        writer.writeheader()
        writer.writerows(sequence_rows)


def decide_read(cache_hits):
    if not cache_hits:
        return ("PASS", 0, "", "")
    if any(hit["cache_class"] == 2 for hit in cache_hits):
        offsets = ",".join(str(hit["offset"]) for hit in cache_hits if hit["cache_class"] == 2)
        return ("DENY", 0, "2", offsets)
    gene_ids = {hit["resolved_gene_idx15"] for hit in cache_hits}
    if len(gene_ids) == 1:
        classes = [hit["cache_class"] for hit in cache_hits]
        offsets = ",".join(str(hit["offset"]) for hit in cache_hits)
        return ("KEEP", next(iter(gene_ids)), str(min(classes)), offsets)
    return ("PASS", 0, "conflict", ",".join(str(hit["offset"]) for hit in cache_hits))


def cmd_scan_fastq(args):
    cache = load_sequence_cache(args.sequence_cache_tsv)
    offsets = [int(part) for part in args.offsets.split(",") if part.strip()]
    summary = defaultdict(int)

    fastq_paths = []
    for item in args.fastq_r2:
        fastq_paths.extend(part for part in item.split(",") if part)

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)
    with open(args.output_tsv, "w", encoding="utf-8", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["qname", "decision", "resolved_gene_idx15", "detail", "offsets"])
        for fastq_path in fastq_paths:
            with open_fastq(fastq_path) as handle:
                for header, seq, _, _ in parse_fastq_handle(handle):
                    qname = normalize_qname(header)
                    hits = []
                    for offset in offsets:
                        if len(seq) < offset + 50:
                            continue
                        window = seq[offset:offset + 50]
                        cached = cache.get(window)
                        if cached is None:
                            continue
                        cache_class, gene_idx15 = cached
                        hits.append(
                            {
                                "offset": offset,
                                "cache_class": cache_class,
                                "resolved_gene_idx15": gene_idx15,
                            }
                        )
                    decision, gene_idx15, detail, hit_offsets = decide_read(hits)
                    writer.writerow([qname, decision, gene_idx15, detail, hit_offsets])
                    summary["total_reads"] += 1
                    summary[f"decision_{decision.lower()}"] += 1
                    if decision == "KEEP":
                        summary["keep_class_0" if detail == "0" else "keep_class_1"] += 1

    with open(args.summary_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for key in sorted(summary):
            writer.writerow([key, summary[key]])


def load_scan_decisions(path):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row["resolved_gene_idx15"] = int(row["resolved_gene_idx15"])
            rows.append(row)
    return rows


def parse_csv_path_list(value):
    return [part for part in value.split(",") if part]


def read_cb_whitelist(path):
    whitelist = set()
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line:
                whitelist.add(line.split("\t")[0])
    return whitelist


def correct_cb_one_mm(cb, whitelist):
    if cb in whitelist:
        return cb, "exact"
    bases = "ACGT"
    matches = []
    chars = list(cb)
    for pos, ref in enumerate(chars):
        for alt in bases:
            if alt == ref:
                continue
            variant = chars[:]
            variant[pos] = alt
            candidate = "".join(variant)
            if candidate in whitelist:
                matches.append(candidate)
                if len(matches) > 1:
                    return "", "ambiguous_1mm"
    if len(matches) == 1:
        return matches[0], "unique_1mm"
    return "", "no_match"


def iter_paired_fastq_records(r1_paths, r2_paths):
    if len(r1_paths) != len(r2_paths):
        raise SystemExit(f"R1/R2 path count mismatch: {len(r1_paths)} vs {len(r2_paths)}")
    for r1_path, r2_path in zip(r1_paths, r2_paths):
        with open_fastq(r1_path) as r1_handle, open_fastq(r2_path) as r2_handle:
            for r1_rec, r2_rec in zip(parse_fastq_handle(r1_handle), parse_fastq_handle(r2_handle)):
                q1 = normalize_qname(r1_rec[0])
                q2 = normalize_qname(r2_rec[0])
                if q1.split(" ", 1)[0] != q2.split(" ", 1)[0]:
                    raise SystemExit(f"FASTQ pair qname mismatch: {q1} vs {q2}")
                yield q2, r1_rec, r2_rec


def cmd_partition_fastq(args):
    _, idx_to_gene = read_probe_list_index(args.probe_list)
    whitelist = read_cb_whitelist(args.cb_whitelist)
    decisions = {}
    with open(args.scan_decisions_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            qname = normalize_qname(row["qname"])
            decisions[qname] = (row["decision"], int(row["resolved_gene_idx15"]))

    r1_paths = parse_csv_path_list(args.reads_r1)
    r2_paths = parse_csv_path_list(args.reads_r2)
    summary = defaultdict(int)
    os.makedirs(os.path.dirname(args.keep_records_tsv), exist_ok=True)
    pass_r1 = open_fastq_writer(args.pass_r1_fastq)
    pass_r2 = open_fastq_writer(args.pass_r2_fastq)
    keep_r1 = open_fastq_writer(args.keep_r1_fastq) if args.keep_r1_fastq else None
    keep_r2 = open_fastq_writer(args.keep_r2_fastq) if args.keep_r2_fastq else None
    with pass_r1, pass_r2, open(args.keep_records_tsv, "w", encoding="utf-8", newline="") as keep_handle:
        keep_writer = csv.writer(keep_handle, delimiter="\t")
        keep_writer.writerow(
            [
                "qname",
                "gene_idx15",
                "gene_id",
                "raw_cb",
                "corrected_cb",
                "cb_status",
                "umi",
                "sample_tag",
                "barcode_tag",
            ]
        )
        for qname, r1_rec, r2_rec in iter_paired_fastq_records(r1_paths, r2_paths):
            decision, gene_idx15 = decisions.get(qname, ("PASS", 0))
            _, r1_seq, _, _ = r1_rec
            _, r2_seq, _, _ = r2_rec
            if decision == "PASS":
                pass_r1.write("\n".join(r1_rec) + "\n")
                pass_r2.write("\n".join(r2_rec) + "\n")
                summary["pass_reads"] += 1
                continue
            if decision == "DENY":
                summary["deny_reads"] += 1
                continue
            if keep_r1 is not None and keep_r2 is not None:
                keep_r1.write("\n".join(r1_rec) + "\n")
                keep_r2.write("\n".join(r2_rec) + "\n")
            raw_cb = r1_seq[:16]
            umi = r1_seq[16:28]
            sample_tag = r2_seq[args.sample_offset:args.sample_offset + args.sample_length]
            corrected_cb, cb_status = correct_cb_one_mm(raw_cb, whitelist)
            barcode_tag = corrected_cb + sample_tag if corrected_cb else ""
            keep_writer.writerow(
                [
                    qname,
                    gene_idx15,
                    idx_to_gene.get(gene_idx15, ""),
                    raw_cb,
                    corrected_cb,
                    cb_status,
                    umi,
                    sample_tag,
                    barcode_tag,
                ]
            )
            summary["keep_reads"] += 1
            summary[f"keep_cb_{cb_status}"] += 1

    with open(args.summary_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for key in sorted(summary):
            writer.writerow([key, summary[key]])


def load_selected_matrix_counts(features_tsv, barcodes_tsv, matrix_mtx, selected_genes):
    barcodes = load_barcode_order(barcodes_tsv)
    selected_row_to_gene = {}
    with open(features_tsv, "r", encoding="utf-8") as handle:
        for idx, raw in enumerate(handle, start=1):
            parts = raw.rstrip("\n").split("\t")
            if parts and parts[0] in selected_genes:
                selected_row_to_gene[idx] = parts[0]

    counts = defaultdict(int)
    with open(matrix_mtx, "r", encoding="utf-8") as handle:
        saw_dimensions = False
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) == 3 and not saw_dimensions:
                saw_dimensions = True
                continue
            if len(parts) != 3:
                continue
            row_idx = int(parts[0])
            gene_id = selected_row_to_gene.get(row_idx)
            if gene_id is None:
                continue
            col_idx = int(parts[1])
            value = int(parts[2])
            counts[(gene_id, barcodes[col_idx - 1])] += value
    return counts


def load_keep_counts(keep_records_tsv, selected_genes):
    umi_sets = defaultdict(set)
    cb_status_counts = defaultdict(int)
    with open(keep_records_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cb_status_counts[row["cb_status"]] += 1
            if row["gene_id"] not in selected_genes:
                continue
            if not row["barcode_tag"]:
                continue
            umi_sets[(row["gene_id"], row["barcode_tag"])].add(row["umi"])

    counts = {(gene_id, barcode): len(umis) for (gene_id, barcode), umis in umi_sets.items()}
    return counts, cb_status_counts


def cmd_compare_keep_pass_matrix(args):
    selected_genes = load_selected_gene_ids(args.selected_probes)
    baseline_counts = load_selected_matrix_counts(
        args.baseline_features_tsv,
        args.baseline_barcodes_tsv,
        args.baseline_matrix_mtx,
        selected_genes,
    )
    pass_counts = load_selected_matrix_counts(
        args.pass_features_tsv,
        args.pass_barcodes_tsv,
        args.pass_matrix_mtx,
        selected_genes,
    )
    keep_counts, keep_cb_status = load_keep_counts(args.keep_records_tsv, selected_genes)

    combined = defaultdict(int)
    for key, value in pass_counts.items():
        combined[key] += value
    for key, value in keep_counts.items():
        combined[key] += value

    all_keys = sorted(set(baseline_counts) | set(combined))
    mismatches = []
    gene_totals = defaultdict(lambda: {"baseline": 0, "combined": 0, "pass": 0, "keep": 0})
    exact_entries = 0
    for gene_id, barcode in all_keys:
        key = (gene_id, barcode)
        base_value = baseline_counts.get(key, 0)
        pass_value = pass_counts.get(key, 0)
        keep_value = keep_counts.get(key, 0)
        combined_value = combined.get(key, 0)
        gene_totals[gene_id]["baseline"] += base_value
        gene_totals[gene_id]["combined"] += combined_value
        gene_totals[gene_id]["pass"] += pass_value
        gene_totals[gene_id]["keep"] += keep_value
        if base_value == combined_value:
            exact_entries += 1
        else:
            mismatches.append(
                {
                    "gene_id": gene_id,
                    "barcode": barcode,
                    "baseline_count": base_value,
                    "pass_count": pass_value,
                    "keep_count": keep_value,
                    "combined_count": combined_value,
                }
            )

    with open(args.diff_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "gene_id",
                "barcode",
                "baseline_count",
                "pass_count",
                "keep_count",
                "combined_count",
            ],
        )
        writer.writeheader()
        writer.writerows(mismatches)

    with open(args.gene_summary_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id", "baseline_total", "combined_total", "pass_total", "keep_total"])
        for gene_id in sorted(gene_totals):
            totals = gene_totals[gene_id]
            writer.writerow(
                [
                    gene_id,
                    totals["baseline"],
                    totals["combined"],
                    totals["pass"],
                    totals["keep"],
                ]
            )

    with open(args.summary_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["selected_genes", len(selected_genes)])
        writer.writerow(["baseline_nonzero_entries", len(baseline_counts)])
        writer.writerow(["pass_nonzero_entries", len(pass_counts)])
        writer.writerow(["keep_nonzero_entries", len(keep_counts)])
        writer.writerow(["exact_entries", exact_entries])
        writer.writerow(["mismatch_entries", len(mismatches)])
        writer.writerow(["keep_cb_exact", keep_cb_status.get("exact", 0)])
        writer.writerow(["keep_cb_unique_1mm", keep_cb_status.get("unique_1mm", 0)])
        writer.writerow(["keep_cb_ambiguous_1mm", keep_cb_status.get("ambiguous_1mm", 0)])
        writer.writerow(["keep_cb_no_match", keep_cb_status.get("no_match", 0)])


def cmd_compare_summed_matrices(args):
    selected_genes = load_selected_gene_ids(args.selected_probes)
    baseline_counts = load_selected_matrix_counts(
        args.baseline_features_tsv,
        args.baseline_barcodes_tsv,
        args.baseline_matrix_mtx,
        selected_genes,
    )
    left_counts = load_selected_matrix_counts(
        args.left_features_tsv,
        args.left_barcodes_tsv,
        args.left_matrix_mtx,
        selected_genes,
    )
    right_counts = load_selected_matrix_counts(
        args.right_features_tsv,
        args.right_barcodes_tsv,
        args.right_matrix_mtx,
        selected_genes,
    )

    combined = defaultdict(int)
    for key, value in left_counts.items():
        combined[key] += value
    for key, value in right_counts.items():
        combined[key] += value

    all_keys = sorted(set(baseline_counts) | set(combined))
    mismatches = []
    gene_totals = defaultdict(lambda: {"baseline": 0, "combined": 0, "left": 0, "right": 0})
    exact_entries = 0
    for gene_id, barcode in all_keys:
        key = (gene_id, barcode)
        base_value = baseline_counts.get(key, 0)
        left_value = left_counts.get(key, 0)
        right_value = right_counts.get(key, 0)
        combined_value = combined.get(key, 0)
        gene_totals[gene_id]["baseline"] += base_value
        gene_totals[gene_id]["combined"] += combined_value
        gene_totals[gene_id]["left"] += left_value
        gene_totals[gene_id]["right"] += right_value
        if base_value == combined_value:
            exact_entries += 1
        else:
            mismatches.append(
                {
                    "gene_id": gene_id,
                    "barcode": barcode,
                    "baseline_count": base_value,
                    "left_count": left_value,
                    "right_count": right_value,
                    "combined_count": combined_value,
                }
            )

    with open(args.diff_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "gene_id",
                "barcode",
                "baseline_count",
                "left_count",
                "right_count",
                "combined_count",
            ],
        )
        writer.writeheader()
        writer.writerows(mismatches)

    with open(args.gene_summary_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id", "baseline_total", "combined_total", "left_total", "right_total"])
        for gene_id in sorted(gene_totals):
            totals = gene_totals[gene_id]
            writer.writerow(
                [
                    gene_id,
                    totals["baseline"],
                    totals["combined"],
                    totals["left"],
                    totals["right"],
                ]
            )

    with open(args.summary_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["selected_genes", len(selected_genes)])
        writer.writerow(["baseline_nonzero_entries", len(baseline_counts)])
        writer.writerow(["left_nonzero_entries", len(left_counts)])
        writer.writerow(["right_nonzero_entries", len(right_counts)])
        writer.writerow(["exact_entries", exact_entries])
        writer.writerow(["mismatch_entries", len(mismatches)])


def cmd_compare_baseline(args):
    _, idx_to_gene = read_probe_list_index(args.probe_list)
    baseline_events = parse_reject_log(args.baseline_reject_log)
    baseline_state = {}
    for qname, rows in baseline_events.items():
        state, gene_idx15, _ = classify_qname_event(rows)
        baseline_state[qname] = (state, gene_idx15)

    scanned = load_scan_decisions(args.scan_decisions_tsv)
    mismatches = []
    matches = 0
    false_denies = 0
    keep_checks = 0
    deny_checks = 0
    for row in scanned:
        decision = row["decision"]
        if decision not in {"KEEP", "DENY"}:
            continue
        qname = normalize_qname(row["qname"])
        base_state, base_gene_idx15 = baseline_state.get(qname, ("DENY", 0))
        if decision == "KEEP":
            keep_checks += 1
            if base_state == "KEEP" and base_gene_idx15 == row["resolved_gene_idx15"]:
                matches += 1
            else:
                mismatches.append(
                    {
                        "qname": qname,
                        "cache_decision": decision,
                        "cache_gene_idx15": row["resolved_gene_idx15"],
                        "cache_gene_id": idx_to_gene.get(row["resolved_gene_idx15"], ""),
                        "baseline_state": base_state,
                        "baseline_gene_idx15": base_gene_idx15,
                        "baseline_gene_id": idx_to_gene.get(base_gene_idx15, ""),
                    }
                )
        else:
            deny_checks += 1
            if base_state == "KEEP":
                false_denies += 1
                mismatches.append(
                    {
                        "qname": qname,
                        "cache_decision": decision,
                        "cache_gene_idx15": 0,
                        "cache_gene_id": "",
                        "baseline_state": base_state,
                        "baseline_gene_idx15": base_gene_idx15,
                        "baseline_gene_id": idx_to_gene.get(base_gene_idx15, ""),
                    }
                )
            else:
                matches += 1

    with open(args.mismatches_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "qname",
                "cache_decision",
                "cache_gene_idx15",
                "cache_gene_id",
                "baseline_state",
                "baseline_gene_idx15",
                "baseline_gene_id",
            ],
        )
        writer.writeheader()
        writer.writerows(mismatches)

    with open(args.summary_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["keep_checks", keep_checks])
        writer.writerow(["deny_checks", deny_checks])
        writer.writerow(["matches", matches])
        writer.writerow(["mismatches", len(mismatches)])
        writer.writerow(["false_denies", false_denies])


def build_parser():
    parser = argparse.ArgumentParser(description="External 100-probe H0/H1 pilot for Flex.")
    sub = parser.add_subparsers(dest="command", required=True)

    select_p = sub.add_parser("select-probes")
    select_p.add_argument("--probe-list", required=True)
    select_p.add_argument("--probes-fasta", required=True)
    select_p.add_argument("--features-tsv", required=True)
    select_p.add_argument("--matrix-mtx", required=True)
    select_p.add_argument("--count", type=int, default=100)
    select_p.add_argument("--output", required=True)
    select_p.set_defaults(func=cmd_select_probes)

    synth_p = sub.add_parser("make-synth-fastq")
    synth_p.add_argument("--selected-probes", required=True)
    synth_p.add_argument("--cb-whitelist", required=True)
    synth_p.add_argument("--sample-whitelist")
    synth_p.add_argument("--sample-tag")
    synth_p.add_argument("--synthetic-umi", default="AAAAAAAAAAAA")
    synth_p.add_argument("--probe-offset", type=int, default=0)
    synth_p.add_argument("--sample-offset", type=int, default=68)
    synth_p.add_argument("--read-length", type=int, default=90)
    synth_p.add_argument("--fill-base", default="A")
    synth_p.add_argument("--left-flank", default="ACGT")
    synth_p.add_argument("--right-flank", default="TGCA")
    synth_p.add_argument("--r1-fastq", required=True)
    synth_p.add_argument("--r2-fastq", required=True)
    synth_p.add_argument("--manifest-tsv", required=True)
    synth_p.set_defaults(func=cmd_make_synth_fastq)

    cache_p = sub.add_parser("build-cache")
    cache_p.add_argument("--manifest-tsv", required=True)
    cache_p.add_argument("--reject-log", required=True)
    cache_p.add_argument("--qname-cache-tsv", required=True)
    cache_p.add_argument("--sequence-cache-tsv", required=True)
    cache_p.set_defaults(func=cmd_build_cache)

    mex_cache_p = sub.add_parser("build-cache-from-mex")
    mex_cache_p.add_argument("--manifest-tsv", required=True)
    mex_cache_p.add_argument("--probe-list", required=True)
    mex_cache_p.add_argument("--barcodes-tsv", required=True)
    mex_cache_p.add_argument("--features-tsv", required=True)
    mex_cache_p.add_argument("--matrix-mtx", required=True)
    mex_cache_p.add_argument("--qname-cache-tsv", required=True)
    mex_cache_p.add_argument("--sequence-cache-tsv", required=True)
    mex_cache_p.set_defaults(func=cmd_build_cache_from_mex)

    h2_fastq_p = sub.add_parser("h2-make-synth-fastq")
    h2_fastq_p.add_argument("--binary-cache", required=True, help="FH01SEQ1 cache; H0 rows (class 0) seed H2 variants")
    h2_fastq_p.add_argument("--probe-list", required=True)
    h2_fastq_p.add_argument("--cb-whitelist", required=True)
    h2_fastq_p.add_argument("--sample-whitelist")
    h2_fastq_p.add_argument("--sample-tag")
    h2_fastq_p.add_argument("--synthetic-umi", default="AAAAAAAAAAAA")
    h2_fastq_p.add_argument("--probe-offset", type=int, default=0)
    h2_fastq_p.add_argument("--sample-offset", type=int, default=68)
    h2_fastq_p.add_argument("--read-length", type=int, default=90)
    h2_fastq_p.add_argument("--fill-base", default="A")
    h2_fastq_p.add_argument("--left-flank", default="ACGT")
    h2_fastq_p.add_argument("--right-flank", default="TGCA")
    h2_fastq_p.add_argument("--num-shards", type=int, default=1)
    h2_fastq_p.add_argument("--shard-id", type=int, default=0)
    h2_fastq_p.add_argument(
        "--max-reads-per-shard",
        type=int,
        default=600_000,
        help="Must be <= CB whitelist size; raise num-shards if hit",
    )
    h2_fastq_p.add_argument("--parent-limit", type=int, default=None, help="Cap unique H0 sequences (testing)")
    h2_fastq_p.add_argument("--r1-fastq", required=True)
    h2_fastq_p.add_argument("--r2-fastq", required=True)
    h2_fastq_p.add_argument("--manifest-tsv", required=True)
    h2_fastq_p.set_defaults(func=cmd_h2_make_synth_fastq)

    h2_mex_p = sub.add_parser("h2-build-cache-from-mex")
    h2_mex_p.add_argument("--manifest-tsv", required=True)
    h2_mex_p.add_argument("--probe-list", required=True)
    h2_mex_p.add_argument("--barcodes-tsv", required=True)
    h2_mex_p.add_argument("--features-tsv", required=True)
    h2_mex_p.add_argument("--matrix-mtx", required=True)
    h2_mex_p.add_argument("--qname-cache-tsv", required=True)
    h2_mex_p.add_argument("--sequence-cache-tsv", required=True)
    h2_mex_p.set_defaults(func=cmd_h2_build_cache_from_mex)

    h2_bin_p = sub.add_parser("h2-write-binary-cache")
    h2_bin_p.add_argument("--sequence-cache-tsv", required=True)
    h2_bin_p.add_argument("--output-bin", required=True)
    h2_bin_p.add_argument("--sample-idx", type=int, default=0)
    h2_bin_p.add_argument("--keep-only", action="store_true", help="Only class-3 (H2 KEEP) rows")
    h2_bin_p.add_argument("--min-class", type=int, default=None)
    h2_bin_p.add_argument("--max-class", type=int, default=None)
    h2_bin_p.set_defaults(func=cmd_h2_write_binary_cache)

    scan_p = sub.add_parser("scan-fastq")
    scan_p.add_argument("--sequence-cache-tsv", required=True)
    scan_p.add_argument("--offsets", default="0,1")
    scan_p.add_argument("--output-tsv", required=True)
    scan_p.add_argument("--summary-tsv", required=True)
    scan_p.add_argument("fastq_r2", nargs="+")
    scan_p.set_defaults(func=cmd_scan_fastq)

    partition_p = sub.add_parser("partition-fastq")
    partition_p.add_argument("--scan-decisions-tsv", required=True)
    partition_p.add_argument("--probe-list", required=True)
    partition_p.add_argument("--cb-whitelist", required=True)
    partition_p.add_argument("--sample-offset", type=int, default=68)
    partition_p.add_argument("--sample-length", type=int, default=8)
    partition_p.add_argument("--reads-r1", required=True)
    partition_p.add_argument("--reads-r2", required=True)
    partition_p.add_argument("--pass-r1-fastq", required=True)
    partition_p.add_argument("--pass-r2-fastq", required=True)
    partition_p.add_argument("--keep-r1-fastq")
    partition_p.add_argument("--keep-r2-fastq")
    partition_p.add_argument("--keep-records-tsv", required=True)
    partition_p.add_argument("--summary-tsv", required=True)
    partition_p.set_defaults(func=cmd_partition_fastq)

    compare_matrix_p = sub.add_parser("compare-keep-pass-matrix")
    compare_matrix_p.add_argument("--selected-probes", required=True)
    compare_matrix_p.add_argument("--baseline-features-tsv", required=True)
    compare_matrix_p.add_argument("--baseline-barcodes-tsv", required=True)
    compare_matrix_p.add_argument("--baseline-matrix-mtx", required=True)
    compare_matrix_p.add_argument("--pass-features-tsv", required=True)
    compare_matrix_p.add_argument("--pass-barcodes-tsv", required=True)
    compare_matrix_p.add_argument("--pass-matrix-mtx", required=True)
    compare_matrix_p.add_argument("--keep-records-tsv", required=True)
    compare_matrix_p.add_argument("--summary-tsv", required=True)
    compare_matrix_p.add_argument("--diff-tsv", required=True)
    compare_matrix_p.add_argument("--gene-summary-tsv", required=True)
    compare_matrix_p.set_defaults(func=cmd_compare_keep_pass_matrix)

    compare_sum_p = sub.add_parser("compare-summed-matrices")
    compare_sum_p.add_argument("--selected-probes", required=True)
    compare_sum_p.add_argument("--baseline-features-tsv", required=True)
    compare_sum_p.add_argument("--baseline-barcodes-tsv", required=True)
    compare_sum_p.add_argument("--baseline-matrix-mtx", required=True)
    compare_sum_p.add_argument("--left-features-tsv", required=True)
    compare_sum_p.add_argument("--left-barcodes-tsv", required=True)
    compare_sum_p.add_argument("--left-matrix-mtx", required=True)
    compare_sum_p.add_argument("--right-features-tsv", required=True)
    compare_sum_p.add_argument("--right-barcodes-tsv", required=True)
    compare_sum_p.add_argument("--right-matrix-mtx", required=True)
    compare_sum_p.add_argument("--summary-tsv", required=True)
    compare_sum_p.add_argument("--diff-tsv", required=True)
    compare_sum_p.add_argument("--gene-summary-tsv", required=True)
    compare_sum_p.set_defaults(func=cmd_compare_summed_matrices)

    compare_p = sub.add_parser("compare-baseline")
    compare_p.add_argument("--probe-list", required=True)
    compare_p.add_argument("--scan-decisions-tsv", required=True)
    compare_p.add_argument("--baseline-reject-log", required=True)
    compare_p.add_argument("--summary-tsv", required=True)
    compare_p.add_argument("--mismatches-tsv", required=True)
    compare_p.set_defaults(func=cmd_compare_baseline)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""FASTQ preflight: detect chemistry, library type, and pair GEX/guide libraries.

Samples R1/R2 reads from each FASTQ set, detects barcode chemistry (TRU vs NXT),
infers library type (GEX vs GUIDE), normalizes barcodes to a canonical TRU
namespace, and pairs files by top-barcode Jaccard overlap.

Pipeline stages:
  1. Sample first N reads from each R1/R2 pair (--sample-reads, default 200K)
  2. Detect barcode chemistry (TRU vs NXT) from R1
  3. Detect library type (GEX vs GUIDE) from R2 guide/scaffold evidence
  4. Build barcode signatures (--top-barcodes, default 500)
  4b. Merge duplicate / multi-lane files (same type+chemistry, Jaccard > 0.20)
  5. Pairwise Jaccard scoring on merged libraries (cross-type)
  6. Connected-component detection and greedy bipartite matching
  7. Confidence tables (S/N ratio) and post-hoc name-based sanity checks

Tuning guidelines:
  The two key parameters are --sample-reads and --top-barcodes.

  --sample-reads  Controls how many R1/R2 reads are sampled per file.
                  More reads → more barcodes observed → denser signatures.
                  Default 200K is robust for 16-sample Perturb-seq with 10x
                  Chromium V3 chemistry. For datasets with very low barcode
                  diversity or very large numbers of samples, try 500K–1M.
                  Runtime scales linearly (200K ≈ 1s/file for gzipped FASTQs).

  --top-barcodes  Controls how many top-frequency barcodes form the signature.
                  Wider signatures capture more overlap between GEX and GUIDE
                  (which have different barcode rank distributions), but can
                  dilute the signal if pushed too far. Validated sweet spot:
                    50    — too narrow, misses cross-type overlap
                    200   — usually adequate for clean data
                    500   — robust for noisy / heterogeneous samples (default)
                    1000  — diminishing returns; only if S/N is still thin

  If the PAIRING CONFIDENCE table shows S/N below ~10x, increase both
  parameters (e.g. --sample-reads 500000 --top-barcodes 1000) and rerun.
  The merge confidence table is typically much stronger (S/N > 100x) because
  same-library duplicates share nearly identical barcode profiles.

See docs/RUNBOOK_FASTQ_PREFLIGHT_LIBRARY_PAIRING_20260317.md for the full spec.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import logging
import math
import os
import re
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

log = logging.getLogger("preflight")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_SAMPLE_READS = 200_000
BARCODE_LEN = 16
TOP_N_SIGNATURE = 500
TOP_N_DIAG = [50, 100, 500]

CHEMISTRY_MIN_RATE = 0.05
CHEMISTRY_DOMINANCE_RATIO = 3.0

GUIDE_EVIDENCE_THRESHOLD = 0.05
GUIDE_STRONG_THRESHOLD = 0.20

PAIRING_JACCARD_MIN = 0.02
PAIRING_MARGIN_RATIO = 2.0

MERGE_JACCARD_MIN = 0.20

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ChemistryCall:
    call: str  # TRU, NXT, AMBIGUOUS, UNKNOWN
    tru_exact: int = 0
    nxt_exact: int = 0
    total_scored: int = 0
    tru_rate: float = 0.0
    nxt_rate: float = 0.0
    confidence: float = 0.0


@dataclass
class LibraryTypeCall:
    call: str  # GEX, GUIDE, AMBIGUOUS, UNKNOWN
    guide_anchor_frac: float = 0.0
    guide_offset_histogram: Dict[int, int] = field(default_factory=dict)
    top_anchors: List[str] = field(default_factory=list)


@dataclass
class FileReport:
    label: str
    r1_path: str
    r2_path: str
    reads_sampled: int = 0
    chemistry: ChemistryCall = field(default_factory=lambda: ChemistryCall("UNKNOWN"))
    library_type: LibraryTypeCall = field(default_factory=lambda: LibraryTypeCall("UNKNOWN"))
    top_barcodes_50: List[Tuple[str, int]] = field(default_factory=list)
    top_barcodes_100: List[Tuple[str, int]] = field(default_factory=list)
    top_barcodes_500: List[Tuple[str, int]] = field(default_factory=list)


@dataclass
class MergedLibrary:
    label: str
    member_labels: List[str]
    chemistry: str
    library_type: str
    signature: collections.Counter = field(default_factory=collections.Counter)


@dataclass
class Pairing:
    file_a: str
    file_b: str
    jaccard: float
    weighted_jaccard: float
    top10_overlap: int
    top50_overlap: int
    chemistry_compatible: bool
    type_compatible: bool
    confidence: str  # safe, warning, manual_review


# ---------------------------------------------------------------------------
# FASTQ I/O
# ---------------------------------------------------------------------------

def open_fastq(path: Path):
    """Open a FASTQ, detecting gzip from magic bytes."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def sample_fastq(path: Path, n: int = DEFAULT_SAMPLE_READS) -> List[str]:
    """Return the first *n* sequences from a FASTQ file."""
    seqs: List[str] = []
    fh = open_fastq(path)
    try:
        while len(seqs) < n:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()  # +
            fh.readline()  # qual
            if seq:
                seqs.append(seq)
    finally:
        fh.close()
    return seqs

# ---------------------------------------------------------------------------
# Barcode namespace helpers
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement_base(base: str) -> str:
    return base.translate(_COMPLEMENT)


def translate_nxt_to_tru(barcode: str) -> str:
    """NXT→TRU: complement bases at 0-based positions 7 and 8."""
    if len(barcode) < 9:
        return barcode
    return barcode[:7] + complement_base(barcode[7:9]) + barcode[9:]


def load_whitelist(path: Path) -> Tuple[Set[str], Optional[Dict[str, str]]]:
    """Load a barcode whitelist.

    Single-column: returns (set_of_barcodes, None).
    Two-column (NXT translation): returns (set_of_col1, {col1: col2}).
    """
    barcodes: Set[str] = set()
    translation: Optional[Dict[str, str]] = None
    with open_fastq(path) as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            barcodes.add(parts[0])
            if len(parts) >= 2:
                if translation is None:
                    translation = {}
                translation[parts[0]] = parts[1]
    return barcodes, translation


# ---------------------------------------------------------------------------
# Chemistry detection
# ---------------------------------------------------------------------------

def detect_chemistry(
    r1_seqs: List[str],
    tru_whitelist: Set[str],
    nxt_whitelist: Set[str],
) -> ChemistryCall:
    """Score R1 barcodes against TRU and NXT whitelists."""
    tru_hits = 0
    nxt_hits = 0
    scored = 0

    for seq in r1_seqs:
        bc = seq[:BARCODE_LEN]
        if len(bc) < BARCODE_LEN:
            continue
        scored += 1
        if bc in tru_whitelist:
            tru_hits += 1
        if bc in nxt_whitelist:
            nxt_hits += 1

    if scored == 0:
        return ChemistryCall("UNKNOWN")

    tru_rate = tru_hits / scored
    nxt_rate = nxt_hits / scored

    if tru_rate < CHEMISTRY_MIN_RATE and nxt_rate < CHEMISTRY_MIN_RATE:
        call = "UNKNOWN"
        confidence = 0.0
    elif tru_rate >= CHEMISTRY_DOMINANCE_RATIO * max(nxt_rate, 1e-9):
        call = "TRU"
        confidence = tru_rate / (tru_rate + nxt_rate) if (tru_rate + nxt_rate) else 0
    elif nxt_rate >= CHEMISTRY_DOMINANCE_RATIO * max(tru_rate, 1e-9):
        call = "NXT"
        confidence = nxt_rate / (tru_rate + nxt_rate) if (tru_rate + nxt_rate) else 0
    else:
        call = "AMBIGUOUS"
        confidence = abs(tru_rate - nxt_rate) / (tru_rate + nxt_rate) if (tru_rate + nxt_rate) else 0

    return ChemistryCall(
        call=call,
        tru_exact=tru_hits,
        nxt_exact=nxt_hits,
        total_scored=scored,
        tru_rate=tru_rate,
        nxt_rate=nxt_rate,
        confidence=confidence,
    )


# ---------------------------------------------------------------------------
# Library type detection
# ---------------------------------------------------------------------------

def load_feature_sequences(feature_ref_path: Path) -> List[str]:
    """Extract guide sequences from a Cell Ranger feature reference CSV."""
    seqs: List[str] = []
    with open(feature_ref_path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            seq = row.get("sequence", "").strip().upper()
            if seq:
                seqs.append(seq)
    return seqs


_BASES = "ACGT"


def _build_hamming1_set(seqs: Set[str]) -> Set[str]:
    """Precompute all Hamming-distance-0-or-1 variants for a set of sequences."""
    h1: Set[str] = set(seqs)
    for s in seqs:
        arr = list(s)
        for i in range(len(arr)):
            orig = arr[i]
            for b in _BASES:
                if b != orig:
                    arr[i] = b
                    h1.add("".join(arr))
            arr[i] = orig
    return h1


def detect_library_type(
    r2_seqs: List[str],
    guide_seqs: Optional[List[str]] = None,
    scaffold_motifs: Optional[List[str]] = None,
) -> LibraryTypeCall:
    """Score R2 reads for guide-capture evidence."""

    if not r2_seqs:
        return LibraryTypeCall("UNKNOWN")

    guide_exact: Set[str] = set()
    guide_h1: Set[str] = set()
    guide_len = 0
    if guide_seqs:
        guide_exact = {s.upper() for s in guide_seqs}
        lengths = {len(s) for s in guide_exact}
        if len(lengths) == 1:
            guide_len = lengths.pop()
            guide_h1 = _build_hamming1_set(guide_exact)

    if scaffold_motifs is None:
        scaffold_motifs = []

    offset_histogram: Dict[int, int] = collections.defaultdict(int)
    anchor_hits = 0
    scaffold_hits = 0
    total = len(r2_seqs)

    for seq in r2_seqs:
        seq_upper = seq.upper()

        for motif in scaffold_motifs:
            if motif in seq_upper:
                scaffold_hits += 1
                break

        if guide_h1 and guide_len > 0:
            max_start = len(seq_upper) - guide_len
            for offset in range(0, min(max_start + 1, 60)):
                candidate = seq_upper[offset : offset + guide_len]
                if candidate in guide_h1:
                    anchor_hits += 1
                    offset_histogram[offset] += 1
                    break

    evidence_frac = (anchor_hits + scaffold_hits) / total if total else 0.0
    anchor_frac = anchor_hits / total if total else 0.0

    sorted_offsets = sorted(offset_histogram.items(), key=lambda x: -x[1])
    top_anchors = [f"offset={o} n={c}" for o, c in sorted_offsets[:5]]

    if evidence_frac >= GUIDE_STRONG_THRESHOLD:
        call = "GUIDE"
    elif evidence_frac >= GUIDE_EVIDENCE_THRESHOLD:
        call = "AMBIGUOUS"
    else:
        call = "GEX"

    return LibraryTypeCall(
        call=call,
        guide_anchor_frac=anchor_frac,
        guide_offset_histogram=dict(offset_histogram),
        top_anchors=top_anchors,
    )


# ---------------------------------------------------------------------------
# Barcode signature construction
# ---------------------------------------------------------------------------

def build_barcode_signature(
    r1_seqs: List[str],
    chemistry: ChemistryCall,
    tru_whitelist: Set[str],
    nxt_whitelist: Set[str],
    top_n: int = TOP_N_SIGNATURE,
) -> collections.Counter:
    """Correct barcodes, normalize to TRU, return top-N counter."""
    counter: collections.Counter = collections.Counter()

    for seq in r1_seqs:
        bc = seq[:BARCODE_LEN]
        if len(bc) < BARCODE_LEN:
            continue

        if chemistry.call == "NXT":
            bc_tru = translate_nxt_to_tru(bc)
        else:
            bc_tru = bc

        if bc_tru in tru_whitelist:
            counter[bc_tru] += 1

    return collections.Counter(dict(counter.most_common(top_n)))


# ---------------------------------------------------------------------------
# Jaccard and pairing
# ---------------------------------------------------------------------------

def jaccard_index(a: Set[str], b: Set[str]) -> float:
    union = a | b
    if not union:
        return 0.0
    return len(a & b) / len(union)


def weighted_jaccard(ca: collections.Counter, cb: collections.Counter) -> float:
    """Min-sum / max-sum weighted Jaccard."""
    all_keys = set(ca) | set(cb)
    if not all_keys:
        return 0.0
    min_sum = sum(min(ca.get(k, 0), cb.get(k, 0)) for k in all_keys)
    max_sum = sum(max(ca.get(k, 0), cb.get(k, 0)) for k in all_keys)
    return min_sum / max_sum if max_sum else 0.0


def rank_overlap(ca: collections.Counter, cb: collections.Counter, n: int) -> int:
    """Overlap count among top-n barcodes of each counter."""
    top_a = {bc for bc, _ in ca.most_common(n)}
    top_b = {bc for bc, _ in cb.most_common(n)}
    return len(top_a & top_b)


# ---------------------------------------------------------------------------
# Overlap graph and component detection
# ---------------------------------------------------------------------------

def build_overlap_graph(
    reports: List[FileReport],
    signatures: Dict[str, collections.Counter],
) -> List[Pairing]:
    """Pairwise Jaccard between all file reports; return scored edges."""
    pairings: List[Pairing] = []

    for i in range(len(reports)):
        for j in range(i + 1, len(reports)):
            ra, rb = reports[i], reports[j]
            sa = signatures[ra.label]
            sb = signatures[rb.label]

            j_idx = jaccard_index(set(sa), set(sb))
            w_idx = weighted_jaccard(sa, sb)
            t10 = rank_overlap(sa, sb, 10)
            t50 = rank_overlap(sa, sb, 50)

            chem_compat = True
            if ra.chemistry.call != "UNKNOWN" and rb.chemistry.call != "UNKNOWN":
                chem_compat = True  # already normalized to TRU

            type_compat = (
                {ra.library_type.call, rb.library_type.call} == {"GEX", "GUIDE"}
            )

            if j_idx >= PAIRING_JACCARD_MIN * PAIRING_MARGIN_RATIO:
                conf = "safe"
            elif j_idx >= PAIRING_JACCARD_MIN:
                conf = "warning"
            else:
                conf = "manual_review"

            pairings.append(Pairing(
                file_a=ra.label,
                file_b=rb.label,
                jaccard=j_idx,
                weighted_jaccard=w_idx,
                top10_overlap=t10,
                top50_overlap=t50,
                chemistry_compatible=chem_compat,
                type_compatible=type_compat,
                confidence=conf,
            ))

    return pairings


def find_connected_components(
    labels: List[str],
    pairings: List[Pairing],
) -> List[Set[str]]:
    """Connected components on edges above the minimum Jaccard threshold."""
    adj: Dict[str, Set[str]] = {lbl: set() for lbl in labels}
    for p in pairings:
        if p.jaccard >= PAIRING_JACCARD_MIN and p.chemistry_compatible:
            adj[p.file_a].add(p.file_b)
            adj[p.file_b].add(p.file_a)

    visited: Set[str] = set()
    components: List[Set[str]] = []

    for start in labels:
        if start in visited:
            continue
        comp: Set[str] = set()
        stack = [start]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            comp.add(node)
            stack.extend(adj[node] - visited)
        components.append(comp)

    return components


def match_within_component(
    comp: Set[str],
    report_map: Dict[str, FileReport],
    pairings: List[Pairing],
) -> List[Pairing]:
    """Greedy bipartite matching of GEX to GUIDE within a component."""
    gex = [lbl for lbl in comp if report_map[lbl].library_type.call in ("GEX", "AMBIGUOUS", "UNKNOWN")]
    guide = [lbl for lbl in comp if report_map[lbl].library_type.call == "GUIDE"]

    if not gex or not guide:
        return []

    pair_lookup: Dict[Tuple[str, str], Pairing] = {}
    for p in pairings:
        pair_lookup[(p.file_a, p.file_b)] = p
        pair_lookup[(p.file_b, p.file_a)] = p

    candidates = []
    for g in gex:
        for gu in guide:
            key = (g, gu)
            if key in pair_lookup:
                candidates.append(pair_lookup[key])
            key2 = (gu, g)
            if key2 in pair_lookup and key2 not in pair_lookup:
                candidates.append(pair_lookup[key2])

    candidates.sort(key=lambda p: -p.jaccard)

    used: Set[str] = set()
    matched: List[Pairing] = []
    for p in candidates:
        if p.file_a in used or p.file_b in used:
            continue
        if p.jaccard < PAIRING_JACCARD_MIN:
            continue
        used.add(p.file_a)
        used.add(p.file_b)
        matched.append(p)

    return matched


# ---------------------------------------------------------------------------
# Phase A: merge duplicate / multi-lane libraries
# ---------------------------------------------------------------------------

def merge_duplicate_libraries(
    reports: List[FileReport],
    signatures: Dict[str, collections.Counter],
    top_n: int = TOP_N_SIGNATURE,
) -> List[MergedLibrary]:
    """Identify same-type, same-chemistry files with high Jaccard and pool them.

    Returns one MergedLibrary per logical group, with pooled barcode counter
    trimmed to *top_n*.
    """
    n = len(reports)
    if n <= 1:
        return [_single_merged(r, signatures[r.label]) for r in reports]

    same_bucket: Dict[Tuple[str, str], List[int]] = collections.defaultdict(list)
    for i, r in enumerate(reports):
        key = (r.chemistry.call, r.library_type.call)
        same_bucket[key].append(i)

    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for bucket_indices in same_bucket.values():
        if len(bucket_indices) < 2:
            continue
        for ii in range(len(bucket_indices)):
            for jj in range(ii + 1, len(bucket_indices)):
                ai, bi = bucket_indices[ii], bucket_indices[jj]
                sa = signatures[reports[ai].label]
                sb = signatures[reports[bi].label]
                j = jaccard_index(set(sa), set(sb))
                if j >= MERGE_JACCARD_MIN:
                    union(ai, bi)

    groups: Dict[int, List[int]] = collections.defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)

    merged: List[MergedLibrary] = []
    for members in groups.values():
        member_reports = [reports[i] for i in members]
        member_labels = [r.label for r in member_reports]

        if len(members) == 1:
            r = member_reports[0]
            merged.append(_single_merged(r, signatures[r.label]))
            continue

        pooled = collections.Counter()
        for r in member_reports:
            pooled.update(signatures[r.label])
        pooled = collections.Counter(dict(pooled.most_common(top_n)))

        label = _merge_label(member_labels)
        merged.append(MergedLibrary(
            label=label,
            member_labels=member_labels,
            chemistry=member_reports[0].chemistry.call,
            library_type=member_reports[0].library_type.call,
            signature=pooled,
        ))

    return merged


def _single_merged(r: FileReport, sig: collections.Counter) -> MergedLibrary:
    return MergedLibrary(
        label=r.label,
        member_labels=[r.label],
        chemistry=r.chemistry.call,
        library_type=r.library_type.call,
        signature=sig,
    )


def _merge_label(labels: List[str]) -> str:
    """Construct a readable label for a group of merged files."""
    prefix = os.path.commonprefix(labels)
    if prefix and not prefix.endswith("_"):
        prefix = prefix.rsplit("_", 1)[0] + "_" if "_" in prefix else prefix
    if prefix and len(prefix) > 3:
        return prefix.rstrip("_") + f"_x{len(labels)}"
    return labels[0] + f"_x{len(labels)}"


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def _sanitize_for_json(obj):
    """Recursively convert an object tree so json.dump handles it cleanly.

    Fixes: tuple dict-keys from dataclasses.asdict, Counter objects,
    int dict-keys (JSON requires string keys), set/Path/tuple values.
    """
    if isinstance(obj, dict):
        return {str(k): _sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize_for_json(v) for v in obj]
    if isinstance(obj, set):
        return sorted(obj)
    if isinstance(obj, Path):
        return str(obj)
    return obj


def write_summary_tsv(reports: List[FileReport], path: Path) -> None:
    fields = [
        "label", "r1_path", "r2_path", "reads_sampled",
        "chemistry_call", "chemistry_confidence", "tru_rate", "nxt_rate",
        "library_type_call", "guide_anchor_frac",
        "top_guide_offsets",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for r in reports:
            top_offsets = ""
            if r.library_type.guide_offset_histogram:
                sorted_h = sorted(r.library_type.guide_offset_histogram.items(), key=lambda x: -x[1])
                top_offsets = ";".join(f"{o}:{c}" for o, c in sorted_h[:5])
            writer.writerow({
                "label": r.label,
                "r1_path": r.r1_path,
                "r2_path": r.r2_path,
                "reads_sampled": r.reads_sampled,
                "chemistry_call": r.chemistry.call,
                "chemistry_confidence": f"{r.chemistry.confidence:.4f}",
                "tru_rate": f"{r.chemistry.tru_rate:.4f}",
                "nxt_rate": f"{r.chemistry.nxt_rate:.4f}",
                "library_type_call": r.library_type.call,
                "guide_anchor_frac": f"{r.library_type.guide_anchor_frac:.4f}",
                "top_guide_offsets": top_offsets,
            })


def write_pairing_matrix_tsv(
    reports: List[FileReport],
    pairings: List[Pairing],
    path: Path,
) -> None:
    labels = [r.label for r in reports]
    jac: Dict[Tuple[str, str], float] = {}
    for p in pairings:
        jac[(p.file_a, p.file_b)] = p.jaccard
        jac[(p.file_b, p.file_a)] = p.jaccard

    with open(path, "w") as fh:
        fh.write("\t" + "\t".join(labels) + "\n")
        for la in labels:
            vals = []
            for lb in labels:
                if la == lb:
                    vals.append("1.0000")
                else:
                    vals.append(f"{jac.get((la, lb), 0.0):.4f}")
            fh.write(la + "\t" + "\t".join(vals) + "\n")


def write_pairing_report_json(
    reports: List[FileReport],
    pairings: List[Pairing],
    components: List[Set[str]],
    matched: List[Pairing],
    path: Path,
    extra: Optional[Dict] = None,
) -> None:
    report = {
        "files": [asdict(r) for r in reports],
        "pairwise_scores": [asdict(p) for p in pairings],
        "inferred_components": [sorted(c) for c in components],
        "matched_pairs": [asdict(p) for p in matched],
        "version": "0.2.0",
    }
    if extra:
        report.update(extra)
    with open(path, "w") as fh:
        json.dump(_sanitize_for_json(report), fh, indent=2)
        fh.write("\n")


# ---------------------------------------------------------------------------
# Multi-lane aggregation
# ---------------------------------------------------------------------------

_LANE_RE = re.compile(r"_L(\d{3})_")


def guess_logical_group(path: str) -> str:
    """Strip lane tag to group files from the same logical library."""
    return _LANE_RE.sub("_LXXX_", os.path.basename(path))


# ---------------------------------------------------------------------------
# Post-hoc sanity checks (name-based, after barcode-driven pairing)
# ---------------------------------------------------------------------------

_ILLUMINA_TOKENS_RE = re.compile(
    r"^(?:GEX_|guides_)?"           # strip our prefixed source-dir tag
    r"(?P<sample>.+?)"              # capture sample id
    r"(?:_S\d+)?(?:_L\d{3})?$"     # optional Illumina _S##, _L###
)


def _extract_sample_token(label: str) -> str:
    """Best-effort extraction of a sample identifier from a label.

    Strips common prefixes (GEX_, guides_) and Illumina suffixes (_S##, _L###)
    to get the core sample+index portion for comparison.
    """
    m = _ILLUMINA_TOKENS_RE.match(label)
    return m.group("sample") if m else label


def validate_preflight(
    reports: List[FileReport],
    merged_libs: List["MergedLibrary"],
    components: List[Set[str]],
    matched: List[Pairing],
    merged_report_map: Dict[str, FileReport],
) -> List[str]:
    """Run name-based and completeness sanity checks. Returns warning strings."""
    warnings: List[str] = []

    # 1. Empty / missing FASTQ files (zero reads sampled)
    for rpt in reports:
        if rpt.reads_sampled == 0:
            warnings.append(
                f"EMPTY_FASTQ: {rpt.label} yielded 0 reads "
                f"(R1={rpt.r1_path}, R2={rpt.r2_path})"
            )

    # 2. Merged groups: check filename consistency among members
    for ml in merged_libs:
        if len(ml.member_labels) < 2:
            continue
        tokens = {lbl: _extract_sample_token(lbl) for lbl in ml.member_labels}
        unique_tokens = set(tokens.values())
        if len(unique_tokens) > 1:
            warnings.append(
                f"NAME_MISMATCH: merged group '{ml.label}' has inconsistent "
                f"sample tokens: {dict(sorted(tokens.items()))}"
            )

    # 3. Paired components: flag if GEX or GUIDE side is empty
    for comp in components:
        types_in_comp = {merged_report_map[lbl].library_type.call for lbl in comp}
        comp_str = ", ".join(sorted(comp))
        if "GEX" not in types_in_comp:
            warnings.append(
                f"NO_GEX: component {{{comp_str}}} has no GEX library"
            )
        if "GUIDE" not in types_in_comp:
            warnings.append(
                f"NO_GUIDE: component {{{comp_str}}} has no GUIDE library"
            )

    # 4. Matched pairs: cross-check that paired labels share a sample token
    for m in matched:
        tok_a = _extract_sample_token(m.file_a)
        tok_b = _extract_sample_token(m.file_b)
        if tok_a != tok_b:
            warnings.append(
                f"PAIR_NAME_MISMATCH: paired {m.file_a} ({tok_a}) <-> "
                f"{m.file_b} ({tok_b}) — names suggest different samples"
            )

    return warnings


@dataclass
class ConfidenceRow:
    step: str          # "merge" or "pair"
    label: str
    partner: str
    best_jaccard: float
    second_jaccard: float
    sn_ratio: float
    margin: float


def compute_confidence_tables(
    merged_libs: List[MergedLibrary],
    individual_sigs: Dict[str, collections.Counter],
    pairings: List[Pairing],
    merged_report_map: Dict[str, FileReport],
    matched: List[Pairing],
) -> List[ConfidenceRow]:
    """Compute S/N confidence for both merge and pairing decisions."""
    rows: List[ConfidenceRow] = []

    # --- Merge confidence ---
    # For each multi-member group: the weakest intra-group Jaccard vs the
    # strongest same-type/chemistry external Jaccard.
    for ml in merged_libs:
        if len(ml.member_labels) < 2:
            continue
        intra_jaccards = []
        for i, a in enumerate(ml.member_labels):
            for b in ml.member_labels[i + 1:]:
                j = jaccard_index(set(individual_sigs[a]), set(individual_sigs[b]))
                intra_jaccards.append(j)
        min_intra = min(intra_jaccards)

        best_ext = 0.0
        for other in merged_libs:
            if other.label == ml.label:
                continue
            if other.library_type != ml.library_type:
                continue
            if other.chemistry != ml.chemistry:
                continue
            for a in ml.member_labels:
                for b in other.member_labels:
                    j = jaccard_index(
                        set(individual_sigs[a]), set(individual_sigs[b])
                    )
                    best_ext = max(best_ext, j)

        sn = min_intra / best_ext if best_ext > 0 else float("inf")
        rows.append(ConfidenceRow(
            step="merge",
            label=ml.label,
            partner=f"{len(ml.member_labels)} members",
            best_jaccard=min_intra,
            second_jaccard=best_ext,
            sn_ratio=sn,
            margin=min_intra - best_ext,
        ))

    # --- Pairing confidence ---
    # For each matched pair: best cross-type Jaccard vs second-best
    # cross-type Jaccard from any alternative partner.
    gex_labels = [
        lbl for lbl, r in merged_report_map.items()
        if r.library_type.call == "GEX"
    ]
    guide_labels = [
        lbl for lbl, r in merged_report_map.items()
        if r.library_type.call == "GUIDE"
    ]
    jac_lookup: Dict[Tuple[str, str], float] = {}
    for p in pairings:
        jac_lookup[(p.file_a, p.file_b)] = p.jaccard
        jac_lookup[(p.file_b, p.file_a)] = p.jaccard

    for m in matched:
        a, b = m.file_a, m.file_b
        best_j = m.jaccard
        if merged_report_map[a].library_type.call == "GEX":
            gex, guide = a, b
        else:
            gex, guide = b, a

        alt_gex = max(
            (jac_lookup.get((gex, g), 0.0) for g in guide_labels if g != guide),
            default=0.0,
        )
        alt_guide = max(
            (jac_lookup.get((guide, g), 0.0) for g in gex_labels if g != gex),
            default=0.0,
        )
        second_j = max(alt_gex, alt_guide)
        sn = best_j / second_j if second_j > 0 else float("inf")

        rows.append(ConfidenceRow(
            step="pair",
            label=gex,
            partner=guide,
            best_jaccard=best_j,
            second_jaccard=second_j,
            sn_ratio=sn,
            margin=best_j - second_j,
        ))

    return rows


def _print_confidence_table(rows: List[ConfidenceRow]) -> None:
    """Print a human-readable confidence table."""
    if not rows:
        return

    merge_rows = [r for r in rows if r.step == "merge"]
    pair_rows = [r for r in rows if r.step == "pair"]

    col = "{:<40s} {:>9s} {:>9s} {:>9s} {:>9s}"
    dat = "{:<40s} {:>9.4f} {:>9.4f} {:>9s} {:>9.4f}"

    def _sn_str(v: float) -> str:
        return "inf" if math.isinf(v) else f"{v:.1f}x"

    if merge_rows:
        print()
        print("MERGE CONFIDENCE (why files were grouped):")
        print(col.format("Group", "Intra-J", "Ext-J", "S/N", "Margin"))
        print("-" * 80)
        for r in merge_rows:
            label = r.label if len(r.label) <= 40 else r.label[:37] + "..."
            print(dat.format(label, r.best_jaccard, r.second_jaccard,
                             _sn_str(r.sn_ratio), r.margin))
        print(f"  Worst merge S/N: {_sn_str(min(r.sn_ratio for r in merge_rows))}")

    if pair_rows:
        print()
        print("PAIRING CONFIDENCE (why GEX matched to GUIDE):")
        print(col.format("GEX <-> GUIDE", "Best-J", "2nd-J", "S/N", "Margin"))
        print("-" * 80)
        for r in pair_rows:
            lbl = f"{r.label} <-> {r.partner}"
            if len(lbl) > 40:
                lbl = lbl[:37] + "..."
            print(dat.format(lbl, r.best_jaccard, r.second_jaccard,
                             _sn_str(r.sn_ratio), r.margin))
        print(f"  Worst pairing S/N: {_sn_str(min(r.sn_ratio for r in pair_rows))}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def parse_fastq_sets(raw_args: List[str]) -> List[Tuple[str, Path, Path]]:
    """Parse --fastqs label:R1:R2 arguments."""
    sets = []
    for i, arg in enumerate(raw_args):
        parts = arg.split(":")
        if len(parts) == 3:
            label, r1, r2 = parts
        elif len(parts) == 2:
            r1, r2 = parts
            label = f"lib{i}"
        else:
            raise ValueError(
                f"Expected label:R1:R2 or R1:R2, got {arg!r}"
            )
        sets.append((label, Path(r1), Path(r2)))
    return sets


def run_preflight(
    fastq_sets: List[Tuple[str, Path, Path]],
    tru_whitelist_path: Path,
    nxt_whitelist_path: Optional[Path],
    feature_ref_path: Optional[Path],
    sample_reads: int,
    top_n_barcodes: int,
    outdir: Path,
) -> int:
    """Run the full preflight pipeline. Returns 0 on success."""

    outdir.mkdir(parents=True, exist_ok=True)

    log.info("Loading TRU whitelist from %s", tru_whitelist_path)
    tru_wl, _ = load_whitelist(tru_whitelist_path)
    log.info("  %d barcodes", len(tru_wl))

    nxt_wl: Set[str] = set()
    if nxt_whitelist_path:
        log.info("Loading NXT whitelist from %s", nxt_whitelist_path)
        nxt_wl_raw, nxt_trans = load_whitelist(nxt_whitelist_path)
        if nxt_trans:
            nxt_wl = set(nxt_trans.keys())
        else:
            nxt_wl = nxt_wl_raw
        log.info("  %d barcodes", len(nxt_wl))
    else:
        log.info("No NXT whitelist provided; generating from TRU by complement rule")
        nxt_wl = {translate_nxt_to_tru(bc) for bc in tru_wl}
        # That gives us TRU->NXT translated *back*, but we actually need
        # the NXT surface. translate_nxt_to_tru is its own inverse.
        nxt_wl = {translate_nxt_to_tru(bc) for bc in tru_wl}

    guide_seqs: Optional[List[str]] = None
    scaffold_motifs: Optional[List[str]] = None
    if feature_ref_path:
        log.info("Loading feature reference from %s", feature_ref_path)
        guide_seqs = load_feature_sequences(feature_ref_path)
        log.info("  %d guide sequences", len(guide_seqs))
        scaffold_motifs = _extract_scaffold_motifs(feature_ref_path)

    # Phase 1-2: sample, detect chemistry, detect library type
    reports: List[FileReport] = []
    cached_r1: Dict[str, List[str]] = {}
    for label, r1, r2 in fastq_sets:
        log.info("Processing %s: R1=%s R2=%s", label, r1, r2)

        r1_seqs = sample_fastq(r1, sample_reads)
        r2_seqs = sample_fastq(r2, sample_reads)
        cached_r1[label] = r1_seqs
        n_sampled = min(len(r1_seqs), len(r2_seqs))
        log.info("  Sampled %d reads", n_sampled)

        chem = detect_chemistry(r1_seqs, tru_wl, nxt_wl)
        log.info("  Chemistry: %s (TRU=%.3f NXT=%.3f conf=%.3f)",
                 chem.call, chem.tru_rate, chem.nxt_rate, chem.confidence)

        ltype = detect_library_type(r2_seqs, guide_seqs, scaffold_motifs)
        log.info("  Library type: %s (anchor_frac=%.3f)",
                 ltype.call, ltype.guide_anchor_frac)

        rpt = FileReport(
            label=label,
            r1_path=str(r1),
            r2_path=str(r2),
            reads_sampled=n_sampled,
            chemistry=chem,
            library_type=ltype,
        )
        reports.append(rpt)

    # Phase 4: build barcode signatures (reuse cached R1 sequences)
    log.info("Building barcode signatures (top_n=%d, sample_reads=%d)", top_n_barcodes, sample_reads)
    signatures: Dict[str, collections.Counter] = {}
    for rpt in reports:
        r1_seqs = cached_r1[rpt.label]
        full_sig = build_barcode_signature(r1_seqs, rpt.chemistry, tru_wl, nxt_wl, top_n_barcodes)
        signatures[rpt.label] = full_sig
        rpt.top_barcodes_50 = full_sig.most_common(50)
        rpt.top_barcodes_100 = full_sig.most_common(100)
        rpt.top_barcodes_500 = full_sig.most_common(min(500, top_n_barcodes))

        n_bc = len(signatures[rpt.label])
        log.info("  %s: %d barcodes in top-%d signature", rpt.label, n_bc, top_n_barcodes)

    del cached_r1
    report_map_raw = {r.label: r for r in reports}

    # Phase 4b: merge duplicate libraries (same-type, same-chemistry,
    # high Jaccard) into logical groups before cross-type pairing.
    merged_libs = merge_duplicate_libraries(reports, signatures, top_n=top_n_barcodes)
    if len(merged_libs) < len(reports):
        log.info("Merged %d files into %d logical libraries:",
                 len(reports), len(merged_libs))
        for ml in merged_libs:
            log.info("  %s (%s/%s): %s",
                     ml.label, ml.chemistry, ml.library_type,
                     ml.member_labels)
    else:
        log.info("No duplicate libraries detected; %d files remain", len(reports))

    merged_reports: List[FileReport] = []
    merged_signatures: Dict[str, collections.Counter] = {}
    merged_report_map: Dict[str, FileReport] = {}
    for ml in merged_libs:
        first_member = report_map_raw[ml.member_labels[0]]
        rpt = FileReport(
            label=ml.label,
            r1_path=", ".join(report_map_raw[m].r1_path for m in ml.member_labels),
            r2_path=", ".join(report_map_raw[m].r2_path for m in ml.member_labels),
            reads_sampled=sum(report_map_raw[m].reads_sampled for m in ml.member_labels),
            chemistry=first_member.chemistry,
            library_type=first_member.library_type,
            top_barcodes_50=ml.signature.most_common(min(50, top_n_barcodes)),
            top_barcodes_100=ml.signature.most_common(min(100, top_n_barcodes)),
            top_barcodes_500=ml.signature.most_common(min(500, top_n_barcodes)),
        )
        merged_reports.append(rpt)
        merged_signatures[ml.label] = ml.signature
        merged_report_map[ml.label] = rpt

    # Phase 5: pairwise scoring on merged libraries
    pairings = build_overlap_graph(merged_reports, merged_signatures)
    for p in pairings:
        log.info("  %s <-> %s: Jaccard=%.4f wJaccard=%.4f top10=%d top50=%d type_compat=%s conf=%s",
                 p.file_a, p.file_b, p.jaccard, p.weighted_jaccard,
                 p.top10_overlap, p.top50_overlap, p.type_compatible, p.confidence)

    # Phase 6: overlap graph components and matching
    merged_labels = [r.label for r in merged_reports]
    components = find_connected_components(merged_labels, pairings)
    log.info("Inferred %d component(s): %s",
             len(components), [sorted(c) for c in components])

    all_matched: List[Pairing] = []
    for comp in components:
        matched = match_within_component(comp, merged_report_map, pairings)
        all_matched.extend(matched)
        for m in matched:
            log.info("  PAIRED: %s <-> %s (Jaccard=%.4f, confidence=%s)",
                     m.file_a, m.file_b, m.jaccard, m.confidence)

    # Phase 7: confidence metrics and post-hoc sanity checks
    confidence_rows = compute_confidence_tables(
        merged_libs, signatures, pairings, merged_report_map, all_matched,
    )
    validation_warnings = validate_preflight(
        reports, merged_libs, components, all_matched, merged_report_map,
    )
    for w in validation_warnings:
        log.warning("SANITY: %s", w)

    # Write outputs (individual-level reports for provenance, merged for pairing)
    summary_path = outdir / "preflight_summary.tsv"
    matrix_path = outdir / "preflight_pairing_matrix.tsv"
    report_path = outdir / "preflight_pairing_report.json"

    write_summary_tsv(reports, summary_path)
    write_pairing_matrix_tsv(merged_reports, pairings, matrix_path)

    merged_info = []
    for ml in merged_libs:
        d = asdict(ml)
        del d["signature"]
        merged_info.append(d)
    write_pairing_report_json(
        merged_reports, pairings, components, all_matched, report_path,
        extra={"parameters": {"sample_reads": sample_reads,
                               "top_n_barcodes": top_n_barcodes},
               "individual_files": [asdict(r) for r in reports],
               "merged_libraries": merged_info,
               "confidence": [asdict(c) for c in confidence_rows],
               "validation_warnings": validation_warnings},
    )

    log.info("Wrote %s", summary_path)
    log.info("Wrote %s", matrix_path)
    log.info("Wrote %s", report_path)

    # Print human-readable summary
    print()
    print("=" * 72)
    print("PREFLIGHT SUMMARY")
    print("=" * 72)

    print(f"  Parameters: sample_reads={sample_reads:,}  top_barcodes={top_n_barcodes}")
    print(f"  {len(reports)} individual files -> {len(merged_libs)} logical libraries")
    print()
    for rpt in reports:
        print(f"  {rpt.label}: chemistry={rpt.chemistry.call} "
              f"type={rpt.library_type.call} "
              f"reads={rpt.reads_sampled}")

    if len(merged_libs) < len(reports):
        print()
        print("MERGED LIBRARIES:")
        for ml in merged_libs:
            if len(ml.member_labels) > 1:
                print(f"  {ml.label} ({ml.library_type}/{ml.chemistry}): "
                      f"{' + '.join(ml.member_labels)}")

    print()
    if all_matched:
        print("PAIRED:")
        for m in all_matched:
            tag = "AUTO-PAIR" if m.confidence == "safe" else m.confidence.upper()
            print(f"  [{tag}] {m.file_a} <-> {m.file_b}  "
                  f"Jaccard={m.jaccard:.4f}")
    else:
        print("NO CONFIDENT PAIRINGS FOUND")

    orphans = set(merged_labels)
    for m in all_matched:
        orphans.discard(m.file_a)
        orphans.discard(m.file_b)
    if orphans:
        print()
        print("UNPAIRED:")
        for o in sorted(orphans):
            r = merged_report_map[o]
            print(f"  {o}: chemistry={r.chemistry.call} "
                  f"type={r.library_type.call}")

    _print_confidence_table(confidence_rows)

    if validation_warnings:
        print()
        print("SANITY WARNINGS:")
        for w in validation_warnings:
            print(f"  ⚠ {w}")
    else:
        print()
        print("SANITY CHECKS: all passed")

    print()
    print(f"Full report: {report_path}")
    print("=" * 72)

    return 0


def _extract_scaffold_motifs(feature_ref_path: Path) -> List[str]:
    """Extract scaffold motifs from the pattern column of a feature reference."""
    motifs: List[str] = []
    try:
        with open(feature_ref_path, "r", newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                pattern = row.get("pattern", "").strip()
                if not pattern:
                    continue
                # Pattern format: 5P(BC)SCAFFOLD or (BC)SCAFFOLD
                # Extract the constant part after (BC)
                m = re.search(r"\(BC\)(.+)", pattern, re.IGNORECASE)
                if m:
                    scaffold = m.group(1).upper().replace("N", "")
                    if len(scaffold) >= 8:
                        motifs.append(scaffold)
    except (KeyError, csv.Error):
        pass
    # Deduplicate
    return list(dict.fromkeys(motifs))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="FASTQ preflight: detect chemistry and pair GEX/guide libraries.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:

  # Two FASTQ sets, TRU whitelist only:
  %(prog)s \\
    --fastqs GEX:/data/GEX_R1.fastq.gz:/data/GEX_R2.fastq.gz \\
    --fastqs GUIDE:/data/Guide_R1.fastq.gz:/data/Guide_R2.fastq.gz \\
    --tru-whitelist /path/to/3M-february-2018.txt \\
    --outdir /tmp/preflight

  # With NXT whitelist and feature reference:
  %(prog)s \\
    --fastqs GEX:/data/GEX_R1.fastq.gz:/data/GEX_R2.fastq.gz \\
    --fastqs GUIDE:/data/Guide_R1.fastq.gz:/data/Guide_R2.fastq.gz \\
    --tru-whitelist /path/to/3M-february-2018.txt \\
    --nxt-whitelist /path/to/3M-february-2018_NXT.txt \\
    --feature-ref /path/to/feature_reference.csv \\
    --outdir /tmp/preflight

  # Increase sensitivity if S/N is low:
  %(prog)s \\
    --fastqs ... \\
    --tru-whitelist ... \\
    --sample-reads 500000 --top-barcodes 1000 \\
    --outdir /tmp/preflight_deep

Tuning:
  If the PAIRING CONFIDENCE table shows S/N < 10, increase --sample-reads
  and/or --top-barcodes. Defaults (200K reads, 500 barcodes) are robust for
  typical 10x Chromium V3 Perturb-seq (validated on 16-sample UCSF dataset:
  worst-case pairing S/N = 50x, merge S/N = 288x).
""",
    )

    parser.add_argument(
        "--fastqs", action="append", required=True, metavar="LABEL:R1:R2",
        help="FASTQ set as label:R1_path:R2_path (may be repeated)",
    )
    parser.add_argument(
        "--tru-whitelist", required=True, type=Path,
        help="Path to the TRU barcode whitelist (single-column)",
    )
    parser.add_argument(
        "--nxt-whitelist", type=Path, default=None,
        help="Path to the NXT barcode whitelist or NXT->TRU translation file "
             "(two-column). If omitted, NXT surface is derived from TRU by "
             "complement rule.",
    )
    parser.add_argument(
        "--feature-ref", type=Path, default=None,
        help="Path to Cell Ranger feature reference CSV (for guide detection)",
    )
    parser.add_argument(
        "--sample-reads", type=int, default=DEFAULT_SAMPLE_READS,
        help="Number of reads to sample from each FASTQ "
             f"(default: {DEFAULT_SAMPLE_READS:,}). Increase if pairing S/N "
             "is low; 200K–500K is usually sufficient.",
    )
    parser.add_argument(
        "--top-barcodes", type=int, default=TOP_N_SIGNATURE,
        help="Number of top barcodes to retain per library for pairing "
             f"(default: {TOP_N_SIGNATURE}). Use 200–500 for typical "
             "Perturb-seq; increase to 1000 if pairing S/N is thin.",
    )
    parser.add_argument(
        "--outdir", type=Path, required=True,
        help="Output directory for preflight reports",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    try:
        fastq_sets = parse_fastq_sets(args.fastqs)
    except ValueError as exc:
        parser.error(str(exc))

    for _, r1, r2 in fastq_sets:
        if not r1.exists():
            parser.error(f"R1 not found: {r1}")
        if not r2.exists():
            parser.error(f"R2 not found: {r2}")
    if not args.tru_whitelist.exists():
        parser.error(f"TRU whitelist not found: {args.tru_whitelist}")
    if args.nxt_whitelist and not args.nxt_whitelist.exists():
        parser.error(f"NXT whitelist not found: {args.nxt_whitelist}")
    if args.feature_ref and not args.feature_ref.exists():
        parser.error(f"Feature reference not found: {args.feature_ref}")

    return run_preflight(
        fastq_sets=fastq_sets,
        tru_whitelist_path=args.tru_whitelist,
        nxt_whitelist_path=args.nxt_whitelist,
        feature_ref_path=args.feature_ref,
        sample_reads=args.sample_reads,
        top_n_barcodes=args.top_barcodes,
        outdir=args.outdir,
    )


if __name__ == "__main__":
    raise SystemExit(main())

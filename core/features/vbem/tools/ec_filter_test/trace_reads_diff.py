#!/usr/bin/env python3
"""
Trace-reads diff harness to compare Salmon vs CLI decisions per read (alignment-mode).

This script compares trace outputs from Salmon and the EC filter CLI to pinpoint
where alignment-mode decisions diverge. It tracks:
- Transcript IDs (txpIDs)
- Alignment scores (as)
- Best alignment score (bestAS)
- Log fragment probability (logFragProb)
- Error likelihood (errLike)
- Log compatibility probability (logCompat)
- Orphan status
- Auxiliary probabilities (auxProb)
- Library format fields (obsFmt, expFmt, isCompat, mateFields, mateStatus)

Usage:
    python3 trace_reads_diff.py --salmon-trace <salmon_trace.txt> --cli-trace <cli_trace.txt> --output <diff_report.txt>
"""

import argparse
import sys
from collections import defaultdict
from typing import Dict, List, Tuple, Optional


def _parse_trace_line(line: str) -> Tuple[str, Dict]:
    """Parse a single trace line and return (qname, fields_dict)."""
    line = line.strip()
    if not line or line.startswith('#') or '\t' not in line:
        return None, None
    
    qname, payload = line.split('\t', 1)
    fields = {}
    for field in payload.split(';'):
        if '=' not in field:
            continue
        key, val = field.split('=', 1)
        fields[key] = val
    
    return qname, fields


def _build_alignment_records(fields: Dict) -> List[Dict]:
    """Build a list of per-alignment records from parsed fields.
    
    Each record contains: txpID, as, logFragProb, errLike, logCompat, 
    obsFmt, expFmt, isCompat, mateFields, mateStatus, orphan, auxProb, pos, matePos, fwd, mateFwd
    """
    records = []
    
    # Parse all per-alignment arrays
    txp_ids = [int(x) for x in fields.get('txpIDs', '').split(',') if x]
    as_scores = [int(x) for x in fields.get('as', '').split(',') if x]
    log_frag = [float(x) for x in fields.get('logFragProb', '').split(',') if x]
    err_like = [float(x) for x in fields.get('errLike', '').split(',') if x]
    log_compat = [float(x) for x in fields.get('logCompat', '').split(',') if x]
    orphan = [bool(int(x)) for x in fields.get('orphan', '').split(',') if x]
    aux_prob = [float(x) for x in fields.get('auxProb', '').split(',') if x]
    dropped = [bool(int(x)) for x in fields.get('dropped', '').split(',') if x]
    
    # Library format fields
    obs_fmt = [int(x) for x in fields.get('obsFmt', '').split(',') if x]
    exp_fmt = [int(x) for x in fields.get('expFmt', '').split(',') if x]
    is_compat = [bool(int(x)) for x in fields.get('isCompat', '').split(',') if x]
    mate_fields = [bool(int(x)) for x in fields.get('mateFields', '').split(',') if x]
    mate_status = [int(x) for x in fields.get('mateStatus', '').split(',') if x]
    
    # Position fields
    pos = [int(x) for x in fields.get('pos', '').split(',') if x]
    mate_pos = [int(x) for x in fields.get('matePos', '').split(',') if x]
    fwd = [bool(int(x)) for x in fields.get('fwd', '').split(',') if x]
    mate_fwd = [bool(int(x)) for x in fields.get('mateFwd', '').split(',') if x]
    
    # Determine which indices to keep (filter out dropped alignments)
    if dropped:
        keep_idx = [i for i, d in enumerate(dropped) if not d]
    else:
        keep_idx = list(range(len(txp_ids)))
    
    for idx in keep_idx:
        record = {
            'txpID': txp_ids[idx] if idx < len(txp_ids) else None,
            'as': as_scores[idx] if idx < len(as_scores) else None,
            'logFragProb': log_frag[idx] if idx < len(log_frag) else None,
            'errLike': err_like[idx] if idx < len(err_like) else None,
            'logCompat': log_compat[idx] if idx < len(log_compat) else None,
            'orphan': orphan[idx] if idx < len(orphan) else False,
            'auxProb': aux_prob[idx] if idx < len(aux_prob) else None,
            'obsFmt': obs_fmt[idx] if idx < len(obs_fmt) else None,
            'expFmt': exp_fmt[idx] if idx < len(exp_fmt) else None,
            'isCompat': is_compat[idx] if idx < len(is_compat) else None,
            'mateFields': mate_fields[idx] if idx < len(mate_fields) else None,
            'mateStatus': mate_status[idx] if idx < len(mate_status) else None,
            'pos': pos[idx] if idx < len(pos) else None,
            'matePos': mate_pos[idx] if idx < len(mate_pos) else None,
            'fwd': fwd[idx] if idx < len(fwd) else None,
            'mateFwd': mate_fwd[idx] if idx < len(mate_fwd) else None,
        }
        records.append(record)
    
    return records


def _parse_trace_file_detailed(trace_file: str) -> Dict[str, Dict]:
    """Parse trace file and return per-qname data with per-alignment records."""
    reads = {}
    with open(trace_file, 'r') as f:
        for line in f:
            qname, fields = _parse_trace_line(line)
            if qname is None:
                continue
            
            # Build alignment records
            records = _build_alignment_records(fields)
            
            # Parse bestAS
            best_as = fields.get('bestAS')
            best_as = int(best_as) if best_as else None
            
            # Store raw line for debugging
            raw_line = line.strip()
            
            reads[qname] = {
                'qname': qname,
                'bestAS': best_as,
                'records': records,
                'raw_line': raw_line,
                # Summarized fields for compatibility comparison
                'txpIDs': sorted(set(r['txpID'] for r in records if r['txpID'] is not None)),
                'auxProb': [r['auxProb'] for r in records],
                'logFragProb': [r['logFragProb'] for r in records],
                'errLike': [r['errLike'] for r in records],
                'logCompat': [r['logCompat'] for r in records],
            }
    
    return reads


def _list_close(a: List[Optional[float]], b: List[Optional[float]], tol: float) -> bool:
    if len(a) != len(b):
        return False
    for i in range(len(a)):
        if a[i] is None and b[i] is None:
            continue
        if a[i] is None or b[i] is None:
            return False
        if abs(a[i] - b[i]) > tol:
            return False
    return True


def _compute_component_deltas(salmon_read: Dict, cli_read: Dict) -> Dict:
    """Compute deltas between components to identify which explains the auxProb mismatch."""
    deltas = {
        'auxProb': [],
        'logFragProb': [],
        'errLike': [],
        'logCompat': [],
    }
    
    salmon_aux = salmon_read.get('auxProb', [])
    cli_aux = cli_read.get('auxProb', [])
    salmon_frag = salmon_read.get('logFragProb', [])
    cli_frag = cli_read.get('logFragProb', [])
    salmon_err = salmon_read.get('errLike', [])
    cli_err = cli_read.get('errLike', [])
    salmon_compat = salmon_read.get('logCompat', [])
    cli_compat = cli_read.get('logCompat', [])
    
    # Compute per-alignment deltas
    n = max(len(salmon_aux), len(cli_aux), 1)
    for i in range(n):
        aux_s = salmon_aux[i] if i < len(salmon_aux) and salmon_aux[i] is not None else 0.0
        aux_c = cli_aux[i] if i < len(cli_aux) and cli_aux[i] is not None else 0.0
        deltas['auxProb'].append(aux_c - aux_s)
        
        frag_s = salmon_frag[i] if i < len(salmon_frag) and salmon_frag[i] is not None else 0.0
        frag_c = cli_frag[i] if i < len(cli_frag) and cli_frag[i] is not None else 0.0
        deltas['logFragProb'].append(frag_c - frag_s)
        
        err_s = salmon_err[i] if i < len(salmon_err) and salmon_err[i] is not None else 0.0
        err_c = cli_err[i] if i < len(cli_err) and cli_err[i] is not None else 0.0
        deltas['errLike'].append(err_c - err_s)
        
        compat_s = salmon_compat[i] if i < len(salmon_compat) and salmon_compat[i] is not None else 0.0
        compat_c = cli_compat[i] if i < len(cli_compat) and cli_compat[i] is not None else 0.0
        deltas['logCompat'].append(compat_c - compat_s)
    
    return deltas


def _identify_mismatch_cause(deltas: Dict, tol: float = 1e-4) -> str:
    """Identify which component(s) explain the auxProb mismatch."""
    causes = []
    
    # Check if logFragProb delta explains auxProb delta
    aux_delta = sum(abs(d) for d in deltas['auxProb']) / max(len(deltas['auxProb']), 1)
    frag_delta = sum(abs(d) for d in deltas['logFragProb']) / max(len(deltas['logFragProb']), 1)
    err_delta = sum(abs(d) for d in deltas['errLike']) / max(len(deltas['errLike']), 1)
    compat_delta = sum(abs(d) for d in deltas['logCompat']) / max(len(deltas['logCompat']), 1)
    
    if frag_delta > tol:
        causes.append(f"logFragProb (delta={frag_delta:.4f})")
    if err_delta > tol:
        causes.append(f"errLike (delta={err_delta:.4f})")
    if compat_delta > tol:
        causes.append(f"logCompat (delta={compat_delta:.4f})")
    
    if not causes:
        causes.append("unknown (all component deltas < tolerance)")
    
    return "; ".join(causes)


def compare_reads_detailed(salmon_reads: Dict[str, Dict], cli_reads: Dict[str, Dict]) -> Tuple[List[Dict], Dict]:
    """Compare Salmon and CLI reads with detailed component analysis."""
    differences = []
    component_stats = defaultdict(int)
    
    all_qnames = set(salmon_reads.keys()) | set(cli_reads.keys())
    
    for qname in sorted(all_qnames):
        salmon_read = salmon_reads.get(qname)
        cli_read = cli_reads.get(qname)
        
        if not salmon_read:
            differences.append({
                'qname': qname,
                'type': 'missing_in_salmon',
                'details': 'Read present in CLI but not in Salmon'
            })
            component_stats['missing_in_salmon'] += 1
            continue
        
        if not cli_read:
            differences.append({
                'qname': qname,
                'type': 'missing_in_cli',
                'details': 'Read present in Salmon but not in CLI'
            })
            component_stats['missing_in_cli'] += 1
            continue
        
        # Compare bestAS
        if salmon_read.get('bestAS') != cli_read.get('bestAS'):
            differences.append({
                'qname': qname,
                'type': 'bestAS_mismatch',
                'salmon': salmon_read.get('bestAS'),
                'cli': cli_read.get('bestAS'),
                'details': f"BestAS mismatch: Salmon={salmon_read.get('bestAS')}, CLI={cli_read.get('bestAS')}"
            })
            component_stats['bestAS_mismatch'] += 1
        
        # Compare txpIDs
        if salmon_read.get('txpIDs') != cli_read.get('txpIDs'):
            differences.append({
                'qname': qname,
                'type': 'txpIDs_mismatch',
                'salmon': salmon_read.get('txpIDs'),
                'cli': cli_read.get('txpIDs'),
                'details': f"Transcript IDs mismatch: Salmon={salmon_read.get('txpIDs')}, CLI={cli_read.get('txpIDs')}"
            })
            component_stats['txpIDs_mismatch'] += 1
        
        # Compare auxProb with tolerance (floating precision)
        if not _list_close(salmon_read.get('auxProb', []), cli_read.get('auxProb', []), tol=1e-6):
            deltas = _compute_component_deltas(salmon_read, cli_read)
            cause = _identify_mismatch_cause(deltas)
            
            diff = {
                'qname': qname,
                'type': 'auxProb_mismatch',
                'salmon': salmon_read.get('auxProb'),
                'cli': cli_read.get('auxProb'),
                'details': f"Auxiliary probability mismatch: Salmon={salmon_read.get('auxProb')}, CLI={cli_read.get('auxProb')}",
                'deltas': deltas,
                'cause': cause,
                'salmon_read': salmon_read,
                'cli_read': cli_read,
            }
            differences.append(diff)
            component_stats['auxProb_mismatch'] += 1
            
            # Track cause breakdown
            if 'logFragProb' in cause:
                component_stats['auxProb_caused_by_logFragProb'] += 1
            if 'errLike' in cause:
                component_stats['auxProb_caused_by_errLike'] += 1
            if 'logCompat' in cause:
                component_stats['auxProb_caused_by_logCompat'] += 1
            if 'unknown' in cause:
                component_stats['auxProb_cause_unknown'] += 1
        
        # Compare logCompat directly
        if not _list_close(salmon_read.get('logCompat', []), cli_read.get('logCompat', []), tol=1e-6):
            component_stats['logCompat_direct_mismatch'] += 1
        
        # Compare logFragProb directly  
        if not _list_close(salmon_read.get('logFragProb', []), cli_read.get('logFragProb', []), tol=1e-6):
            component_stats['logFragProb_direct_mismatch'] += 1
        
        # Compare errLike directly
        if not _list_close(salmon_read.get('errLike', []), cli_read.get('errLike', []), tol=1e-6):
            component_stats['errLike_direct_mismatch'] += 1
    
    return differences, component_stats


def write_detailed_report(differences: List[Dict], component_stats: Dict, output_file: str):
    """Write enhanced comparison report with component analysis."""
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("Salmon vs CLI Trace Comparison Report (Enhanced)\n")
        f.write("=" * 80 + "\n\n")
        
        if not differences:
            f.write("✅ No differences found! All reads match.\n")
            return
        
        # Group differences by type
        by_type = defaultdict(list)
        for diff in differences:
            by_type[diff['type']].append(diff)
        
        f.write(f"Total differences: {len(differences)}\n")
        f.write(f"Unique reads with differences: {len(set(d['qname'] for d in differences))}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("Difference Summary by Type:\n")
        f.write("=" * 80 + "\n")
        for diff_type, diffs in sorted(by_type.items()):
            f.write(f"  {diff_type}: {len(diffs)}\n")
        f.write("\n")
        
        f.write("=" * 80 + "\n")
        f.write("Component Statistics:\n")
        f.write("=" * 80 + "\n")
        for stat, count in sorted(component_stats.items()):
            f.write(f"  {stat}: {count}\n")
        f.write("\n")
        
        # First divergence
        if differences:
            first_diff = differences[0]
            f.write("=" * 80 + "\n")
            f.write("FIRST DIVERGENCE:\n")
            f.write("=" * 80 + "\n")
            f.write(f"Read: {first_diff['qname']}\n")
            f.write(f"Type: {first_diff['type']}\n")
            f.write(f"Details: {first_diff['details']}\n")
            if 'cause' in first_diff:
                f.write(f"Cause: {first_diff['cause']}\n")
            if 'salmon' in first_diff:
                f.write(f"Salmon value: {first_diff['salmon']}\n")
            if 'cli' in first_diff:
                f.write(f"CLI value: {first_diff['cli']}\n")
            f.write("\n")
        
        # Top 5 auxProb mismatches with full context
        auxprob_mismatches = [d for d in differences if d['type'] == 'auxProb_mismatch']
        if auxprob_mismatches:
            f.write("=" * 80 + "\n")
            f.write("TOP 5 auxProb MISMATCHES WITH FULL CONTEXT:\n")
            f.write("=" * 80 + "\n\n")
            
            for i, diff in enumerate(auxprob_mismatches[:5], 1):
                f.write(f"--- Mismatch {i}: {diff['qname']} ---\n")
                f.write(f"Cause: {diff.get('cause', 'unknown')}\n")
                f.write(f"auxProb delta: {diff.get('deltas', {}).get('auxProb', [])}\n")
                f.write(f"logFragProb delta: {diff.get('deltas', {}).get('logFragProb', [])}\n")
                f.write(f"errLike delta: {diff.get('deltas', {}).get('errLike', [])}\n")
                f.write(f"logCompat delta: {diff.get('deltas', {}).get('logCompat', [])}\n\n")
                
                # Per-alignment breakdown
                salmon_read = diff.get('salmon_read', {})
                cli_read = diff.get('cli_read', {})
                
                f.write("Per-alignment breakdown:\n")
                salmon_records = salmon_read.get('records', [])
                cli_records = cli_read.get('records', [])
                
                f.write("  Salmon alignments:\n")
                for j, rec in enumerate(salmon_records[:10]):  # Limit to first 10
                    f.write(f"    [{j}] txpID={rec.get('txpID')}, pos={rec.get('pos')}, matePos={rec.get('matePos')}, ")
                    f.write(f"logFragProb={rec.get('logFragProb')}, errLike={rec.get('errLike')}, ")
                    f.write(f"logCompat={rec.get('logCompat')}, auxProb={rec.get('auxProb')}, ")
                    f.write(f"obsFmt={rec.get('obsFmt')}, expFmt={rec.get('expFmt')}, ")
                    f.write(f"isCompat={rec.get('isCompat')}, mateFields={rec.get('mateFields')}\n")
                
                f.write("  CLI alignments:\n")
                for j, rec in enumerate(cli_records[:10]):  # Limit to first 10
                    f.write(f"    [{j}] txpID={rec.get('txpID')}, pos={rec.get('pos')}, matePos={rec.get('matePos')}, ")
                    f.write(f"logFragProb={rec.get('logFragProb')}, errLike={rec.get('errLike')}, ")
                    f.write(f"logCompat={rec.get('logCompat')}, auxProb={rec.get('auxProb')}, ")
                    f.write(f"obsFmt={rec.get('obsFmt')}, expFmt={rec.get('expFmt')}, ")
                    f.write(f"isCompat={rec.get('isCompat')}, mateFields={rec.get('mateFields')}\n")
                
                f.write("\n  Raw Salmon trace line:\n")
                f.write(f"    {salmon_read.get('raw_line', 'N/A')[:500]}...\n")
                f.write("\n  Raw CLI trace line:\n")
                f.write(f"    {cli_read.get('raw_line', 'N/A')[:500]}...\n")
                f.write("\n")
        
        # Detailed differences (first 50)
        f.write("=" * 80 + "\n")
        f.write("Detailed Differences (first 50):\n")
        f.write("=" * 80 + "\n\n")
        
        for i, diff in enumerate(differences[:50], 1):
            f.write(f"{i}. Read: {diff['qname']}\n")
            f.write(f"   Type: {diff['type']}\n")
            f.write(f"   Details: {diff['details']}\n")
            if 'cause' in diff:
                f.write(f"   Cause: {diff['cause']}\n")
            if 'salmon' in diff:
                f.write(f"   Salmon: {diff['salmon']}\n")
            if 'cli' in diff:
                f.write(f"   CLI: {diff['cli']}\n")
            f.write("\n")
        
        if len(differences) > 50:
            f.write(f"... and {len(differences) - 50} more differences (truncated)\n")


def main():
    parser = argparse.ArgumentParser(description='Compare Salmon vs CLI trace outputs')
    parser.add_argument('--salmon-trace', required=True, help='Salmon trace file')
    parser.add_argument('--cli-trace', required=True, help='CLI trace file')
    parser.add_argument('--output', required=True, help='Output diff report file')
    
    args = parser.parse_args()
    
    print(f"Parsing Salmon trace: {args.salmon_trace}")
    salmon_reads = _parse_trace_file_detailed(args.salmon_trace)
    print(f"  Found {len(salmon_reads)} reads")
    
    print(f"Parsing CLI trace: {args.cli_trace}")
    cli_reads = _parse_trace_file_detailed(args.cli_trace)
    print(f"  Found {len(cli_reads)} reads")
    
    print("Comparing reads...")
    differences, component_stats = compare_reads_detailed(salmon_reads, cli_reads)
    
    print(f"Writing report to: {args.output}")
    write_detailed_report(differences, component_stats, args.output)
    
    if differences:
        print(f"\n⚠️  Found {len(differences)} differences")
        print(f"   First divergence: {differences[0]['qname']} ({differences[0]['type']})")
        print("\nComponent breakdown:")
        for stat, count in sorted(component_stats.items()):
            print(f"  {stat}: {count}")
        return 1
    else:
        print("\n✅ No differences found!")
        return 0


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python3
"""
Compare equivalence class files for Salmon parity testing.

Parses both Salmon and our EC files, canonicalizes by (labels, weights, count),
and reports mismatches with detailed diffs.
"""

import argparse
import gzip
import sys
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

@dataclass
class EquivalenceClass:
    """Represents a single equivalence class."""
    transcript_ids: List[int]      # Sorted transcript indices
    weights: List[float]           # Weights corresponding to transcript_ids
    count: float                   # Read count for this EC
    original_line: int = 0         # Line number in source file
    
    def canonical_key(self) -> Tuple:
        """Return a hashable canonical form for comparison."""
        # Sort by transcript ID, keeping weights aligned
        pairs = sorted(zip(self.transcript_ids, self.weights))
        sorted_ids = tuple(p[0] for p in pairs)
        sorted_weights = tuple(p[1] for p in pairs)
        return (sorted_ids, sorted_weights, self.count)
    
    def label_key(self) -> Tuple[int, ...]:
        """Return just the sorted transcript IDs."""
        return tuple(sorted(self.transcript_ids))


def parse_salmon_ec_file(filepath: str) -> Tuple[List[str], List[EquivalenceClass]]:
    """
    Parse Salmon eq_classes.txt format.
    
    Format:
      Line 1: num_transcripts
      Line 2: num_ecs
      Lines 3 to (3+num_transcripts-1): transcript names
      Remaining: <k> <idx1> ... <idxk> [<w1> ... <wk>] <count>
    
    Returns (transcript_names, list of ECs)
    """
    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'
    
    with open_func(filepath, mode) as f:
        lines = [line.strip() for line in f if line.strip()]
    
    num_transcripts = int(lines[0])
    num_ecs = int(lines[1])
    
    transcript_names = lines[2:2+num_transcripts]
    ec_lines = lines[2+num_transcripts:]
    
    ecs = []
    for line_num, line in enumerate(ec_lines):
        parts = line.split()
        k = int(parts[0])
        
        # Determine if weighted format (has 2k+2 parts) or unweighted (k+2 parts)
        if len(parts) == 2 * k + 2:
            # Weighted: <k> <idx1>...<idxk> <w1>...<wk> <count>
            transcript_ids = [int(parts[i]) for i in range(1, k+1)]
            weights = [float(parts[i]) for i in range(k+1, 2*k+1)]
            count = float(parts[-1])
        elif len(parts) == k + 2:
            # Unweighted: <k> <idx1>...<idxk> <count>
            transcript_ids = [int(parts[i]) for i in range(1, k+1)]
            weights = [1.0 / k] * k  # Uniform weights
            count = float(parts[-1])
        else:
            raise ValueError(f"Invalid EC format at line {line_num}: {line}")
        
        ecs.append(EquivalenceClass(
            transcript_ids=transcript_ids,
            weights=weights,
            count=count,
            original_line=line_num + 2 + num_transcripts + 1
        ))
    
    return transcript_names, ecs


def compare_ecs(
    salmon_ecs: List[EquivalenceClass],
    our_ecs: List[EquivalenceClass],
    tolerance: float = 1e-6
) -> Dict:
    """
    Compare two lists of equivalence classes.
    
    Returns a dict with comparison results.
    """
    results = {
        'salmon_ec_count': len(salmon_ecs),
        'our_ec_count': len(our_ecs),
        'matching_labels': 0,
        'matching_weights': 0,
        'matching_counts': 0,
        'fully_matching': 0,
        'label_mismatches': [],
        'weight_mismatches': [],
        'count_mismatches': [],
        'salmon_only': [],
        'ours_only': [],
    }
    
    # Build lookup by canonical label key
    salmon_by_label = {ec.label_key(): ec for ec in salmon_ecs}
    our_by_label = {ec.label_key(): ec for ec in our_ecs}
    
    salmon_labels = set(salmon_by_label.keys())
    our_labels = set(our_by_label.keys())
    
    # Find common and unique labels
    common_labels = salmon_labels & our_labels
    results['salmon_only'] = list(salmon_labels - our_labels)
    results['ours_only'] = list(our_labels - salmon_labels)
    
    # Compare ECs with matching labels
    for label in common_labels:
        salmon_ec = salmon_by_label[label]
        our_ec = our_by_label[label]
        
        results['matching_labels'] += 1
        
        # Compare weights (after sorting by transcript ID)
        s_pairs = sorted(zip(salmon_ec.transcript_ids, salmon_ec.weights))
        o_pairs = sorted(zip(our_ec.transcript_ids, our_ec.weights))
        
        weights_match = True
        for (s_id, s_w), (o_id, o_w) in zip(s_pairs, o_pairs):
            if abs(s_w - o_w) > tolerance:
                weights_match = False
                results['weight_mismatches'].append({
                    'label': label,
                    'transcript_id': s_id,
                    'salmon_weight': s_w,
                    'our_weight': o_w,
                    'diff': abs(s_w - o_w)
                })
        
        if weights_match:
            results['matching_weights'] += 1
        
        # Compare counts
        if abs(salmon_ec.count - our_ec.count) <= tolerance:
            results['matching_counts'] += 1
        else:
            results['count_mismatches'].append({
                'label': label,
                'salmon_count': salmon_ec.count,
                'our_count': our_ec.count,
                'diff': abs(salmon_ec.count - our_ec.count)
            })
        
        # Fully matching?
        if weights_match and abs(salmon_ec.count - our_ec.count) <= tolerance:
            results['fully_matching'] += 1
    
    return results


def print_report(results: Dict, salmon_names: List[str], out_file=None):
    """Print a human-readable comparison report."""
    out = out_file or sys.stdout
    
    print("=" * 60, file=out)
    print("SALMON PARITY TEST REPORT", file=out)
    print("=" * 60, file=out)
    print(file=out)
    
    print(f"Salmon ECs:     {results['salmon_ec_count']}", file=out)
    print(f"Our ECs:        {results['our_ec_count']}", file=out)
    print(file=out)
    
    print("MATCHING STATISTICS:", file=out)
    print(f"  Labels match:   {results['matching_labels']}/{results['salmon_ec_count']}", file=out)
    print(f"  Weights match:  {results['matching_weights']}/{results['matching_labels']}", file=out)
    print(f"  Counts match:   {results['matching_counts']}/{results['matching_labels']}", file=out)
    print(f"  Fully matching: {results['fully_matching']}/{results['matching_labels']}", file=out)
    print(file=out)
    
    if results['salmon_only']:
        print(f"SALMON-ONLY ECs ({len(results['salmon_only'])}):", file=out)
        for label in results['salmon_only'][:5]:
            print(f"  {label}", file=out)
        if len(results['salmon_only']) > 5:
            print(f"  ... and {len(results['salmon_only']) - 5} more", file=out)
        print(file=out)
    
    if results['ours_only']:
        print(f"OUR-ONLY ECs ({len(results['ours_only'])}):", file=out)
        for label in results['ours_only'][:5]:
            print(f"  {label}", file=out)
        if len(results['ours_only']) > 5:
            print(f"  ... and {len(results['ours_only']) - 5} more", file=out)
        print(file=out)
    
    if results['weight_mismatches']:
        print(f"WEIGHT MISMATCHES ({len(results['weight_mismatches'])}):", file=out)
        for m in results['weight_mismatches'][:10]:
            txp_name = salmon_names[m['transcript_id']] if m['transcript_id'] < len(salmon_names) else f"txp_{m['transcript_id']}"
            print(f"  EC {m['label'][:3]}...: {txp_name} salmon={m['salmon_weight']:.6f} ours={m['our_weight']:.6f} diff={m['diff']:.2e}", file=out)
        if len(results['weight_mismatches']) > 10:
            print(f"  ... and {len(results['weight_mismatches']) - 10} more", file=out)
        print(file=out)
    
    if results['count_mismatches']:
        print(f"COUNT MISMATCHES ({len(results['count_mismatches'])}):", file=out)
        for m in results['count_mismatches'][:10]:
            print(f"  EC {m['label'][:3]}...: salmon={m['salmon_count']:.1f} ours={m['our_count']:.1f} diff={m['diff']:.2f}", file=out)
        if len(results['count_mismatches']) > 10:
            print(f"  ... and {len(results['count_mismatches']) - 10} more", file=out)
        print(file=out)
    
    # Summary
    parity_pct = 100.0 * results['fully_matching'] / max(1, results['matching_labels'])
    print("=" * 60, file=out)
    print(f"PARITY: {parity_pct:.1f}% ({results['fully_matching']}/{results['matching_labels']} ECs)", file=out)
    
    if parity_pct >= 99.0:
        print("STATUS: PASS (>=99% parity)", file=out)
    elif parity_pct >= 95.0:
        print("STATUS: WARN (>=95% but <99% parity)", file=out)
    else:
        print("STATUS: FAIL (<95% parity)", file=out)
    print("=" * 60, file=out)


def main():
    parser = argparse.ArgumentParser(description='Compare EC files for Salmon parity')
    parser.add_argument('--salmon', required=True, help='Salmon eq_classes.txt')
    parser.add_argument('--ours', required=True, help='Our eq_classes.txt')
    parser.add_argument('--tolerance', type=float, default=1e-6, help='Weight comparison tolerance')
    parser.add_argument('--report', help='Output report file (default: stdout)')
    args = parser.parse_args()
    
    # Parse both files
    salmon_names, salmon_ecs = parse_salmon_ec_file(args.salmon)
    our_names, our_ecs = parse_salmon_ec_file(args.ours)
    
    # Compare
    results = compare_ecs(salmon_ecs, our_ecs, args.tolerance)
    
    # Print report
    if args.report:
        with open(args.report, 'w') as f:
            print_report(results, salmon_names, f)
        print(f"Report written to: {args.report}")
    else:
        print_report(results, salmon_names)
    
    # Exit with error if parity < 95%
    parity_pct = 100.0 * results['fully_matching'] / max(1, results['matching_labels'])
    sys.exit(0 if parity_pct >= 95.0 else 1)


if __name__ == '__main__':
    main()

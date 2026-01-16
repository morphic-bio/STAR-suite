#!/usr/bin/env python3
"""
Verify EC filtering unit test results.
"""

import json
import sys
import math

def verify_hard_filter(result, expected):
    """Verify hardFilter test results."""
    kept_tids = sorted([a['tid'] for a in result['kept_alignments']])
    expected_tids = sorted(expected['expected_tids'])
    
    if kept_tids != expected_tids:
        print(f"FAIL: hardFilter test - Expected {expected_tids}, got {kept_tids}")
        return False
    print("PASS: hardFilter test")
    return True

def verify_min_aln_prob(result, expected):
    """Verify minAlnProb test results."""
    kept_tids = sorted([a['tid'] for a in result['kept_alignments']])
    expected_tids = sorted(expected['expected_tids'])
    
    if kept_tids != expected_tids:
        print(f"FAIL: minAlnProb test - Expected {expected_tids}, got {kept_tids}")
        return False
    
    # Verify estAlnProb calculation
    best_score = expected.get('best_score', 100)
    for aln in result['kept_alignments']:
        v = best_score - aln['score']
        expected_prob = math.exp(-expected['params']['score_exp'] * v)
        if abs(aln['est_aln_prob'] - expected_prob) > 1e-6:
            print(f"FAIL: minAlnProb test - estAlnProb mismatch for tid={aln['tid']}")
            return False
    
    print("PASS: minAlnProb test")
    return True

def verify_per_transcript_best_hit(result, expected):
    """Verify per-transcript best hit test results."""
    kept_tids = sorted([a['tid'] for a in result['kept_alignments']])
    expected_tids = sorted(expected['expected_tids'])
    
    if kept_tids != expected_tids:
        print(f"FAIL: per-transcript best hit test - Expected {expected_tids}, got {kept_tids}")
        return False
    
    # Verify only one alignment per transcript
    tids_seen = {}
    for aln in result['kept_alignments']:
        tid = aln['tid']
        if tid in tids_seen:
            print(f"FAIL: per-transcript best hit test - Multiple alignments for tid={tid}")
            return False
        tids_seen[tid] = True
    
    print("PASS: per-transcript best hit test")
    return True

def verify_decoy_threshold(result, expected):
    """Verify decoy threshold test results."""
    kept_tids = sorted([a['tid'] for a in result['kept_alignments']])
    expected_tids = sorted(expected['expected_tids'])
    
    if kept_tids != expected_tids:
        print(f"FAIL: decoy threshold test - Expected {expected_tids}, got {kept_tids}")
        return False
    print("PASS: decoy threshold test")
    return True

def verify_range_factorization(result, expected):
    """Verify range factorization test results."""
    ec = result['ec']
    range_count = result['range_count']
    expected_range_count = expected['expected_range_count']
    
    if range_count != expected_range_count:
        print(f"FAIL: range factorization test - Expected range_count={expected_range_count}, got {range_count}")
        return False
    
    # Check bin IDs (should be appended to transcript_ids)
    if len(ec['transcript_ids']) < len(expected['ec']['transcript_ids']) + len(expected['expected_bin_ids']):
        print(f"FAIL: range factorization test - Missing bin IDs")
        return False
    
    # Extract bin IDs (last N values)
    original_count = len(expected['ec']['transcript_ids'])
    bin_ids = ec['transcript_ids'][original_count:]
    
    if bin_ids != expected['expected_bin_ids']:
        print(f"FAIL: range factorization test - Expected bin_ids={expected['expected_bin_ids']}, got {bin_ids}")
        return False
    
    print("PASS: range factorization test")
    return True

def verify_rank_eq_classes(result, expected):
    """Verify rankEqClasses test results."""
    ec = result['ec']
    sorted_ids = ec['transcript_ids']
    sorted_probs = ec['aux_probs']
    
    expected_ids = expected['expected_sorted_ids']
    expected_probs = expected['expected_sorted_probs']
    
    if sorted_ids != expected_ids:
        print(f"FAIL: rankEqClasses test - Expected ids={expected_ids}, got {sorted_ids}")
        return False
    
    for i, (got, exp) in enumerate(zip(sorted_probs, expected_probs)):
        if abs(got - exp) > 1e-6:
            print(f"FAIL: rankEqClasses test - Prob mismatch at index {i}: expected={exp}, got={got}")
            return False
    
    print("PASS: rankEqClasses test")
    return True

def main():
    """Run verification on test results."""
    if len(sys.argv) < 2:
        print("Usage: verify_filtering.py <results.json>")
        sys.exit(1)
    
    with open(sys.argv[1], 'r') as f:
        results = json.load(f)
    
    with open('fixtures.json', 'r') as f:
        fixtures = json.load(f)
    
    all_passed = True
    
    for test_result in results['test_results']:
        test_id = test_result['test_id']
        
        # Find corresponding fixture
        fixture = None
        for f in fixtures['tests']:
            if f['test_id'] == test_id:
                fixture = f
                break
        
        if not fixture:
            print(f"ERROR: No fixture found for test_id={test_id}")
            all_passed = False
            continue
        
        # Verify based on test type
        if test_id == 'hard_filter':
            all_passed &= verify_hard_filter(test_result, fixture)
        elif test_id == 'min_aln_prob':
            all_passed &= verify_min_aln_prob(test_result, fixture)
        elif test_id == 'per_transcript_best_hit':
            all_passed &= verify_per_transcript_best_hit(test_result, fixture)
        elif test_id == 'decoy_threshold':
            all_passed &= verify_decoy_threshold(test_result, fixture)
        elif test_id == 'range_factorization':
            all_passed &= verify_range_factorization(test_result, fixture)
        elif test_id == 'rank_eq_classes':
            all_passed &= verify_rank_eq_classes(test_result, fixture)
        else:
            print(f"WARN: Unknown test_id={test_id}")
    
    if all_passed:
        print("\nAll tests PASSED")
        sys.exit(0)
    else:
        print("\nSome tests FAILED")
        sys.exit(1)

if __name__ == '__main__':
    main()

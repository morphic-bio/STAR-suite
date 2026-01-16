#!/usr/bin/env python3
"""
Generate test fixtures for EC filtering unit tests.
Creates Salmon-specific test scenarios.
"""

import json
import math

def generate_hard_filter_test():
    """Test 1: hardFilter Behavior"""
    return {
        "test_id": "hard_filter",
        "description": "With hardFilter=true, only best score kept",
        "alignments": [
            {"tid": 0, "score": 100, "log_frag_prob": 0.0, "log_compat_prob": 0.0, "is_decoy": False},
            {"tid": 1, "score": 95, "log_frag_prob": 0.0, "log_compat_prob": 0.0, "is_decoy": False},
            {"tid": 2, "score": 90, "log_frag_prob": 0.0, "log_compat_prob": 0.0, "is_decoy": False}
        ],
        "params": {
            "hard_filter": True,
            "score_exp": 1.0,
            "min_aln_prob": 1e-5,
            "decoy_threshold": 1.0
        },
        "expected_tids": [0]  # Only best score kept
    }

def generate_min_aln_prob_test():
    """Test 2: minAlnProb Gating"""
    return {
        "test_id": "min_aln_prob",
        "description": "Filter alignments below minAlnProb threshold",
        "alignments": [
            {"tid": 0, "score": 100, "log_frag_prob": 0.0, "log_compat_prob": 0.0, "is_decoy": False},
            {"tid": 1, "score": 95, "log_frag_prob": 0.0, "log_compat_prob": 0.0, "is_decoy": False},
            {"tid": 2, "score": 50, "log_frag_prob": 0.0, "log_compat_prob": 0.0, "is_decoy": False}
        ],
        "params": {
            "hard_filter": False,
            "score_exp": 1.0,
            "min_aln_prob": 1e-5,
            "decoy_threshold": 1.0
        },
        "expected_tids": [0, 1],  # tid=2 filtered (estAlnProb = exp(-50) too low)
        "best_score": 100
    }

def generate_per_transcript_best_hit_test():
    """Test 3: Per-Transcript Best Hit"""
    return {
        "test_id": "per_transcript_best_hit",
        "description": "Keep only best alignment per transcript",
        "alignments": [
            {"tid": 0, "score": 95, "log_frag_prob": 0.0, "log_compat_prob": 0.0, "is_decoy": False},
            {"tid": 0, "score": 100, "log_frag_prob": 0.0, "log_compat_prob": 0.0, "is_decoy": False},
            {"tid": 1, "score": 90, "log_frag_prob": 0.0, "log_compat_prob": 0.0, "is_decoy": False}
        ],
        "params": {
            "hard_filter": False,
            "score_exp": 1.0,
            "min_aln_prob": 1e-5,
            "decoy_threshold": 1.0
        },
        "expected_tids": [0, 1],  # Only best hit per transcript kept
        "best_score": 100
    }

def generate_decoy_threshold_test():
    """Test 4: Decoy Threshold"""
    return {
        "test_id": "decoy_threshold",
        "description": "Filter alignments below decoy threshold",
        "alignments": [
            {"tid": 0, "score": 100, "log_frag_prob": 0.0, "log_compat_prob": 0.0, "is_decoy": False},
            {"tid": 1, "score": 75, "log_frag_prob": 0.0, "log_compat_prob": 0.0, "is_decoy": False},
            {"tid": 2, "score": 80, "log_frag_prob": 0.0, "log_compat_prob": 0.0, "is_decoy": True}
        ],
        "params": {
            "hard_filter": False,
            "score_exp": 1.0,
            "min_aln_prob": 1e-5,
            "decoy_threshold": 1.0
        },
        "expected_tids": [0],  # tid=1 filtered (75 < 1.0 * 80), tid=2 is decoy
        "best_decoy_score": 80,
        "best_score": 100
    }

def generate_range_factorization_test():
    """Test 5: Range Factorization Bins"""
    return {
        "test_id": "range_factorization",
        "description": "Append range bin IDs to EC labels",
        "ec": {
            "transcript_ids": [0, 1, 2, 3],
            "aux_probs": [0.5, 0.3, 0.15, 0.05]
        },
        "params": {
            "range_factorization_bins": 4
        },
        "expected_range_count": 6,  # sqrt(4) + 4 = 6
        "expected_bin_ids": [3, 1, 0, 0]  # floor(auxProb * rangeCount)
    }

def generate_rank_eq_classes_test():
    """Test 6: rankEqClasses"""
    return {
        "test_id": "rank_eq_classes",
        "description": "Sort transcripts by conditional probability",
        "ec": {
            "transcript_ids": [0, 1, 2],
            "aux_probs": [0.5, 0.1, 0.4]
        },
        "params": {
            "use_rank_eq_classes": True
        },
        "expected_sorted_ids": [1, 2, 0],  # Sorted by aux_probs ascending
        "expected_sorted_probs": [0.1, 0.4, 0.5]
    }

def main():
    """Generate all test fixtures."""
    fixtures = {
        "tests": [
            generate_hard_filter_test(),
            generate_min_aln_prob_test(),
            generate_per_transcript_best_hit_test(),
            generate_decoy_threshold_test(),
            generate_range_factorization_test(),
            generate_rank_eq_classes_test()
        ]
    }
    
    print(json.dumps(fixtures, indent=2))

if __name__ == '__main__':
    main()

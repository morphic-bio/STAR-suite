#!/usr/bin/env python3
"""
Generate test fixtures using cutadapt 3.2 for STAR-Flex compatibility testing.

This script must be run inside a Docker container with cutadapt 3.2 installed:
    docker run --rm -v $(pwd):/work -w /work python:3.9-slim bash -c \
        "pip install -q cutadapt==3.2 && python3 tests/fixtures/generate_cutadapt32_fixtures.py"
"""

import json
import sys

# Check cutadapt version
try:
    import cutadapt
    from cutadapt._align import Aligner
    from cutadapt import align
    version = cutadapt.__version__
    if not version.startswith("3.2"):
        print(f"WARNING: Expected cutadapt 3.2, got {version}", file=sys.stderr)
except ImportError as e:
    print(f"ERROR: Could not import cutadapt: {e}", file=sys.stderr)
    sys.exit(1)

# Flags for 3' adapter (BackAdapter)
# START_WITHIN_SEQ2 | STOP_WITHIN_SEQ2 | STOP_WITHIN_SEQ1
BACK_FLAGS = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ1

# Standard Illumina adapter
ILLUMINA_ADAPTER = "AGATCGGAAGAGC"

def generate_test_cases():
    """Generate a variety of test cases for adapter matching."""
    test_cases = []
    
    # Test case 1: Exact full adapter at end
    test_cases.append({
        "name": "exact_full_adapter_at_end",
        "read": "ACGTACGTACGTACGT" + ILLUMINA_ADAPTER,
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.1,
    })
    
    # Test case 2: Exact partial adapter (3bp prefix)
    test_cases.append({
        "name": "exact_3bp_prefix",
        "read": "ACGTACGTACGTACGTACGTACGT" + "AGA",
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.1,
    })
    
    # Test case 3: Exact partial adapter (1bp prefix)  
    test_cases.append({
        "name": "exact_1bp_prefix",
        "read": "ACGTACGTACGTACGTACGTACGT" + "A",
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.1,
    })
    
    # Test case 4: Adapter with 1 mismatch
    test_cases.append({
        "name": "adapter_1_mismatch",
        "read": "ACGTACGTACGTACGT" + "AGATCGGTAGAGC",  # A->T at position 7
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.1,
    })
    
    # Test case 5: No adapter
    test_cases.append({
        "name": "no_adapter",
        "read": "ACGTACGTACGTACGTACGTACGTACGT",
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.1,
    })
    
    # Test case 6: Adapter in middle (should find it)
    test_cases.append({
        "name": "adapter_in_middle",
        "read": "ACGTACGT" + ILLUMINA_ADAPTER + "ACGTACGT",
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.1,
    })
    
    # Test case 7: Multiple possible matches - critical test for 3.2 vs 5.x difference
    # Read has adapter-like sequence at multiple positions
    test_cases.append({
        "name": "multiple_possible_matches",
        "read": "AGATCGGAAGAGCACGTACGTAAGAGC",  # Full adapter early, partial later
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.1,
    })
    
    # Test case 8: Short overlap (2bp)
    test_cases.append({
        "name": "short_overlap_2bp",
        "read": "ACGTACGTACGTACGTACGTACGTAG",  # AG at end matches adapter prefix
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.1,
    })
    
    # Test case 9: Edge case - entire read is adapter
    test_cases.append({
        "name": "entire_read_is_adapter",
        "read": ILLUMINA_ADAPTER,
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.1,
    })
    
    # Test case 10: Adapter with insertion
    test_cases.append({
        "name": "adapter_with_insertion",
        "read": "ACGTACGTACGTACGT" + "AGAXTCGGAAGAGC",  # X inserted
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.15,
    })
    
    # Test case 11: Adapter with deletion
    test_cases.append({
        "name": "adapter_with_deletion",
        "read": "ACGTACGTACGTACGT" + "AGATCGAAGAGC",  # G deleted
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.15,
    })
    
    # Test case 12: min_overlap=3 filtering
    test_cases.append({
        "name": "min_overlap_3_filter",
        "read": "ACGTACGTACGTACGTACGTACGTAG",  # AG at end (2bp)
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 3,
        "max_error_rate": 0.1,
    })
    
    # Test case 13: Ambiguous adapter position - key differentiator
    # This tests the matches vs score comparison
    test_cases.append({
        "name": "ambiguous_position_test",
        "read": "AGATCGGACGTAAGAGCAGATCGGAAGAGC",  # Two similar matches
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.1,
    })
    
    # Test case 14: Very short read with partial adapter
    test_cases.append({
        "name": "short_read_partial_adapter",
        "read": "ACGTAGA",  # 4bp + AGA
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.1,
    })
    
    # Test case 15: Lowercase handling
    test_cases.append({
        "name": "lowercase_adapter",
        "read": "ACGTACGTACGTACGT" + "agatcggaagagc",
        "adapter": ILLUMINA_ADAPTER,
        "min_overlap": 1,
        "max_error_rate": 0.1,
    })
    
    return test_cases


def run_cutadapt_alignment(read, adapter, min_overlap, max_error_rate):
    """
    Run cutadapt 3.2 aligner and return the result.
    Returns (ref_start, ref_stop, query_start, query_stop, matches, errors) or None
    """
    aligner = Aligner(
        adapter.upper(),
        max_error_rate,
        flags=BACK_FLAGS,
        wildcard_ref=False,
        wildcard_query=False,
        indel_cost=1,
        min_overlap=min_overlap,
    )
    
    result = aligner.locate(read.upper())
    return result


def main():
    print(f"Generating fixtures with cutadapt {cutadapt.__version__}")
    
    test_cases = generate_test_cases()
    fixtures = []
    
    for tc in test_cases:
        result = run_cutadapt_alignment(
            tc["read"],
            tc["adapter"],
            tc["min_overlap"],
            tc["max_error_rate"],
        )
        
        fixture = {
            "name": tc["name"],
            "read": tc["read"],
            "adapter": tc["adapter"],
            "min_overlap": tc["min_overlap"],
            "max_error_rate": tc["max_error_rate"],
            "cutadapt_version": cutadapt.__version__,
        }
        
        if result is None:
            fixture["expected_trim_pos"] = len(tc["read"])  # No trim
            fixture["match"] = None
        else:
            ref_start, ref_stop, query_start, query_stop, matches, errors = result
            fixture["expected_trim_pos"] = query_start
            fixture["match"] = {
                "ref_start": ref_start,
                "ref_stop": ref_stop,
                "query_start": query_start,
                "query_stop": query_stop,
                "matches": matches,
                "errors": errors,
            }
        
        fixtures.append(fixture)
        
        # Print summary
        if fixture["match"]:
            m = fixture["match"]
            print(f"  {tc['name']}: trim at {fixture['expected_trim_pos']} "
                  f"(matches={m['matches']}, errors={m['errors']})")
        else:
            print(f"  {tc['name']}: no match (no trim)")
    
    # Write fixtures to JSON
    output_path = "tests/fixtures/cutadapt32_adapter_fixtures.json"
    with open(output_path, "w") as f:
        json.dump(fixtures, f, indent=2)
    
    print(f"\nWrote {len(fixtures)} fixtures to {output_path}")


if __name__ == "__main__":
    main()

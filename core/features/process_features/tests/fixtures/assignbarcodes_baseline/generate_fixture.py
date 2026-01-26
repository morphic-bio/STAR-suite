#!/usr/bin/env python3
"""
Generate synthetic FASTQ fixture for assignBarcodes regression testing.
Creates R1 (barcode+UMI) and R2 (feature sequence) FASTQs with known data.
"""

import os

# Input directory
INPUT_DIR = os.path.dirname(__file__) + "/input"

# Feature sequences from features.csv
FEATURES = {
    "FeatureA": "ATCGATCGATCGATCG",
    "FeatureB": "GCTAGCTAGCTAGCTA",
    "FeatureC": "TTAATTAATTAATTAA",
    "FeatureD": "CCGGCCGGCCGGCCGG",
    "FeatureE": "AACCGGTTAACCGGTT",
}

# Barcodes from whitelist (first 5)
BARCODES = [
    "AAACCCAAGAAACCAT",
    "AAACCCAAGAAACCCA",
    "AAACCCAAGAAACCCT",
    "AAACCCAAGAAACCGG",
    "AAACCCAAGAAACCTA",
]

# Test data: (barcode_idx, feature_name, umi_suffix) tuples
# umi_suffix determines uniqueness for deduplication
TEST_READS = [
    # Barcode 0 gets: FeatureA x3 (2 unique UMIs), FeatureB x2 (1 unique)
    (0, "FeatureA", "AAAAAAAAAAAA"),
    (0, "FeatureA", "AAAAAAAAAAAA"),  # duplicate UMI
    (0, "FeatureA", "AAAAAAAAAAAB"),
    (0, "FeatureB", "BBBBBBBBBBBB"),
    (0, "FeatureB", "BBBBBBBBBBBB"),  # duplicate
    
    # Barcode 1 gets: FeatureC x4 (3 unique)
    (1, "FeatureC", "CCCCCCCCCCCC"),
    (1, "FeatureC", "CCCCCCCCCCCC"),  # duplicate
    (1, "FeatureC", "CCCCCCCCCCD1"),
    (1, "FeatureC", "CCCCCCCCCCD2"),
    
    # Barcode 2 gets: FeatureD x2, FeatureE x2 (all unique)
    (2, "FeatureD", "DDDDDDDDDDDD"),
    (2, "FeatureD", "DDDDDDDDDDD2"),
    (2, "FeatureE", "EEEEEEEEEEEE"),
    (2, "FeatureE", "EEEEEEEEEEE2"),
    
    # Barcode 3 gets: FeatureA x1
    (3, "FeatureA", "AAAAAAAAAA03"),
    
    # Barcode 4 gets: FeatureB x1, FeatureC x1, FeatureD x1
    (4, "FeatureB", "BBBBBBBBB004"),
    (4, "FeatureC", "CCCCCCCCC004"),
    (4, "FeatureD", "DDDDDDDDD004"),
]

def generate_quality(length):
    """Generate quality string of given length (all high quality)."""
    return "I" * length

def main():
    r1_path = os.path.join(INPUT_DIR, "test_R1.fastq")
    r2_path = os.path.join(INPUT_DIR, "test_R2.fastq")
    
    with open(r1_path, "w") as r1_file, open(r2_path, "w") as r2_file:
        for i, (bc_idx, feature_name, umi) in enumerate(TEST_READS):
            barcode = BARCODES[bc_idx]
            feature_seq = FEATURES[feature_name]
            
            # R1: barcode (16bp) + UMI (12bp) = 28bp
            r1_seq = barcode + umi
            r1_qual = generate_quality(len(r1_seq))
            
            # R2: feature sequence (16bp) + poly-A padding (no N's to avoid filter)
            r2_seq = feature_seq + "AAAAAAAAAAAA"  # 28bp total
            r2_qual = generate_quality(len(r2_seq))
            
            read_name = f"read_{i:04d}"
            
            # Write R1
            r1_file.write(f"@{read_name}\n")
            r1_file.write(f"{r1_seq}\n")
            r1_file.write("+\n")
            r1_file.write(f"{r1_qual}\n")
            
            # Write R2
            r2_file.write(f"@{read_name}\n")
            r2_file.write(f"{r2_seq}\n")
            r2_file.write("+\n")
            r2_file.write(f"{r2_qual}\n")
    
    print(f"Generated {len(TEST_READS)} reads")
    print(f"  R1: {r1_path}")
    print(f"  R2: {r2_path}")
    
    # Print expected counts for verification
    print("\nExpected deduped counts per barcode:")
    from collections import defaultdict
    deduped = defaultdict(lambda: defaultdict(set))
    for bc_idx, feature_name, umi in TEST_READS:
        deduped[bc_idx][feature_name].add(umi)
    
    for bc_idx in sorted(deduped.keys()):
        bc = BARCODES[bc_idx]
        counts = {f: len(umis) for f, umis in deduped[bc_idx].items()}
        total = sum(counts.values())
        print(f"  {bc}: {counts} (total: {total})")

if __name__ == "__main__":
    main()

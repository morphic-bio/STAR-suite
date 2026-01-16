# GC Bias Correction - Implementation Handoff

## Status: ✅ Phase 1 Complete - Expected GC Computation Tool

### What Has Been Implemented

#### 1. Standalone C Tool: `compute_expected_gc`

**Location**: `/mnt/pikachu/STAR-Flex/tools/compute_expected_gc/`

**Features**:
- Pure C implementation with zero dependencies (only libc and libm)
- Two input modes:
  - `--transcriptome FILE`: Direct transcriptome FASTA input
  - `--genome FILE --gtf FILE`: Extract transcriptome from genome using GTF annotations
- Computes expected GC distribution (101 bins, 0-100% GC)
- Default fragment length distribution (normal-like, mean=200, sd=80)
- Optional custom FLD via `--fld FILE`

**Output Format**: `expected_gc.tsv`
```
0	0.0001234567
1	0.0002345678
...
100	0.0003456789
```
101 rows total, tab-separated: `gc_percent<TAB>probability`

**Build**: `cd tools/compute_expected_gc && make`

**Test**: 
```bash
./compute_expected_gc --transcriptome test/test_transcriptome.fa --output test.tsv
```

#### 2. STAR Integration

**Modified Files**:
- `source/ParametersGenome.h`: Added `string transcriptomeFasta;`
- `source/Parameters.cpp`: Added parameter parsing for `transcriptomeFasta`
- `source/parametersDefault`: Added default value `-` for `transcriptomeFasta`
- `source/Genome_genomeGenerate.cpp`: Added automatic call to `compute_expected_gc` during genome generation

**Integration Logic**:
- After genome files are written, checks if GTF or transcriptome FASTA is provided
- If `--transcriptomeFasta` provided: uses it directly
- If `--sjdbGTFfile` provided: extracts transcriptome from genome+GTF
- Calls `compute_expected_gc` and stores `expected_gc.tsv` in genome directory
- Tool path detection with multiple fallback locations

**Usage**:
```bash
# Option 1: Provide transcriptome FASTA
STAR --runMode genomeGenerate --genomeDir genome_index \
     --genomeFastaFiles genome.fa --transcriptomeFasta transcripts.fa

# Option 2: Extract from genome+GTF (automatic)
STAR --runMode genomeGenerate --genomeDir genome_index \
     --genomeFastaFiles genome.fa --sjdbGTFfile annotations.gtf
```

**Output**: `genome_index/expected_gc.tsv` is automatically created

## What Remains: Phase 2 - GC Bias Application

### Goal
Use the pre-computed `expected_gc.tsv` during quantification to:
1. Collect observed GC distribution from aligned fragments
2. Compute bias weights: `bias[gc] = observed[gc] / expected[gc]`
3. Apply bias weights to effective length calculations
4. Use corrected effective lengths in EC filter and VB optimization

### Implementation Steps

#### Step 1: Load Expected GC Distribution

**Location**: Wherever transcriptome/EC building is initialized

**Code**:
```cpp
#include <fstream>
#include <vector>

std::vector<double> load_expected_gc(const std::string& genome_dir) {
    std::vector<double> expected_gc(101, 0.0);
    std::string filename = genome_dir + "/expected_gc.tsv";
    
    std::ifstream f(filename);
    if (!f.is_open()) {
        // File doesn't exist - GC bias correction not available
        return expected_gc;  // Return zeros
    }
    
    int gc_pct;
    double prob;
    while (f >> gc_pct >> prob) {
        if (gc_pct >= 0 && gc_pct < 101) {
            expected_gc[gc_pct] = prob;
        }
    }
    return expected_gc;
}
```

#### Step 2: Collect Observed GC During Alignment Processing

**Location**: During EC building or alignment processing loop

**Code**:
```cpp
std::vector<double> observed_gc(101, 0.0);

// For each aligned fragment
for (auto& alignment : alignments) {
    // Extract fragment sequence from BAM
    std::string fragment_seq = extract_fragment_sequence(alignment);
    
    // Compute GC percentage
    int gc_count = 0;
    for (char c : fragment_seq) {
        if (c == 'G' || c == 'C' || c == 'g' || c == 'c') gc_count++;
    }
    int gc_pct = (int)round(100.0 * gc_count / fragment_seq.length());
    
    // Accumulate with fragment weight
    double weight = alignment.weight;  // From error model or uniform
    observed_gc[gc_pct] += weight;
}

// Normalize observed GC
double total = 0.0;
for (double w : observed_gc) total += w;
if (total > 0) {
    for (double& w : observed_gc) w /= total;
}
```

#### Step 3: Compute Bias Weights

**Code**:
```cpp
std::vector<double> compute_gc_bias_weights(
    const std::vector<double>& observed,
    const std::vector<double>& expected) {
    
    std::vector<double> bias(101, 1.0);
    double prior = 1e-5;  // Small prior to avoid division by zero
    
    for (int i = 0; i < 101; i++) {
        if (expected[i] > 0) {
            bias[i] = observed[i] / (prior + expected[i]);
            // Clamp to reasonable range (Salmon uses 0.01 to 100)
            if (bias[i] > 100.0) bias[i] = 100.0;
            if (bias[i] < 0.01) bias[i] = 0.01;
        }
    }
    return bias;
}
```

#### Step 4: Apply to Effective Lengths

**Location**: Effective length computation (similar to Salmon's `SalmonUtils.cpp:3080-3125`)

**Code**:
```cpp
double compute_effective_length_with_gc_bias(
    const Transcript& txp,
    const FragmentLengthDistribution& fld,
    const std::vector<double>& gc_bias) {
    
    double eff_length = 0.0;
    
    // For each fragment length
    for (int fl = fld.min(); fl <= fld.max(); fl++) {
        double fl_weight = exp(fld.pmf(fl));
        
        // For each start position
        for (int start = 0; start < txp.length - fl; start++) {
            int end = start + fl;
            
            // Compute GC percentage for this fragment
            int gc_count = 0;
            for (int pos = start; pos < end; pos++) {
                char c = toupper(txp.sequence[pos]);
                if (c == 'G' || c == 'C') gc_count++;
            }
            int gc_pct = (int)round(100.0 * gc_count / fl);
            
            // Apply GC bias weight
            double gc_weight = gc_bias[gc_pct];
            
            eff_length += fl_weight * gc_weight;
        }
    }
    
    return eff_length;
}
```

### Files to Modify

1. **EC Builder** (`source/libem/ec_builder.cpp` or similar):
   - Load `expected_gc.tsv` during initialization
   - Collect observed GC during read processing
   - Compute bias weights
   - Store bias weights for use in effective length calculations

2. **Effective Length Computation** (wherever transcript effective lengths are computed):
   - Add GC bias weight application
   - Use GC bias weights when computing effective lengths

3. **VB Optimization** (if effective lengths affect VB):
   - Ensure GC-corrected effective lengths are used

### Testing Strategy

1. **Unit Test Expected GC Loading**:
   ```cpp
   auto expected = load_expected_gc("/path/to/genome");
   assert(expected.size() == 101);
   assert(expected[50] > 0);  // Some probability at 50% GC
   ```

2. **Integration Test**:
   - Run STAR genome generation with GTF
   - Verify `expected_gc.tsv` is created
   - Run quantification
   - Verify observed GC is collected
   - Verify bias weights are computed
   - Compare effective lengths with/without GC bias

3. **Parity Test with Salmon**:
   - Run Salmon with `--gcBias` on same data
   - Compare GC bias weights
   - Compare effective lengths
   - Compare final quantification results

### Reference: Salmon's Implementation

**Key Files**:
- `/mnt/pikachu/salmon/src/SalmonUtils.cpp:3091-3110` - Effective length with GC bias
- `/mnt/pikachu/salmon/src/SalmonUtils.cpp:3957-3964` - Bias weight computation
- `/mnt/pikachu/salmon/include/GCFragModel.hpp` - GC bias model class

**Salmon's Approach**:
- Uses `GCFragModel` class (Eigen matrix) for conditional GC bias (by context)
- We can simplify to unconditional GC bias (just fragment GC%)
- Salmon computes expected GC during genome indexing (similar to our approach)
- Applies bias during effective length computation

### Questions to Resolve

1. **When to collect observed GC?**
   - During initial alignment processing?
   - During EC building?
   - During VB optimization?

2. **Fragment sequence extraction**:
   - Extract from BAM sequence field?
   - Reconstruct from CIGAR + reference?
   - Use transcript sequence directly?

3. **Conditional vs Unconditional GC bias**:
   - Salmon uses conditional GC bias (by 5'/3' context)
   - Simpler: unconditional (just fragment GC%)
   - Decision: Start simple, add conditional later if needed?

4. **Integration point**:
   - Modify existing effective length computation?
   - Add new GC-aware effective length function?
   - Where is effective length currently computed?

### Next Steps

1. **Locate effective length computation** in STAR-Flex codebase
2. **Identify where to collect observed GC** (alignment processing loop)
3. **Implement GC bias weight computation**
4. **Modify effective length computation** to use GC bias
5. **Test with small dataset**
6. **Compare with Salmon** for parity

### Files Created/Modified

**New Files**:
- `tools/compute_expected_gc/compute_expected_gc.c` - Main tool
- `tools/compute_expected_gc/Makefile` - Build file
- `tools/compute_expected_gc/README.md` - Documentation
- `tools/compute_expected_gc/test/test_transcriptome.fa` - Test data

**Modified Files**:
- `source/ParametersGenome.h` - Added `transcriptomeFasta` parameter
- `source/Parameters.cpp` - Added parameter parsing
- `source/parametersDefault` - Added default value
- `source/Genome_genomeGenerate.cpp` - Added tool call integration

### Verification

**Test the tool**:
```bash
cd /mnt/pikachu/STAR-Flex/tools/compute_expected_gc
make
./compute_expected_gc --transcriptome test/test_transcriptome.fa --output test.tsv
head test.tsv  # Should show 101 rows
```

**Parity test with Salmon**:
```bash
# Run parity test script
./test_gc_parity.sh

# Or manually:
# 1. Run Salmon with --gcBias
salmon quant -t transcriptome.fa -l A -a alignments.bam --gcBias -o salmon_out

# 2. Run our tool
./compute_expected_gc --transcriptome transcriptome.fa --output our_expected_gc.tsv

# 3. Compare
python3 compare_with_salmon.py salmon_out/aux_info/exp_gc.gz our_expected_gc.tsv
```

**Note**: Parity check has NOT been performed yet. The test scripts are ready but need to be run with actual data.

**Test STAR integration**:
```bash
# Generate genome index with GTF
STAR --runMode genomeGenerate \
     --genomeDir test_genome \
     --genomeFastaFiles /path/to/genome.fa \
     --sjdbGTFfile /path/to/annotations.gtf

# Verify expected_gc.tsv was created
ls test_genome/expected_gc.tsv
head test_genome/expected_gc.tsv
```

### Notes

- The tool is **pure C** with zero dependencies - very portable
- Tool path detection in STAR integration tries multiple locations
- If tool not found, STAR continues without GC bias (warning only)
- Expected GC computation is **reference-specific** - computed once per genome
- Observed GC is **sample-specific** - computed per sample during quantification

### Contact

For questions about the implementation, refer to:
- This document
- Salmon source code for reference implementation
- `tools/compute_expected_gc/README.md` for tool usage

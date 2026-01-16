# Debug Instrumentation Status

## ✅ Completed: em_quant Instrumentation

### Implementation
- ✅ Enhanced debug tracing in `source/libem/vb_engine.cpp`
- ✅ Per-iteration logging: alpha, logNorm, expTheta, expected_count
- ✅ Per-EC logging: ec_id, denom, expTheta, aux, contribution
- ✅ Single-transcript EC detection and logging
- ✅ Multi-transcript support

### Helper Tools Created
- ✅ `compare_traces.py` - Compare em_quant and Salmon traces
- ✅ `extract_trace.py` - Extract specific iterations/transcripts/ECs
- ✅ `validate_trace.py` - Validate trace file format
- ✅ `README_TRACING.md` - Comprehensive documentation

### Testing
- ✅ Single transcript tracing works
- ✅ Multi-transcript tracing works
- ✅ Trace format validated
- ✅ Helper scripts tested

## ⏳ Pending: Salmon Instrumentation

### Required Changes
1. **Add debug parameters to SalmonOpts**
   - `bool debugTrace`
   - `std::string debugTraceFile`
   - `std::vector<std::string> debugTranscripts`

2. **Add command-line flags**
   - `--debugTrace <file>`
   - `--debugTranscripts <comma-separated IDs>`

3. **Modify VBEMUpdate_ function**
   - Add debug parameters
   - Log per-iteration values (before and after E-step)
   - Log per-EC details

4. **Modify optimize function**
   - Setup debug stream
   - Parse debug transcript IDs
   - Pass to VBEMUpdate_
   - Log initial state (iter -1)

### Location
- **File**: `/mnt/pikachu/salmon/src/CollapsedEMOptimizer.cpp`
- **Function**: `VBEMUpdate_` (line 104 for single-threaded, line 241 for multi-threaded)
- **Main Loop**: `optimize` function (line 738)

### Guide Created
- ✅ `SALMON_INSTRUMENTATION_PATCH.md` - Detailed patch guide with code examples

## Next Steps

1. **Apply Salmon instrumentation** using the guide in `SALMON_INSTRUMENTATION_PATCH.md`
2. **Generate Salmon trace** with same fixture and transcript IDs
3. **Compare traces** using `compare_traces.py` to find first divergence
4. **Analyze divergence** to identify root cause (weights parsing, digamma, truncation, or convergence)

## Usage Example

Once Salmon instrumentation is added:

```bash
# Generate em_quant trace
export EM_QUANT_DEBUG_TXPS=ENST00000465752
em_quant --vb --threads 1 \
    -e eq_classes.txt -l quant.sf \
    -o em_output.tsv \
    --debug-trace /tmp/em_trace.txt

# Generate Salmon trace (after adding instrumentation)
salmon quant -i index -l A -1 r1.fq -2 r2.fq \
    --dumpEqWeights --threads 1 \
    --debugTrace /tmp/salmon_trace.txt \
    --debugTranscripts ENST00000465752 \
    -o salmon_out

# Compare traces
python3 compare_traces.py /tmp/em_trace.txt /tmp/salmon_trace.txt
```

# Salmon Instrumentation Debug - Quick Summary

## Problem

Salmon's error model instrumentation was implemented (per `INSTRUMENTATION_IMPLEMENTATION_SUMMARY.md`), but **trace files and matrix dumps are not being generated** when running with environment variables.

## What Should Happen

```bash
SALMON_ERROR_MODEL_TRACE=/tmp/trace.txt \
SALMON_TRACE_LEVEL=2 \
SALMON_DUMP_MATRICES=/tmp/matrices \
salmon quant -t transcriptome.fa -l A -a alignments.bam \
  --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
  --noFragLengthDist -p 1 -o output
```

**Expected**: 
- ✅ `/tmp/trace.txt` created with trace entries
- ✅ `/tmp/matrices_*.tsv` files created

**Actual**: 
- ❌ No trace file created
- ❌ No matrix files created

## Quick Debug Steps

1. **Check if binary includes instrumentation**:
   ```bash
   nm /mnt/pikachu/salmon/build/src/salmon | grep -i "trace\|dumpMatrices"
   ```

2. **Verify environment variables are read** in `AlignmentModel` constructor:
   - Check `/mnt/pikachu/salmon/src/AlignmentModel.cpp` constructor
   - Should read `SALMON_ERROR_MODEL_TRACE`, `SALMON_TRACE_LEVEL`, `SALMON_DUMP_MATRICES`

3. **Verify error model is enabled** in `SalmonQuantifyAlignments.cpp`:
   - Check `emTraceEnabled` flag is true
   - Check `setTraceContext()` is called
   - Check `traceRead()` is called

4. **Test with minimal case**:
   ```bash
   SALMON_ERROR_MODEL_TRACE=/tmp/test.txt SALMON_TRACE_LEVEL=1 \
   salmon quant -t transcriptome.fa -l A -a minimal.bam \
     --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
     --noFragLengthDist -p 1 -o out
   ls -lh /tmp/test.txt  # Should exist
   ```

## Expected Trace Format

If working, `/tmp/trace.txt` should contain:
```
READ SRR6357070.6710070 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
READ SRR6357070.5843839 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
```

For level 2, also:
```
ALN SRR6357070.6710070 tid=0 fg=-123.456 bg=-125.789 errLike=2.333 bin=2
```

## Reference

- **Full Debug Guide**: `SALMON_INSTRUMENTATION_DEBUG_HANDOFF.md`
- **Implementation Summary**: `/mnt/pikachu/salmon/INSTRUMENTATION_IMPLEMENTATION_SUMMARY.md`
- **Original Handoff**: `SALMON_HANDOFF.md`
- **Working Example**: `/tmp/error_model_test/our_trace.txt` (29,726 lines, working correctly)

## Files to Check

1. `/mnt/pikachu/salmon/src/AlignmentModel.cpp` - Constructor, trace methods
2. `/mnt/pikachu/salmon/include/AlignmentModel.hpp` - Public methods
3. `/mnt/pikachu/salmon/src/SalmonQuantifyAlignments.cpp` - Trace integration

## Success Criteria

✅ Trace file created at `SALMON_ERROR_MODEL_TRACE` path  
✅ Trace file contains READ entries (level 1+)  
✅ Matrix files created (if `SALMON_DUMP_MATRICES` set)  
✅ Format matches specification in `SALMON_HANDOFF.md`

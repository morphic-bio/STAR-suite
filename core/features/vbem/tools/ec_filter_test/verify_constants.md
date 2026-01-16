# Error Model Constants Verification

## State Constants

**Our code** (`source/libem/alignment_model.h:13-22`):
```cpp
enum AlignmentState : uint8_t {
    ALN_A = 0,
    ALN_C = 1,
    ALN_G = 2,
    ALN_T = 3,
    ALN_DASH = 4,
    ALN_SOFT_CLIP = 5,
    ALN_HARD_CLIP = 6,
    ALN_PAD = 7,
    ALN_REF_SKIP = 8
};
```

**Salmon** (`salmon/include/AlignmentCommon.hpp:42-52`):
```cpp
enum AlignmentModelChar {
    ALN_A = 0,
    ALN_C = 1,
    ALN_G = 2,
    ALN_T = 3,
    ALN_DASH = 4,
    ALN_SOFT_CLIP = 5,
    ALN_HARD_CLIP = 6,
    ALN_PAD = 7,
    ALN_REF_SKIP = 8
};
```

**Status**: ✅ **MATCH** - All values identical

## Start State Index

**Our code** (`source/libem/alignment_model.h:28`):
```cpp
constexpr uint32_t START_STATE_IDX = 81;
```

**Salmon** (`salmon/include/AlignmentModel.hpp:57`):
```cpp
constexpr static uint32_t startStateIdx = 81;
```

**Status**: ✅ **MATCH** - Both use 81

## Two-Bit Encoding

**Our code** (`source/libem/alignment_model.h:32-34`):
```cpp
constexpr uint8_t samToTwoBit[] = {
    0, /*A*/ 0, /*C*/ 1, 0, /*G*/ 2, 0, 0, 0, /*T*/ 3, 0, 0, 0, 0, 0, 0, 0
};
```

**Salmon** (`salmon/include/SalmonStringUtils.hpp:106-107`):
```cpp
constexpr uint8_t samToTwoBit[] = {
    0, /*A*/ 0, /*C*/ 1, 0, /*G*/ 2, 0, 0, 0, /*T*/ 3, 0, 0, 0, 0, 0, 0, 0};
```

**Status**: ✅ **MATCH** - Arrays identical

## Left/Right Read Determination

**Our code** (`source/libem/alignment_model.cpp:496`):
```cpp
bool read1IsLeft = (refPos1 <= refPos2);
```

**Salmon** (`salmon/src/AlignmentModel.cpp:247-250`):
```cpp
bam_seq_t* leftRead =
    (bam_pos(hit.read1) < bam_pos(hit.read2)) ? hit.read1 : hit.read2;
```

**Status**: ⚠️ **POTENTIAL DIFFERENCE** - We use `<=`, Salmon uses `<`

**Impact**: When `refPos1 == refPos2`, our code treats read1 as left, Salmon treats read2 as left. This could cause different bin assignments for paired reads starting at the same position.

**Recommendation**: Test with reads where both mates start at the same position to verify if this causes divergence.

## Number of States

**Our code** (`source/libem/alignment_model.h:26-27`):
```cpp
constexpr uint32_t NUM_ALIGNMENT_STATES = 82;
constexpr uint32_t NUM_BASE_STATES = 9;
```

**Salmon** (`salmon/include/AlignmentModel.hpp:55-56`):
```cpp
constexpr static uint32_t numAlignmentStates() { return 82; }
constexpr static uint32_t numStates = 9;
```

**Status**: ✅ **MATCH** - Both use 82 total states (81 + 1 start state), 9 base states

## Alpha Initialization

**Our code**: Default `alpha = 0.001` in default constructor, `alpha = 1.0` when explicitly passed (as in CLI)

**Salmon**: Need to verify default alpha value

**Status**: ⚠️ **NEEDS VERIFICATION** - Check Salmon's default alpha

## Read Bins

**Our code**: Default `readBins = 4`, CLI uses `readBins = 6` (Salmon default)

**Salmon**: Default `readBins = 4` in constructor signature

**Status**: ⚠️ **NEEDS VERIFICATION** - Check if Salmon uses 6 bins by default in practice

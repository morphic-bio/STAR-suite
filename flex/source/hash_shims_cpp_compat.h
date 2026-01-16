#ifndef HASH_SHIMS_CPP_COMPAT_H
#define HASH_SHIMS_CPP_COMPAT_H

// C++ compatibility wrapper for hash_shims.h
// We only need cr_hash_t and cg_agg types, not bitset functions
// This avoids C++ compatibility issues with bitset_new (which uses malloc without casts)
// Note: Include paths are set via Makefile (-I flags)

#include <cstdlib>
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

extern "C" {
// Include klib first (needed for hash types)
// Path is set via Makefile: -I../tools/flex_debug/flex/third_party
#include "klib/khash.h"

// Forward declare bitset_t structure (matches hash_shims.h definition)
// We don't use bitset functions, but some headers might reference the type
typedef struct {
    uint64_t *words;
    size_t nbits;
    size_t nwords;
} bitset_t;

// Only define the hash types we actually use
// 1. cr_hash: Cell Ranger compatible molecule hash (compact 64-bit key -> uint16_t count)
KHASH_MAP_INIT_INT64(cr, uint16_t)

// 2. cg_agg: Map uint64_t -> uint32_t (cell|gene key -> count)
KHASH_MAP_INIT_INT64(cg_agg, uint32_t)

// Type aliases
typedef khash_t(cr) cr_hash_t;

// CR hash functions (C++ compatible with explicit casts)
static inline cr_hash_t* cr_hash_new(void) {
    return kh_init(cr);
}

static inline void cr_hash_free(cr_hash_t* map) {
    if (map) kh_destroy(cr, map);
}

static inline void cr_hash_reserve(cr_hash_t* map, size_t n) {
    if (map) kh_resize(cr, map, n);
}

static inline int cr_hash_key_is_unmatched(uint64_t key) {
    return (key & (1ULL << 63)) != 0;
}

static inline void cr_hash_unpack(uint64_t key,
                                  uint32_t *cell_idx,
                                  uint32_t *umi24,
                                  uint16_t *gene_idx) {
    // Mask out bit 63 (unmatched flag) before unpacking
    key &= ~(1ULL << 63);
    if (cell_idx) *cell_idx = (uint32_t)((key >> 40) & 0xFFFFFF);
    if (umi24) *umi24 = (uint32_t)((key >> 16) & 0xFFFFFF);
    if (gene_idx) *gene_idx = (uint16_t)(key & 0x7FFF);
}

static inline void cr_hash_increment(cr_hash_t* map, uint64_t key) {
    if (!map) return;
    int absent;
    khiter_t iter = kh_put(cr, map, key, &absent);
    if (absent) {
        kh_val(map, iter) = 1;
    } else if (kh_val(map, iter) < UINT16_MAX) {
        kh_val(map, iter)++;
    }
    // Silently cap at UINT16_MAX
}

static inline uint16_t cr_hash_get(cr_hash_t* map, uint64_t key) {
    if (!map) return 0;
    khiter_t iter = kh_get(cr, map, key);
    return (iter != kh_end(map)) ? kh_val(map, iter) : 0;
}

static inline size_t cr_hash_size(const cr_hash_t* map) {
    return map ? kh_size(map) : 0;
}

// Iterator macro for cr_hash
#define cr_hash_foreach(map, key, val) \
    for (khiter_t __i = kh_begin(map); __i != kh_end(map); ++__i) \
        if (kh_exist(map, __i) && ((key) = kh_key(map, __i), (val) = kh_val(map, __i), 1))

} // extern "C"

// Inline hash key pack/unpack functions for cg_agg hash
// Key format: [CB20][UMI24][GENE15][TAG5] MSB→LSB
// CB: 20 bits (0-based whitelist index, covers 737k WL)
// UMI: 24 bits (packed UMI12)
// Gene: 15 bits (gene index, 0 = no gene)
// Tag: 5 bits (tag whitelist index, 0 = no tag, 1-31 = whitelist index)

static inline uint64_t packCgAggKey(uint32_t cbIdx, uint32_t umi24, uint16_t geneIdx, uint8_t tagIdx) {
    return ((uint64_t)(cbIdx & 0xFFFFF) << 44) |
           ((uint64_t)(umi24 & 0xFFFFFF) << 20) |
           ((uint64_t)(geneIdx & 0x7FFF) << 5) |
           ((uint64_t)(tagIdx & 0x1F));
}

static inline void unpackCgAggKey(uint64_t key, uint32_t *cbIdx, uint32_t *umi24, uint16_t *geneIdx, uint8_t *tagIdx) {
    if (cbIdx) *cbIdx = (uint32_t)((key >> 44) & 0xFFFFF);
    if (umi24) *umi24 = (uint32_t)((key >> 20) & 0xFFFFFF);
    if (geneIdx) *geneIdx = (uint16_t)((key >> 5) & 0x7FFF);
    if (tagIdx) *tagIdx = (uint8_t)(key & 0x1F);
}

#endif // HASH_SHIMS_CPP_COMPAT_H

#ifndef KHASH_WRAPPER_H
#define KHASH_WRAPPER_H

#include <stdint.h>
#include <string.h>
#include "../../../../flex/source/klib/khash.h"

/* ============================================================================
 * Common hash table types
 * ============================================================================ */

/* Hash table: uint32_t -> void* (for direct pointers) */
KHASH_INIT(u32ptr, uint32_t, void*, 1, kh_int_hash_func, kh_int_hash_equal)

/* Hash table: uint64_t -> void* */
KHASH_INIT(u64ptr, uint64_t, void*, 1, kh_int64_hash_func, kh_int64_hash_equal)

/* Hash table: uint32_t -> uint32_t */
KHASH_INIT(u32u32, uint32_t, uint32_t, 1, kh_int_hash_func, kh_int_hash_equal)

/* Hash table: uint64_t -> uint64_t */
KHASH_INIT(u64u64, uint64_t, uint64_t, 1, kh_int64_hash_func, kh_int64_hash_equal)

/* Hash table: char* -> void* (string keys) */
KHASH_INIT(strptr, const char*, void*, 1, kh_str_hash_func, kh_str_hash_equal)

/* Hash table: char* -> char* (string -> string) */
KHASH_INIT(strstr, const char*, char*, 1, kh_str_hash_func, kh_str_hash_equal)

/* Hash table: char* -> uint32_t (string -> uint32_t, used by prehash) */
KHASH_INIT(stru32, const char*, uint32_t, 1, kh_str_hash_func, kh_str_hash_equal)

/* ============================================================================
 * Integer-key hash tables for feature sequence encoding
 *
 * 2-bit encoding: A=00, C=01, G=10, T=11
 * Fixed-length mode (prehash, exact hash when all features same length):
 *   pure 2-bit encoding, no length bits.
 * Variable-length mode (exact hash when features have mixed lengths):
 *   top 6 bits encode (len - 1), remaining bits encode sequence.
 * ============================================================================ */

/* Hash table: uint64_t -> uint32_t (for <=32bp fixed or <=29bp variable) */
KHASH_INIT(u64u32, uint64_t, uint32_t, 1, kh_int64_hash_func, kh_int64_hash_equal)

/* 128-bit sequence key (for 33-64bp fixed or 30-61bp variable) */
typedef struct { uint64_t w[2]; } seq128_t;

static inline khint_t seq128_hash_func(seq128_t k) {
    khint_t h = (khint_t)(k.w[0] ^ (k.w[0] >> 33));
    h *= 0xff51afd7ed558ccdULL;
    h ^= (khint_t)(k.w[1] ^ (k.w[1] >> 33));
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 16;
    return h;
}

static inline int seq128_hash_equal(seq128_t a, seq128_t b) {
    return a.w[0] == b.w[0] && a.w[1] == b.w[1];
}

#define seq128_hash(k) seq128_hash_func(k)
#define seq128_eq(a, b) seq128_hash_equal(a, b)
KHASH_INIT(seq128u32, seq128_t, uint32_t, 1, seq128_hash, seq128_eq)

/* Runtime mode selection and dispatch union */
typedef enum { SEQ_KEY_64, SEQ_KEY_128 } seq_key_mode_t;

typedef union {
    khash_t(u64u32)    *h64;
    khash_t(seq128u32) *h128;
} seq_hash_t;

#define MAX_VARIABLE_KEY_SEQ_LENGTH_64  29
#define MAX_VARIABLE_KEY_SEQ_LENGTH_128 61
#define MAX_FIXED_KEY_SEQ_LENGTH_64     32
#define MAX_FIXED_KEY_SEQ_LENGTH_128    64

/* ============================================================================
 * Sequence encoding helpers
 * ============================================================================ */

/*
 * These require seq2code[] to be initialized (via initseq2Code()).
 * The seq2code table is declared extern in globals.h.
 */

static inline uint64_t seq_encode_64_fixed(const char *seq, int len) {
    uint64_t key = 0;
    for (int i = 0; i < len; i++) {
        key = (key << 2) | (uint64_t)(((const unsigned char *)seq)[i] == 'C' ? 1 :
                                       ((const unsigned char *)seq)[i] == 'G' ? 2 :
                                       ((const unsigned char *)seq)[i] == 'T' ? 3 : 0);
    }
    return key;
}

static inline uint64_t seq_encode_64_var(const char *seq, int len) {
    uint64_t key = (uint64_t)(len - 1) << 58;
    for (int i = 0; i < len; i++) {
        key |= ((uint64_t)(((const unsigned char *)seq)[i] == 'C' ? 1 :
                            ((const unsigned char *)seq)[i] == 'G' ? 2 :
                            ((const unsigned char *)seq)[i] == 'T' ? 3 : 0)) << (2 * (28 - i));
    }
    return key;
}

static inline seq128_t seq_encode_128_fixed(const char *seq, int len) {
    seq128_t key = {{0, 0}};
    for (int i = 0; i < len; i++) {
        uint64_t bits = ((const unsigned char *)seq)[i] == 'C' ? 1 :
                        ((const unsigned char *)seq)[i] == 'G' ? 2 :
                        ((const unsigned char *)seq)[i] == 'T' ? 3 : 0;
        int bit_pos = 2 * (63 - i);
        if (bit_pos >= 64) {
            key.w[0] |= bits << (bit_pos - 64);
        } else {
            key.w[1] |= bits << bit_pos;
        }
    }
    return key;
}

static inline seq128_t seq_encode_128_var(const char *seq, int len) {
    seq128_t key = {{0, 0}};
    key.w[0] = (uint64_t)(len - 1) << 58;
    for (int i = 0; i < len; i++) {
        uint64_t bits = ((const unsigned char *)seq)[i] == 'C' ? 1 :
                        ((const unsigned char *)seq)[i] == 'G' ? 2 :
                        ((const unsigned char *)seq)[i] == 'T' ? 3 : 0;
        /* Bits go into positions [57..0] of w[0] then [63..0] of w[1] */
        int bit_pos = 2 * (60 - i);  /* 60 because 6 bits reserved in w[0] top */
        if (bit_pos >= 64) {
            key.w[0] |= bits << (bit_pos - 64);
        } else if (bit_pos >= 0) {
            key.w[1] |= bits << bit_pos;
        }
    }
    return key;
}

/* ============================================================================
 * Dispatch wrappers for seq_hash_t
 * ============================================================================ */

static inline void seq_hash_init(seq_hash_t *sh, seq_key_mode_t mode) {
    if (mode == SEQ_KEY_64) {
        sh->h64 = kh_init(u64u32);
    } else {
        sh->h128 = kh_init(seq128u32);
    }
}

static inline void seq_hash_destroy(seq_hash_t *sh, seq_key_mode_t mode) {
    if (mode == SEQ_KEY_64) {
        if (sh->h64) { kh_destroy(u64u32, sh->h64); sh->h64 = NULL; }
    } else {
        if (sh->h128) { kh_destroy(seq128u32, sh->h128); sh->h128 = NULL; }
    }
}

static inline int seq_hash_put_64(seq_hash_t *sh, uint64_t key, uint32_t val) {
    int ret;
    khint_t k = kh_put(u64u32, sh->h64, key, &ret);
    if (ret < 0) return -1;
    kh_val(sh->h64, k) = val;
    return ret;
}

static inline int seq_hash_put_128(seq_hash_t *sh, seq128_t key, uint32_t val) {
    int ret;
    khint_t k = kh_put(seq128u32, sh->h128, key, &ret);
    if (ret < 0) return -1;
    kh_val(sh->h128, k) = val;
    return ret;
}

static inline uint32_t seq_hash_get_64(const seq_hash_t *sh, uint64_t key) {
    khint_t k = kh_get(u64u32, sh->h64, key);
    return (k != kh_end(sh->h64)) ? kh_val(sh->h64, k) : 0;
}

static inline uint32_t seq_hash_get_128(const seq_hash_t *sh, seq128_t key) {
    khint_t k = kh_get(seq128u32, sh->h128, key);
    return (k != kh_end(sh->h128)) ? kh_val(sh->h128, k) : 0;
}

static inline khint_t seq_hash_size(const seq_hash_t *sh, seq_key_mode_t mode) {
    if (mode == SEQ_KEY_64) {
        return sh->h64 ? kh_size(sh->h64) : 0;
    } else {
        return sh->h128 ? kh_size(sh->h128) : 0;
    }
}

/* ============================================================================
 * Variable-length key structure (replaces GBytes)
 * ============================================================================ */
typedef struct {
    uint8_t *ptr;
    uint16_t len;
} var_key_t;

/* Hash function for variable-length keys */
static inline khint_t var_key_hash_func(const var_key_t *k) {
    khint_t h = k->len;
    const uint8_t *p = k->ptr;
    for (uint16_t i = 0; i < k->len; ++i) {
        h = (h << 5) - h + (khint_t)p[i];
    }
    return h;
}

/* Equality function for variable-length keys */
static inline int var_key_equal(const var_key_t *a, const var_key_t *b) {
    if (a->len != b->len) return 0;
    return memcmp(a->ptr, b->ptr, a->len) == 0;
}

/* Hash table: var_key_t -> uint32_t (for feature_code_hash) */
#define var_key_hash(k) var_key_hash_func(&(k))
#define var_key_eq(a, b) var_key_equal(&(a), &(b))
KHASH_INIT(codeu32, var_key_t, uint32_t, 1, var_key_hash, var_key_eq)

/* Hash table: var_key_t -> void* */
KHASH_INIT(codeptr, var_key_t, void*, 1, var_key_hash, var_key_eq)

/* ============================================================================
 * Helper macros for common operations
 * ============================================================================ */

/* Lookup and return value, or NULL if not found */
#define kh_get_val(name, h, k, def) \
    ({ khint_t __k = kh_get(name, h, k); \
       (__k != kh_end(h)) ? kh_val(h, __k) : (def); })

/* Insert or replace */
#define kh_put_replace(name, h, k, v) \
    ({ int __ret; \
       khint_t __k = kh_put(name, h, k, &__ret); \
       kh_val(h, __k) = (v); \
       __ret; })

/* ============================================================================
 * Simple dynamic array (kvec-style, replaces GArray)
 * ============================================================================ */
typedef struct {
    uint32_t *a;
    size_t n;
    size_t m;
} vec_u32_t;

static inline vec_u32_t* vec_u32_init(void) {
    vec_u32_t *v = (vec_u32_t*)calloc(1, sizeof(vec_u32_t));
    return v;
}

static inline void vec_u32_destroy(vec_u32_t *v) {
    if (v) {
        free(v->a);
        free(v);
    }
}

static inline void vec_u32_clear(vec_u32_t *v) {
    if (v) v->n = 0;
}

static inline void vec_u32_resize(vec_u32_t *v, size_t new_size) {
    if (new_size > v->m) {
        size_t new_m = v->m ? v->m * 2 : 8;
        while (new_m < new_size) new_m *= 2;
        v->a = (uint32_t*)realloc(v->a, new_m * sizeof(uint32_t));
        memset(v->a + v->m, 0, (new_m - v->m) * sizeof(uint32_t));
        v->m = new_m;
    }
    if (new_size > v->n) {
        memset(v->a + v->n, 0, (new_size - v->n) * sizeof(uint32_t));
    }
    v->n = new_size;
}

static inline void vec_u32_set_size(vec_u32_t *v, size_t size) {
    vec_u32_resize(v, size);
}

static inline uint32_t vec_u32_get(vec_u32_t *v, size_t i) {
    return v->a[i];
}

static inline void vec_u32_set(vec_u32_t *v, size_t i, uint32_t val) {
    if (i >= v->n) vec_u32_resize(v, i + 1);
    v->a[i] = val;
}

static inline void vec_u32_inc(vec_u32_t *v, size_t i) {
    if (i >= v->n) vec_u32_resize(v, i + 1);
    v->a[i]++;
}

static inline size_t vec_u32_size(vec_u32_t *v) {
    return v ? v->n : 0;
}

/* Simple dynamic pointer array (replaces GPtrArray) */
typedef struct {
    void **a;
    size_t n;
    size_t m;
    void (*free_func)(void*);
} vec_ptr_t;

static inline vec_ptr_t* vec_ptr_init(void (*free_func)(void*)) {
    vec_ptr_t *v = (vec_ptr_t*)calloc(1, sizeof(vec_ptr_t));
    v->free_func = free_func;
    return v;
}

static inline void vec_ptr_destroy(vec_ptr_t *v) {
    if (v) {
        if (v->free_func) {
            for (size_t i = 0; i < v->n; i++) {
                if (v->a[i]) v->free_func(v->a[i]);
            }
        }
        free(v->a);
        free(v);
    }
}

static inline void vec_ptr_add(vec_ptr_t *v, void *item) {
    if (v->n >= v->m) {
        size_t new_m = v->m ? v->m * 2 : 8;
        v->a = (void**)realloc(v->a, new_m * sizeof(void*));
        v->m = new_m;
    }
    v->a[v->n++] = item;
}

static inline void* vec_ptr_get(vec_ptr_t *v, size_t i) {
    return (i < v->n) ? v->a[i] : NULL;
}

static inline size_t vec_ptr_size(vec_ptr_t *v) {
    return v ? v->n : 0;
}

/* Simple barcode list for feature->barcodes mapping (replaces GList) */
typedef struct {
    uint32_t *keys;
    size_t n;
    size_t m;
} barcode_list_t;

static inline barcode_list_t* barcode_list_init(void) {
    barcode_list_t *bl = (barcode_list_t*)calloc(1, sizeof(barcode_list_t));
    return bl;
}

static inline void barcode_list_destroy(barcode_list_t *bl) {
    if (bl) {
        free(bl->keys);
        free(bl);
    }
}

static inline void barcode_list_add(barcode_list_t *bl, uint32_t key) {
    if (bl->n >= bl->m) {
        size_t new_m = bl->m ? bl->m * 2 : 8;
        bl->keys = (uint32_t*)realloc(bl->keys, new_m * sizeof(uint32_t));
        bl->m = new_m;
    }
    bl->keys[bl->n++] = key;
}

#endif /* KHASH_WRAPPER_H */

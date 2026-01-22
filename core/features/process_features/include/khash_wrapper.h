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

/* Hash table: char* -> uint32_t (string -> uint32_t) */
KHASH_INIT(stru32, const char*, uint32_t, 1, kh_str_hash_func, kh_str_hash_equal)

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
    vec_u32_t *v = calloc(1, sizeof(vec_u32_t));
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
        v->a = realloc(v->a, new_m * sizeof(uint32_t));
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
    vec_ptr_t *v = calloc(1, sizeof(vec_ptr_t));
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
        v->a = realloc(v->a, new_m * sizeof(void*));
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
    barcode_list_t *bl = calloc(1, sizeof(barcode_list_t));
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
        bl->keys = realloc(bl->keys, new_m * sizeof(uint32_t));
        bl->m = new_m;
    }
    bl->keys[bl->n++] = key;
}

#endif /* KHASH_WRAPPER_H */

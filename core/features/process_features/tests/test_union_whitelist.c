/**
 * Test expand_hash_union_namespace() and filtered_barcode_hash_contains()
 * with mixed NXT+TRU barcode sets.
 *
 * NXT/TRU translation: complement positions 7 and 8 (0-based).
 *
 * Example barcode pair (16 bases):
 *   TRU: AAACCCAAGTAACCCC   positions 7,8 = A,G
 *   NXT: AAACCCATCTAACCCC   positions 7,8 = T,C (complement of A,G)
 *
 * Tests:
 *   1. expand_hash_union_namespace adds missing NXT/TRU translations.
 *   2. After expansion, filtered_barcode_hash_contains finds both forms.
 *   3. Both forms already present → no expansion needed.
 *   4. Barcodes too short for translation (<9 bases) are left alone.
 *   5. Multiple mixed barcodes expand correctly.
 *   6. Runtime fallback in filtered_barcode_hash_contains works without expansion.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/common.h"
#include "../include/prototypes.h"
#include "../include/io.h"
#include "../include/pf_api.h"

static int failures = 0;

#define ASSERT_EQ_INT(msg, got, expected) do { \
    int _g = (got), _e = (expected); \
    if (_g != _e) { \
        fprintf(stderr, "  FAIL: %s: got %d, expected %d\n", (msg), _g, _e); \
        failures++; \
    } else { \
        fprintf(stderr, "  PASS: %s\n", (msg)); \
    } \
} while(0)

static khash_t(strptr)* make_hash_from_list(const char **barcodes, int n) {
    khash_t(strptr) *h = kh_init(strptr);
    for (int i = 0; i < n; i++) {
        char *dup = strdup(barcodes[i]);
        int ret;
        kh_put(strptr, h, dup, &ret);
    }
    return h;
}

static int hash_size(khash_t(strptr) *h) {
    return (int)kh_size(h);
}

/*
 * Barcode reference for all tests:
 *   TRU1: AAACCCAAGTAACCCC  (pos7=A, pos8=G)
 *   NXT1: AAACCCATCTAACCCC  (pos7=T, pos8=C — complement)
 *
 *   TRU2: TTTTCCAAGTCCCCCC  (pos7=A, pos8=G)
 *   NXT2: TTTTCCATCTCCCCCC  (pos7=T, pos8=C — complement)
 *
 *   BC3:  GGGGCCAACACCCCCC  (pos7=A, pos8=C)
 *   BC3T: GGGGCCATGACCCCCC  (pos7=T, pos8=G — complement)
 */

static void test_expansion_adds_translations(void) {
    fprintf(stderr, "\n--- test_expansion_adds_translations ---\n");

    const char *barcodes[] = {
        "AAACCCAAGTAACCCC",  /* TRU1 only */
    };
    khash_t(strptr) *h = make_hash_from_list(barcodes, 1);
    ASSERT_EQ_INT("initial size 1", hash_size(h), 1);

    int added = expand_hash_union_namespace(h);
    ASSERT_EQ_INT("added 1 translation", added, 1);
    ASSERT_EQ_INT("size after expansion", hash_size(h), 2);

    ASSERT_EQ_INT("contains original TRU1", filtered_barcode_hash_contains(h, "AAACCCAAGTAACCCC"), 1);
    ASSERT_EQ_INT("contains NXT1 translation", filtered_barcode_hash_contains(h, "AAACCCATCTAACCCC"), 1);

    free_strptr_hash(h);
}

static void test_expansion_both_forms_present(void) {
    fprintf(stderr, "\n--- test_expansion_both_forms_present ---\n");

    const char *barcodes[] = {
        "AAACCCAAGTAACCCC",  /* TRU1 */
        "AAACCCATCTAACCCC",  /* NXT1 — correct translation pair */
    };
    khash_t(strptr) *h = make_hash_from_list(barcodes, 2);
    ASSERT_EQ_INT("initial size 2", hash_size(h), 2);

    int added = expand_hash_union_namespace(h);
    ASSERT_EQ_INT("added 0 (both already present)", added, 0);
    ASSERT_EQ_INT("size unchanged", hash_size(h), 2);

    free_strptr_hash(h);
}

static void test_expansion_short_barcode_skipped(void) {
    fprintf(stderr, "\n--- test_expansion_short_barcode_skipped ---\n");

    const char *barcodes[] = {
        "AAACCCAA",  /* 8 chars — too short for NXT translation (needs >= 9) */
    };
    khash_t(strptr) *h = make_hash_from_list(barcodes, 1);

    int added = expand_hash_union_namespace(h);
    ASSERT_EQ_INT("added 0 for short barcode", added, 0);
    ASSERT_EQ_INT("size unchanged", hash_size(h), 1);

    free_strptr_hash(h);
}

static void test_expansion_multiple_barcodes(void) {
    fprintf(stderr, "\n--- test_expansion_multiple_barcodes ---\n");

    const char *barcodes[] = {
        "AAACCCAAGTAACCCC",  /* TRU1 only */
        "TTTTCCAAGTCCCCCC",  /* TRU2 only */
        "GGGGCCAACACCCCCC",  /* BC3 only */
    };
    khash_t(strptr) *h = make_hash_from_list(barcodes, 3);
    ASSERT_EQ_INT("initial size 3", hash_size(h), 3);

    int added = expand_hash_union_namespace(h);
    ASSERT_EQ_INT("added 3 translations", added, 3);
    ASSERT_EQ_INT("size after expansion", hash_size(h), 6);

    ASSERT_EQ_INT("TRU1 present", filtered_barcode_hash_contains(h, "AAACCCAAGTAACCCC"), 1);
    ASSERT_EQ_INT("NXT1 present", filtered_barcode_hash_contains(h, "AAACCCATCTAACCCC"), 1);
    ASSERT_EQ_INT("TRU2 present", filtered_barcode_hash_contains(h, "TTTTCCAAGTCCCCCC"), 1);
    ASSERT_EQ_INT("NXT2 present", filtered_barcode_hash_contains(h, "TTTTCCATCTCCCCCC"), 1);
    ASSERT_EQ_INT("BC3 present",  filtered_barcode_hash_contains(h, "GGGGCCAACACCCCCC"), 1);
    ASSERT_EQ_INT("BC3T present", filtered_barcode_hash_contains(h, "GGGGCCATGACCCCCC"), 1);

    free_strptr_hash(h);
}

static void test_exact_only_without_expansion(void) {
    fprintf(stderr, "\n--- test_exact_only_without_expansion ---\n");

    const char *barcodes[] = {
        "AAACCCAAGTAACCCC",  /* TRU1 only — NXT1 not in hash */
    };
    khash_t(strptr) *h = make_hash_from_list(barcodes, 1);

    /* Direct hit */
    ASSERT_EQ_INT("exact match TRU1", filtered_barcode_hash_contains(h, "AAACCCAAGTAACCCC"), 1);
    /* Strict exact-only: NXT1 is NOT in hash, must NOT be found */
    ASSERT_EQ_INT("NXT1 not found (exact-only)", filtered_barcode_hash_contains(h, "AAACCCATCTAACCCC"), 0);
    /* Unknown barcode */
    ASSERT_EQ_INT("unrelated not found", filtered_barcode_hash_contains(h, "TTTTTTTTTTTTTTTT"), 0);

    free_strptr_hash(h);
}

static void test_idempotent_double_expansion(void) {
    fprintf(stderr, "\n--- test_idempotent_double_expansion ---\n");

    const char *barcodes[] = {
        "AAACCCAAGTAACCCC",
    };
    khash_t(strptr) *h = make_hash_from_list(barcodes, 1);

    int added1 = expand_hash_union_namespace(h);
    ASSERT_EQ_INT("first expansion adds 1", added1, 1);
    ASSERT_EQ_INT("size after first expansion", hash_size(h), 2);

    int added2 = expand_hash_union_namespace(h);
    ASSERT_EQ_INT("second expansion adds 0 (idempotent)", added2, 0);
    ASSERT_EQ_INT("size unchanged after second expansion", hash_size(h), 2);

    free_strptr_hash(h);
}

static void test_translate_involution(void) {
    fprintf(stderr, "\n--- test_translate_involution ---\n");

    char bc[] = "AAACCCAAGTAACCCC";  /* 16 bases */
    char original[17];
    memcpy(original, bc, 17);

    pf_translate_barcode_inplace(bc, 16);
    ASSERT_EQ_INT("pos7 complemented", bc[7] == 'T', 1);
    ASSERT_EQ_INT("pos8 complemented", bc[8] == 'C', 1);

    pf_translate_barcode_inplace(bc, 16);
    ASSERT_EQ_INT("involution restores original", strcmp(bc, original) == 0, 1);

    /* Short barcode: no-op */
    char short_bc[] = "ACGTACGT";  /* 8 bases */
    char short_orig[9];
    memcpy(short_orig, short_bc, 9);
    pf_translate_barcode_inplace(short_bc, 8);
    ASSERT_EQ_INT("short barcode unchanged", strcmp(short_bc, short_orig) == 0, 1);

    /* NULL: no crash */
    pf_translate_barcode_inplace(NULL, 16);
    fprintf(stderr, "  PASS: NULL no crash\n");
}

static void test_namespace_enum_helpers(void) {
    fprintf(stderr, "\n--- test_namespace_enum_helpers ---\n");

    ASSERT_EQ_INT("NXT round-trip", pf_namespace_from_string("NXT"), PF_NS_NXT);
    ASSERT_EQ_INT("TRU round-trip", pf_namespace_from_string("TRU"), PF_NS_TRU);
    ASSERT_EQ_INT("UNION round-trip", pf_namespace_from_string("UNION"), PF_NS_UNION);
    ASSERT_EQ_INT("UNKNOWN round-trip", pf_namespace_from_string("UNKNOWN"), PF_NS_UNKNOWN);
    ASSERT_EQ_INT("junk is UNKNOWN", pf_namespace_from_string("junk"), PF_NS_UNKNOWN);
    ASSERT_EQ_INT("NULL is UNKNOWN", pf_namespace_from_string(NULL), PF_NS_UNKNOWN);
    ASSERT_EQ_INT("lowercase nxt", pf_namespace_from_string("nxt"), PF_NS_NXT);

    ASSERT_EQ_INT("to_string NXT", strcmp(pf_namespace_to_string(PF_NS_NXT), "NXT") == 0, 1);
    ASSERT_EQ_INT("to_string TRU", strcmp(pf_namespace_to_string(PF_NS_TRU), "TRU") == 0, 1);
    ASSERT_EQ_INT("to_string UNION", strcmp(pf_namespace_to_string(PF_NS_UNION), "UNION") == 0, 1);
    ASSERT_EQ_INT("to_string UNKNOWN", strcmp(pf_namespace_to_string(PF_NS_UNKNOWN), "UNKNOWN") == 0, 1);
}

int main(void) {
    test_expansion_adds_translations();
    test_expansion_both_forms_present();
    test_expansion_short_barcode_skipped();
    test_expansion_multiple_barcodes();
    test_exact_only_without_expansion();
    test_idempotent_double_expansion();
    test_translate_involution();
    test_namespace_enum_helpers();

    fprintf(stderr, "\n");
    if (failures > 0) {
        fprintf(stderr, "FAILED: %d assertion(s) failed\n", failures);
        return 1;
    }
    fprintf(stderr, "ALL PASSED\n");
    return 0;
}

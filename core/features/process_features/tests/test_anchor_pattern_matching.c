#include "../include/barcode_match.h"
#include "../include/globals.h"
#include "../include/io.h"
#include "../include/prototypes.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define TEST_ASSERT(cond, msg) do { \
    if (!(cond)) { \
        fprintf(stderr, "FAIL: %s (line %d)\n", msg, __LINE__); \
        return 1; \
    } \
} while (0)

static int write_single_feature_csv(const char *csv_path, const char *pattern, const char *sequence) {
    FILE *f = fopen(csv_path, "w");
    if (!f) {
        return -1;
    }
    fprintf(f, "name,sequence,pattern\n");
    fprintf(f, "Guide1,%s,%s\n", sequence, pattern);
    fclose(f);
    return 0;
}

static int run_anchor_case(const char *tmp_dir, const char *case_name, const char *pattern, const char *feature_seq, const char *read_seq, int expected_match_pos) {
    char csv_path[512];
    snprintf(csv_path, sizeof(csv_path), "%s/%s.csv", tmp_dir, case_name);
    TEST_ASSERT(write_single_feature_csv(csv_path, pattern, feature_seq) == 0, "write feature CSV");

    feature_arrays *features = read_features_file(csv_path);
    TEST_ASSERT(features != NULL, "read features");
    TEST_ASSERT(features->number_of_features == 1, "single feature loaded");

    use_feature_anchor_search = 1;
    require_feature_anchor_match = 1;
    use_feature_offset_array = 0;
    feature_mode_bootstrap_reads = 0;
    feature_mode_reads_seen = 0;
    feature_mode_bootstrap_done = 0;
    limit_search = -1;

    uint32_t feature_index = 0;
    int hamming_distance = -1;
    uint16_t match_position = 0;
    char matching_sequence[LINE_LENGTH];
    char read_buf[LINE_LENGTH];
    strncpy(read_buf, read_seq, LINE_LENGTH - 1);
    read_buf[LINE_LENGTH - 1] = '\0';

    process_feature_sequence(read_buf, features, 0, 1, 0, 1, &feature_index, &hamming_distance, matching_sequence, &match_position, NULL);

    TEST_ASSERT(feature_index == 1, "feature assigned");
    TEST_ASSERT(hamming_distance == 0, "zero hamming");
    TEST_ASSERT((int)match_position == expected_match_pos, "match position");
    TEST_ASSERT(strcmp(matching_sequence, feature_seq) == 0, "matched sequence");

    free_feature_arrays(features);
    feature_code_hash.h64 = NULL;
    unlink(csv_path);
    return 0;
}

int main(void) {
    barcode_match_init();
    initialize_complement();
    initialize_unit_sizes();

    char tmp_template[] = "/tmp/pf_anchor_test_XXXXXX";
    char *tmp_dir = mkdtemp(tmp_template);
    TEST_ASSERT(tmp_dir != NULL, "mkdtemp");

    printf("Test: prefix anchor pattern... ");
    fflush(stdout);
    TEST_ASSERT(run_anchor_case(
        tmp_dir,
        "prefix_only",
        "AAAA(BC)",
        "CCGTAACCGGTTAACCGGTT",
        "TTTTAAAACCGTAACCGGTTAACCGGTTGGGG",
        8
    ) == 0, "prefix anchor case");
    printf("PASS\n");

    printf("Test: suffix anchor at BC=0 with wildcard N... ");
    fflush(stdout);
    TEST_ASSERT(run_anchor_case(
        tmp_dir,
        "suffix_only_bc_zero",
        "(BC)GTTTNAGAGCTAAGC",
        "ACGTACGTACGTACGTACGT",
        "GGACGTACGTACGTACGTACGTGTTTAAGAGCTAAGCTT",
        2
    ) == 0, "suffix anchor case");
    printf("PASS\n");

    printf("Test: prefix+suffix anchor pattern... ");
    fflush(stdout);
    TEST_ASSERT(run_anchor_case(
        tmp_dir,
        "prefix_suffix",
        "AAA(BC)TTN",
        "GGCCTTAAGG",
        "CCAAAGGCCTTAAGGTTAGG",
        5
    ) == 0, "prefix+suffix anchor case");
    printf("PASS\n");

    rmdir(tmp_dir);
    printf("All anchor pattern tests passed.\n");
    return 0;
}

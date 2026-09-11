#include "../include/barcode_match.h"
#include "../include/globals.h"
#include "../include/io.h"
#include "../include/prototypes.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define CHECK(condition, message) do { \
    if (!(condition)) { fprintf(stderr, "FAIL: %s at line %d\n", message, __LINE__); return 1; } \
} while (0)

/* Actual A375 singleton reads: the old maxN=0 allocation could overwrite the
 * input with the match offset and return position zero. No reader or merge is
 * involved in this test. */
static int check_read(feature_arrays *features, const char *original, unsigned expected_pos) {
    char read[LINE_LENGTH], match[LINE_LENGTH];
    strcpy(read, original);
    int distance = -1;
    uint16_t position = 999;
    char ambiguous = 0;
    int feature = checkAndCorrectFeature(read, features, 1, 1, &distance, match,
                                        0, &ambiguous, &position, NULL);
    CHECK(strcmp(read, original) == 0, "fallback preserves the input read");
    CHECK(feature == 1 && distance == 1 && !ambiguous, "unique Hamming-1 match");
    CHECK(position == expected_pos, "fallback reports the actual match position");
    CHECK(strlen(match) == 7 && strncmp(match, original + expected_pos, 7) == 0,
          "matching sequence is copied from the reported position");
    return 0;
}

int main(void) {
    barcode_match_init();
    initialize_complement();
    initialize_unit_sizes();
    use_feature_offset_array = 0;
    char path[] = "/tmp/pf_zero_n_fallback_XXXXXX";
    int fd = mkstemp(path);
    CHECK(fd >= 0, "create feature fixture");
    FILE *file = fdopen(fd, "w");
    CHECK(file != NULL, "open feature fixture");
    fputs("name,sequence\nIL1B_sg2_HEK,TGAACCA\n", file);
    fclose(file);
    feature_arrays *features = read_features_file(path);
    unlink(path);
    CHECK(features != NULL, "load feature fixture");

    CHECK(check_read(features,
        "CAAGTTGATAACGGACTAGCCCATATAAGAAACTTGGGTTAGTCAGTGAGGTGCTGACCCAGATCGGAAGAGCGTCGTGTAGGGAAAGAG\n", 54) == 0,
        "first actual read");
    CHECK(check_read(features,
        "CAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAACTGTGTCATGGCCTCAAATGACCCATATAAGAAATGCC\n", 70) == 0,
        "second actual read");

    char query[] = "TGACCCA\n";
    int distance = -1;
    CHECK(simpleCorrectFeature(query, features, 0, 1, &distance, NULL) == 1 && distance == 1,
          "fixed-position Hamming search supports maxN zero");
    CHECK(strcmp(query, "TGACCCA\n") == 0, "fixed-position search preserves its input");

    char buffer[4 * 5], *alternatives[4];
    char exact[] = "ACGT", one_n[] = "ACNT", two_n[] = "ANNT";
    CHECK(checkSequenceAndCorrectForN(exact, alternatives, buffer, 4, 0) == 1,
          "maxN zero retains the original sequence");
    CHECK(alternatives[0] == exact, "original sequence occupies the required first slot");
    CHECK(checkSequenceAndCorrectForN(one_n, alternatives, buffer, 4, 0) == 0,
          "maxN zero rejects N before indexing its position");
    CHECK(checkSequenceAndCorrectForN(two_n, alternatives, buffer, 4, 1) == 0,
          "too many Ns reject before writing past the index array");
    CHECK(checkSequenceAndCorrectForN(one_n, alternatives, buffer, 4, 1) == 4,
          "allowed N expands to four sequences");
    CHECK(strcmp(alternatives[0], "ACAT") == 0 && strcmp(alternatives[3], "ACTT") == 0,
          "N expansion keeps the original base order");
    CHECK(strcmp(one_n, "ACNT") == 0, "N expansion preserves the input");

    free_feature_arrays(features);
    feature_code_hash.h64 = NULL;
    puts("PASS: zero-N fallback preserves reads and match positions; N bounds checked before writes");
    return 0;
}

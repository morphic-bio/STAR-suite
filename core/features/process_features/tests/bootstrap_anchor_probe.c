/* Differential test driver: emit per-read decisions and learned offset counts.
 * Run with both the preserved and changed PF libraries on identical fixtures. */
#include "barcode_match.h"
#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc != 7) return 2; /* reference, sequences, output, learning, Hamming, maxN */
    barcode_match_init();
    initialize_complement();
    initialize_unit_sizes();
    feature_prehash_max_hamming = atoi(argv[5]);
    feature_arrays *features = read_features_file(argv[1]);
    if (!features) return 3;
    feature_mode_bootstrap_reads = atoi(argv[4]);
    feature_mode_reads_seen = 0;
    feature_mode_bootstrap_done = 0;
    feature_mode_hist = calloc((size_t)features->number_of_features * feature_mode_max_offset,
                               sizeof(*feature_mode_hist));
    feature_mode_offsets = malloc((size_t)features->number_of_features * sizeof(*feature_mode_offsets));
    if (!feature_mode_hist || !feature_mode_offsets) return 4;
    for (int j = 0; j < features->number_of_features; ++j) feature_mode_offsets[j] = -1;
    use_feature_anchor_search = require_feature_anchor_match = 1;
    use_feature_offset_array = 0;
    limit_search = -1;
    FILE *input = fopen(argv[2], "r"), *output = fopen(argv[3], "w");
    if (!input || !output) return 5;
    char sequence[LINE_LENGTH], original[LINE_LENGTH], matched[LINE_LENGTH];
    size_t ordinal = 0;
    while (fgets(sequence, sizeof(sequence), input)) {
        if (!strchr(sequence, '\n')) return 6;
        strcpy(original, sequence);
        uint32_t feature = 0;
        uint16_t position = 0;
        int distance = atoi(argv[5]) + 1;
        matched[0] = '\0';
        process_feature_sequence(sequence, features, atoi(argv[5]), 1, 0, atoi(argv[6]),
            &feature, &distance, matched, &position, NULL);
        if (strcmp(sequence, original)) return 7;
        fprintf(output, "R\t%zu\t%u\t%d\t%u\t%s\n", ordinal++, feature,
                distance, position, feature ? matched : "");
    }
    for (int j = 0; j < features->number_of_features; ++j) {
        if (feature_mode_offsets[j] >= 0)
            fprintf(output, "M\t%d\t%d\n", j + 1, feature_mode_offsets[j]);
        for (int pos = 0; pos < feature_mode_max_offset; ++pos) {
            unsigned count = feature_mode_hist[(size_t)j * feature_mode_max_offset + pos];
            if (count) fprintf(output, "H\t%d\t%d\t%u\n", j + 1, pos, count);
        }
    }
    fprintf(output, "END\t%zu\t%d\n", ordinal, feature_mode_bootstrap_done);
    fclose(input);
    if (fclose(output)) return 8;
    free(feature_mode_hist); feature_mode_hist = NULL;
    free(feature_mode_offsets); feature_mode_offsets = NULL;
    feature_mode_search_offsets_reset();
    free_feature_arrays(features);
    return 0;
}

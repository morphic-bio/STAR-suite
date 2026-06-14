#include "../include/io.h"
#include <sys/stat.h>
#include <limits.h>

static void parse_feature_pattern(const char *pattern, int *offset, const char **prefix_start, int *prefix_len, const char **suffix_start, int *suffix_len) {
    if (offset) *offset = -1;
    if (prefix_start) *prefix_start = NULL;
    if (prefix_len) *prefix_len = 0;
    if (suffix_start) *suffix_start = NULL;
    if (suffix_len) *suffix_len = 0;
    if (!pattern) {
        return;
    }
    const char *marker = strstr(pattern, "(BC)");
    if (!marker) {
        return;
    }
    const char *suffix = marker + 4;
    const int prefix_length = (int)(marker - pattern);
    const int suffix_length = (int)strlen(suffix);
    if (offset) *offset = prefix_length;
    if (prefix_start) *prefix_start = pattern;
    if (prefix_len) *prefix_len = prefix_length;
    if (suffix_start) *suffix_start = suffix;
    if (suffix_len) *suffix_len = suffix_length;
}

static const char *feature_name_from_csv_line(char *line, int nameIndex, const char *filename, int line_number) {
    char *fields[LINE_LENGTH];
    line[strcspn(line, "\r\n")] = '\0';
    int nFields = split_line(line, fields, ",");
    if (nameIndex >= nFields) {
        fprintf(stderr, "Error: Invalid feature line %s:%d - too few fields %d\n",
                filename, line_number, nFields);
        exit(EXIT_FAILURE);
    }
    return fields[nameIndex];
}

static int feature_name_seen_or_add(khash_t(strptr) *seen_names,
                                    const char *feature_name,
                                    const char *filename,
                                    int line_number,
                                    int warn_on_duplicate) {
    khint_t k = kh_get(strptr, seen_names, feature_name);
    if (k != kh_end(seen_names)) {
        if (warn_on_duplicate) {
            fprintf(stderr,
                    "Warning: duplicate feature name '%s' in %s at line %d; "
                    "ignoring later definition and keeping the first\n",
                    feature_name, filename, line_number);
        }
        return 1;
    }

    char *feature_name_copy = strdup(feature_name);
    if (!feature_name_copy) {
        fprintf(stderr, "Failed to allocate memory for feature-name deduplication\n");
        exit(EXIT_FAILURE);
    }
    int ret = 0;
    k = kh_put(strptr, seen_names, feature_name_copy, &ret);
    if (ret < 0) {
        free(feature_name_copy);
        fprintf(stderr, "Failed to allocate feature-name deduplication hash entry\n");
        exit(EXIT_FAILURE);
    }
    if (ret == 0) {
        free(feature_name_copy);
        return 1;
    }
    kh_val(seen_names, k) = NULL;
    return 0;
}

/* Packed payload for cumulative feature prehash entries:
 * bits [0..28]: feature id (1-based)
 * bits [29..30]: best hamming distance (0..2)
 * bit  [31]: ambiguous at best distance
 */
#define FEATURE_PREHASH_AMBIG_MASK   0x80000000u
#define FEATURE_PREHASH_DIST_SHIFT   29u
#define FEATURE_PREHASH_DIST_MASK    0x60000000u
#define FEATURE_PREHASH_FEATURE_MASK 0x1FFFFFFFu

static inline uint32_t feature_prehash_pack(uint32_t feature_id, uint32_t distance, int ambiguous) {
    return (feature_id & FEATURE_PREHASH_FEATURE_MASK) |
           ((distance & 0x3u) << FEATURE_PREHASH_DIST_SHIFT) |
           (ambiguous ? FEATURE_PREHASH_AMBIG_MASK : 0u);
}

static inline uint32_t feature_prehash_feature(uint32_t payload) {
    return payload & FEATURE_PREHASH_FEATURE_MASK;
}

static inline uint32_t feature_prehash_distance(uint32_t payload) {
    return (payload & FEATURE_PREHASH_DIST_MASK) >> FEATURE_PREHASH_DIST_SHIFT;
}

static inline int feature_prehash_ambiguous(uint32_t payload) {
    return (payload & FEATURE_PREHASH_AMBIG_MASK) != 0u;
}

static int prehash_insert_64(khash_t(u64u32) *hash, uint64_t key,
                             uint32_t feature_index_1based, uint32_t distance,
                             unsigned char *no_ambiguity) {
    int ret;
    khint_t k = kh_put(u64u32, hash, key, &ret);
    if (ret < 0) {
        return -1;
    }
    if (ret > 0) {
        kh_val(hash, k) = feature_prehash_pack(feature_index_1based, distance, 0);
        return 0;
    }
    uint32_t existing = kh_val(hash, k);
    uint32_t ex_feat = feature_prehash_feature(existing);
    uint32_t ex_dist = feature_prehash_distance(existing);
    int ex_amb = feature_prehash_ambiguous(existing);

    if (distance < ex_dist) {
        kh_val(hash, k) = feature_prehash_pack(feature_index_1based, distance, 0);
        return 0;
    }
    if (distance > ex_dist || ex_feat == feature_index_1based) {
        return 0;
    }
    if (ex_amb) {
        if (no_ambiguity) no_ambiguity[feature_index_1based - 1] = 0;
        return 0;
    }
    kh_val(hash, k) = feature_prehash_pack(ex_feat, distance, 1);
    if (no_ambiguity) {
        if (ex_feat > 0) no_ambiguity[ex_feat - 1] = 0;
        no_ambiguity[feature_index_1based - 1] = 0;
    }
    return 0;
}

static int prehash_insert_128(khash_t(seq128u32) *hash, seq128_t key,
                              uint32_t feature_index_1based, uint32_t distance,
                              unsigned char *no_ambiguity) {
    int ret;
    khint_t k = kh_put(seq128u32, hash, key, &ret);
    if (ret < 0) {
        return -1;
    }
    if (ret > 0) {
        kh_val(hash, k) = feature_prehash_pack(feature_index_1based, distance, 0);
        return 0;
    }
    uint32_t existing = kh_val(hash, k);
    uint32_t ex_feat = feature_prehash_feature(existing);
    uint32_t ex_dist = feature_prehash_distance(existing);
    int ex_amb = feature_prehash_ambiguous(existing);

    if (distance < ex_dist) {
        kh_val(hash, k) = feature_prehash_pack(feature_index_1based, distance, 0);
        return 0;
    }
    if (distance > ex_dist || ex_feat == feature_index_1based) {
        return 0;
    }
    if (ex_amb) {
        if (no_ambiguity) no_ambiguity[feature_index_1based - 1] = 0;
        return 0;
    }
    kh_val(hash, k) = feature_prehash_pack(ex_feat, distance, 1);
    if (no_ambiguity) {
        if (ex_feat > 0) no_ambiguity[ex_feat - 1] = 0;
        no_ambiguity[feature_index_1based - 1] = 0;
    }
    return 0;
}

static int prehash_insert(seq_hash_t *sh, seq_key_mode_t mode,
                          const char *seq, int slen,
                          uint32_t feature_index_1based, uint32_t distance,
                          unsigned char *no_ambiguity) {
    if (mode == SEQ_KEY_64) {
        uint64_t key = seq_encode_64_fixed(seq, slen);
        return prehash_insert_64(sh->h64, key, feature_index_1based, distance, no_ambiguity);
    }
    seq128_t key = seq_encode_128_fixed(seq, slen);
    return prehash_insert_128(sh->h128, key, feature_index_1based, distance, no_ambiguity);
}

static unsigned long long estimate_prehash_bytes(uint64_t entries, seq_key_mode_t mode) {
    const size_t key_size = (mode == SEQ_KEY_64) ? 8 : 16;
    const unsigned long long bytes_per_entry = (unsigned long long)key_size + 4ULL + 1ULL;
    return entries * bytes_per_entry * 3ULL / 2ULL;
}

static void prehash_tier_cleanup(seq_hash_t *hashp, seq_key_mode_t mode, unsigned char **ambiguity) {
    seq_hash_destroy(hashp, mode);
    if (ambiguity && *ambiguity) {
        free(*ambiguity);
        *ambiguity = NULL;
    }
}

static void build_feature_hamming_variant_hashes(feature_arrays *features) {
    if (!features) {
        return;
    }
    features->feature_hamming_le1_enabled = 0;
    features->feature_hamming_le2_enabled = 0;
    memset(&features->feature_hamming_le1_hash, 0, sizeof(seq_hash_t));
    memset(&features->feature_hamming_le2_hash, 0, sizeof(seq_hash_t));
    features->feature_no_ambiguity_le1 = NULL;
    features->feature_no_ambiguity_le2 = NULL;

    if (feature_prehash_max_hamming <= 0) {
        fprintf(stderr, "Feature prehash disabled (max_hamming=%d)\n", feature_prehash_max_hamming);
        return;
    }
    if (features->number_of_mismatched_features > 0) {
        fprintf(stderr, "Feature prehash disabled due to %d mismatched-length features\n", features->number_of_mismatched_features);
        return;
    }
    if (features->common_length <= 0) {
        fprintf(stderr, "Feature prehash disabled (common_length=%d)\n", features->common_length);
        return;
    }

    uint64_t n_common = 0;
    for (int i = 0; i < features->number_of_features; i++) {
        if ((int)features->feature_lengths[i] == features->common_length) {
            n_common++;
        }
    }
    if (n_common == 0) {
        fprintf(stderr, "Feature prehash disabled (no common-length features)\n");
        return;
    }

    const uint64_t L = (uint64_t)features->common_length;
    const uint64_t est_le1 = n_common * (1ULL + 3ULL * L);
    const uint64_t est_le2 = n_common * (1ULL + 3ULL * L + 9ULL * (L * (L - 1ULL) / 2ULL));
    const uint64_t entry_budget = (feature_prehash_max_entries > 0ULL) ? feature_prehash_max_entries : ULLONG_MAX;
    const seq_key_mode_t mode = features->code_hash_mode;

    const unsigned long long est_le1_bytes = estimate_prehash_bytes(est_le1, mode);
    const unsigned long long est_le2_bytes = estimate_prehash_bytes(est_le2, mode);
    const unsigned long long mem_budget = feature_prehash_memory_budget;

    /* Entry budget check (existing) */
    if (est_le1 > entry_budget) {
        fprintf(stderr,
                "PREHASH_BUDGET: tier=d1 entries=%llu est_bytes=%lluMB budget=%lluGB decision=SKIP (entry budget)\n",
                (unsigned long long)est_le1,
                (unsigned long long)(est_le1_bytes / (1024ULL*1024ULL)),
                mem_budget ? (unsigned long long)(mem_budget / (1024ULL*1024ULL*1024ULL)) : 0ULL);
        return;
    }

    /* Memory budget check for d1 */
    if (mem_budget > 0 && est_le1_bytes > mem_budget) {
        fprintf(stderr,
                "PREHASH_BUDGET: tier=d1 entries=%llu est_bytes=%lluMB budget=%lluGB decision=SKIP\n",
                (unsigned long long)est_le1,
                (unsigned long long)(est_le1_bytes / (1024ULL*1024ULL)),
                (unsigned long long)(mem_budget / (1024ULL*1024ULL*1024ULL)));
        return;
    }

    fprintf(stderr,
            "PREHASH_BUDGET: tier=d1 entries=%llu est_bytes=%lluMB budget=%lluGB decision=BUILD\n",
            (unsigned long long)est_le1,
            (unsigned long long)(est_le1_bytes / (1024ULL*1024ULL)),
            mem_budget ? (unsigned long long)(mem_budget / (1024ULL*1024ULL*1024ULL)) : 0ULL);

    seq_hash_init(&features->feature_hamming_le1_hash, mode);
    if (mode == SEQ_KEY_64) {
        if (!features->feature_hamming_le1_hash.h64) {
            fprintf(stderr, "WARNING: failed to allocate feature_hamming_le1_hash, skipping prehash\n");
            return;
        }
        kh_resize(u64u32, features->feature_hamming_le1_hash.h64,
                  (khint_t)(est_le1 * 13ULL / 10ULL + 1024ULL));
    } else {
        if (!features->feature_hamming_le1_hash.h128) {
            fprintf(stderr, "WARNING: failed to allocate feature_hamming_le1_hash, skipping prehash\n");
            return;
        }
        kh_resize(seq128u32, features->feature_hamming_le1_hash.h128,
                  (khint_t)(est_le1 * 13ULL / 10ULL + 1024ULL));
    }

    features->feature_no_ambiguity_le1 = malloc((size_t)features->number_of_features);
    if (!features->feature_no_ambiguity_le1) {
        fprintf(stderr, "WARNING: failed to allocate feature_no_ambiguity_le1, skipping prehash\n");
        prehash_tier_cleanup(&features->feature_hamming_le1_hash, mode, NULL);
        return;
    }
    memset(features->feature_no_ambiguity_le1, 1, (size_t)features->number_of_features);

    static const char bases[] = "ACGT";
    char variant[MAX_FEATURE_SEQUENCE_LENGTH + 1];
    int oom = 0;

    for (int i = 0; i < features->number_of_features && !oom; i++) {
        if ((int)features->feature_lengths[i] != features->common_length) continue;
        const char *seq = features->feature_sequences[i];
        if (prehash_insert(&features->feature_hamming_le1_hash, mode, seq, features->common_length,
                           (uint32_t)(i + 1), 0, features->feature_no_ambiguity_le1) < 0) { oom = 1; break; }
        memcpy(variant, seq, (size_t)features->common_length);
        variant[features->common_length] = '\0';
        for (int p = 0; p < features->common_length && !oom; p++) {
            const char orig = seq[p];
            for (int b = 0; b < 4; b++) {
                if (bases[b] == orig) continue;
                variant[p] = bases[b];
                if (prehash_insert(&features->feature_hamming_le1_hash, mode, variant, features->common_length,
                                   (uint32_t)(i + 1), 1, features->feature_no_ambiguity_le1) < 0) { oom = 1; break; }
            }
            variant[p] = orig;
        }
    }

    if (oom) {
        fprintf(stderr, "WARNING: OOM during d1 prehash build, cleaning up and skipping prehash\n");
        prehash_tier_cleanup(&features->feature_hamming_le1_hash, mode, &features->feature_no_ambiguity_le1);
        return;
    }
    features->feature_hamming_le1_enabled = 1;

    if (feature_prehash_max_hamming < 2) {
        fprintf(stderr,
                "Feature prehash enabled: <=1 only (entry_budget=%llu)\n",
                (unsigned long long)entry_budget);
        return;
    }

    /* Memory budget check for d2 */
    if (mem_budget > 0 && est_le2_bytes > mem_budget) {
        fprintf(stderr,
                "PREHASH_BUDGET: tier=d2 entries=%llu est_bytes=%lluMB budget=%lluGB decision=SKIP\n",
                (unsigned long long)est_le2,
                (unsigned long long)(est_le2_bytes / (1024ULL*1024ULL)),
                (unsigned long long)(mem_budget / (1024ULL*1024ULL*1024ULL)));
        fprintf(stderr, "Feature prehash enabled: <=1 only (d2 skipped due to memory budget)\n");
        return;
    }

    /* Entry budget check for d2 */
    if (est_le2 > entry_budget) {
        fprintf(stderr,
                "PREHASH_BUDGET: tier=d2 entries=%llu est_bytes=%lluMB budget=%lluGB decision=SKIP (entry budget)\n",
                (unsigned long long)est_le2,
                (unsigned long long)(est_le2_bytes / (1024ULL*1024ULL)),
                mem_budget ? (unsigned long long)(mem_budget / (1024ULL*1024ULL*1024ULL)) : 0ULL);
        fprintf(stderr, "Feature prehash enabled: <=1 only (d2 skipped due to entry budget)\n");
        return;
    }

    fprintf(stderr,
            "PREHASH_BUDGET: tier=d2 entries=%llu est_bytes=%lluMB budget=%lluGB decision=BUILD\n",
            (unsigned long long)est_le2,
            (unsigned long long)(est_le2_bytes / (1024ULL*1024ULL)),
            mem_budget ? (unsigned long long)(mem_budget / (1024ULL*1024ULL*1024ULL)) : 0ULL);

    seq_hash_init(&features->feature_hamming_le2_hash, mode);
    if (mode == SEQ_KEY_64) {
        if (!features->feature_hamming_le2_hash.h64) {
            fprintf(stderr, "WARNING: failed to allocate feature_hamming_le2_hash, d2 skipped\n");
            return;
        }
        kh_resize(u64u32, features->feature_hamming_le2_hash.h64,
                  (khint_t)(est_le2 * 13ULL / 10ULL + 1024ULL));
    } else {
        if (!features->feature_hamming_le2_hash.h128) {
            fprintf(stderr, "WARNING: failed to allocate feature_hamming_le2_hash, d2 skipped\n");
            return;
        }
        kh_resize(seq128u32, features->feature_hamming_le2_hash.h128,
                  (khint_t)(est_le2 * 13ULL / 10ULL + 1024ULL));
    }

    features->feature_no_ambiguity_le2 = malloc((size_t)features->number_of_features);
    if (!features->feature_no_ambiguity_le2) {
        fprintf(stderr, "WARNING: failed to allocate feature_no_ambiguity_le2, d2 skipped\n");
        prehash_tier_cleanup(&features->feature_hamming_le2_hash, mode, NULL);
        return;
    }
    memset(features->feature_no_ambiguity_le2, 1, (size_t)features->number_of_features);

    oom = 0;
    for (int i = 0; i < features->number_of_features && !oom; i++) {
        if ((int)features->feature_lengths[i] != features->common_length) continue;
        const char *seq = features->feature_sequences[i];
        if (prehash_insert(&features->feature_hamming_le2_hash, mode, seq, features->common_length,
                           (uint32_t)(i + 1), 0, features->feature_no_ambiguity_le2) < 0) { oom = 1; break; }
        memcpy(variant, seq, (size_t)features->common_length);
        variant[features->common_length] = '\0';
        for (int p1 = 0; p1 < features->common_length && !oom; p1++) {
            const char orig1 = seq[p1];
            for (int b1 = 0; b1 < 4 && !oom; b1++) {
                if (bases[b1] == orig1) continue;
                variant[p1] = bases[b1];
                if (prehash_insert(&features->feature_hamming_le2_hash, mode, variant, features->common_length,
                                   (uint32_t)(i + 1), 1, features->feature_no_ambiguity_le2) < 0) { oom = 1; break; }
            }
            variant[p1] = orig1;
            for (int p2 = p1 + 1; p2 < features->common_length && !oom; p2++) {
                const char orig2 = seq[p2];
                for (int b1 = 0; b1 < 4 && !oom; b1++) {
                    if (bases[b1] == orig1) continue;
                    variant[p1] = bases[b1];
                    for (int b2 = 0; b2 < 4 && !oom; b2++) {
                        if (bases[b2] == orig2) continue;
                        variant[p2] = bases[b2];
                        if (prehash_insert(&features->feature_hamming_le2_hash, mode, variant, features->common_length,
                                           (uint32_t)(i + 1), 2, features->feature_no_ambiguity_le2) < 0) { oom = 1; break; }
                    }
                }
                variant[p2] = orig2;
            }
            variant[p1] = orig1;
        }
    }

    if (oom) {
        fprintf(stderr, "WARNING: OOM during d2 prehash build, cleaning up d2 tier (d1 still active)\n");
        prehash_tier_cleanup(&features->feature_hamming_le2_hash, mode, &features->feature_no_ambiguity_le2);
        return;
    }

    features->feature_hamming_le2_enabled = 1;
    fprintf(stderr,
            "Feature prehash enabled: <=1+<=2 (entry_budget=%llu)\n",
            (unsigned long long)entry_budget);
}

feature_arrays* read_features_file(const char* filename) {
    //expext a comma separated file with column names at least one with name and sequence fields
    int seq_size=0;
    int name_size=0;
    int id_size=0;
    int type_size=0;
    int code_size=0;
    int anchor_size=0;
    int suffix_anchor_size=0;
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open tags file");
        exit(EXIT_FAILURE);
    }
    char line[LINE_LENGTH];
    clear_feature_lookup_hashes();
    int count = 0;
    //skip the header
    //count the lines and check that the sequences are valid
    int maxFeatureLength=0;
    int seqIndex=-1;
    int nameIndex=-1;
    int patternIndex=-1;
    int idIndex=-1;
    int featureTypeIndex=-1;
    if (!fgets(line, LINE_LENGTH, file)) {
        perror("Failed to read tags header");
        exit(EXIT_FAILURE);
    }
    find_feature_ref_fields(line, &nameIndex, &seqIndex, &patternIndex, &idIndex, &featureTypeIndex);
    
    khash_t(u32u32)* length_counts = kh_init(u32u32);
    khash_t(strptr)* seen_feature_names = kh_init(strptr);
    int line_number = 1;

    while (fgets(line, LINE_LENGTH, file) != NULL) {
        line_number++;
        char name_line[LINE_LENGTH];
        strncpy(name_line, line, sizeof(name_line) - 1);
        name_line[sizeof(name_line) - 1] = '\0';
        const char *feature_name = feature_name_from_csv_line(name_line, nameIndex, filename, line_number);
        if (feature_name_seen_or_add(seen_feature_names, feature_name, filename, line_number, 1)) {
            continue;
        }
        int length = get_feature_line_sizes(line, nameIndex, seqIndex, patternIndex, idIndex, featureTypeIndex,
                                            &name_size, &seq_size, &code_size, &anchor_size, &suffix_anchor_size,
                                            &id_size, &type_size, &maxFeatureLength);
        if (length > 0) {
            khint_t k = kh_get(u32u32, length_counts, length);
            uint32_t current_count = 1;
            if (k != kh_end(length_counts)) {
                current_count = kh_val(length_counts, k) + 1;
                kh_val(length_counts, k) = current_count;
            } else {
                int ret;
                khint_t kh = kh_put(u32u32, length_counts, length, &ret);
                kh_val(length_counts, kh) = current_count;
            }
        }
        count++;
    }
    free_strptr_hash(seen_feature_names);

    int most_common_length = 0;
    int max_count = 0;
    khint_t k;
    for (k = kh_begin(length_counts); k != kh_end(length_counts); ++k) {
        if (!kh_exist(length_counts, k)) continue;
        int length = kh_key(length_counts, k);
        int current_count = kh_val(length_counts, k);
        if (current_count > max_count) {
            max_count = current_count;
            most_common_length = length;
        }
    }
    kh_destroy(u32u32, length_counts);

    fprintf(stderr, "Read %d unique tags with max length %d and most common length %d\n", count, maxFeatureLength, most_common_length);

    feature_arrays *myfeatures = allocate_feature_arrays(name_size, seq_size, code_size, anchor_size, suffix_anchor_size, id_size, type_size, count, maxFeatureLength);
    myfeatures->common_length = most_common_length;
    myfeatures->source_csv_path = strdup(filename);
    myfeatures->source_csv_fingerprint = malloc(17);
    if (myfeatures->source_csv_fingerprint) {
        if (pf_file_fingerprint(filename, myfeatures->source_csv_fingerprint, 17) != 0) {
            strcpy(myfeatures->source_csv_fingerprint, "unknown");
        }
    }

    /* Select integer-key mode based on sequence length and init per-instance hash. */
    if (most_common_length <= MAX_FIXED_KEY_SEQ_LENGTH_64) {
        myfeatures->code_hash_mode = SEQ_KEY_64;
    } else {
        myfeatures->code_hash_mode = SEQ_KEY_128;
    }
    myfeatures->code_hash_fixed_length = most_common_length;
    seq_hash_init(&myfeatures->feature_code_hash, myfeatures->code_hash_mode);

    /* Set globals as non-owning aliases for legacy code */
    feature_code_hash = myfeatures->feature_code_hash;
    feature_code_hash_mode = myfeatures->code_hash_mode;
    feature_code_hash_fixed_length = myfeatures->code_hash_fixed_length;
    //rewind the file and read the sequences
    fseek(file, 0, SEEK_SET);
    if (!fgets(line, LINE_LENGTH, file)) {
        perror("Failed to headers file");
        exit(EXIT_FAILURE);
    }
    khash_t(strptr)* emitted_feature_names = kh_init(strptr);
    count=0;
    line_number = 1;
    while (fgets(line, LINE_LENGTH, file) != NULL) {
        line_number++;
        char name_line[LINE_LENGTH];
        strncpy(name_line, line, sizeof(name_line) - 1);
        name_line[sizeof(name_line) - 1] = '\0';
        const char *feature_name = feature_name_from_csv_line(name_line, nameIndex, filename, line_number);
        if (feature_name_seen_or_add(emitted_feature_names, feature_name, filename, line_number, 0)) {
            continue;
        }
        process_feature_line(line, nameIndex, seqIndex, patternIndex, idIndex, featureTypeIndex, myfeatures, count);
        count++;
    }
    free_strptr_hash(emitted_feature_names);
    fprintf(stderr, "Read %d unique tags\n", count);
    fclose(file);

    int mismatched_count = 0;
    for (int i = 0; i < myfeatures->number_of_features; i++) {
        if (myfeatures->feature_lengths[i] != myfeatures->common_length) {
            myfeatures->mismatched_feature_indices[mismatched_count++] = i;
        }
    }
    myfeatures->number_of_mismatched_features = mismatched_count;
    fprintf(stderr, "Found %d features with length different from common length %d\n", mismatched_count, myfeatures->common_length);
    build_feature_hamming_variant_hashes(myfeatures);

    return myfeatures;
}
int get_feature_line_sizes(char *line, int nameIndex, int seqIndex, int patternIndex, int idIndex, int featureTypeIndex,
                           int *name_size, int *seq_size, int *code_size, int *anchor_size, int *suffix_anchor_size,
                           int *id_size, int *type_size, int *maxFeatureLength) {
    line[strcspn(line, "\r\n")] = '\0';
    char *fields[LINE_LENGTH];
    int nFields = split_line(line, fields, ",");
    if (seqIndex >= nFields || nameIndex >= nFields) {
        fprintf(stderr, "Error: Invalid line - too few fields %d \n",nFields);
        exit(EXIT_FAILURE);
    }
    char *tmpSeq = fields[seqIndex];
            // Remove possible newline character from the sequence
    if (!check_sequence(tmpSeq, strlen(tmpSeq))){
        //fprintf(stderr, "Error: Invalid sequence %s\n", tmpSeq);
        exit(EXIT_FAILURE);
    }
    const int string_length = strlen(tmpSeq);
    *seq_size += string_length + 1; 
    *code_size += string_length / 4;
    if (string_length % 4){
        (*code_size)++;
    } 
    if (string_length > *maxFeatureLength){
        *maxFeatureLength = string_length;
    }
    *name_size += strlen(fields[nameIndex]) + 1;
    if (idIndex >= 0 && idIndex < nFields && id_size) {
        *id_size += strlen(fields[idIndex]) + 1;
    }
    if (featureTypeIndex >= 0 && featureTypeIndex < nFields && type_size) {
        *type_size += strlen(fields[featureTypeIndex]) + 1;
    }
    if (anchor_size || suffix_anchor_size) {
        int prefix_anchor_len = 0;
        int suffix_anchor_len = 0;
        if (patternIndex >= 0 && patternIndex < nFields) {
            parse_feature_pattern(fields[patternIndex], NULL, NULL, &prefix_anchor_len, NULL, &suffix_anchor_len);
        }
        if (anchor_size) {
            *anchor_size += prefix_anchor_len + 1;
        }
        if (suffix_anchor_size) {
            *suffix_anchor_size += suffix_anchor_len + 1;
        }
    }
    return string_length;
}
void process_feature_line(char *line, int nameIndex, int seqIndex, int patternIndex, int idIndex, int featureTypeIndex, feature_arrays *myfeatures, int count) {
    // Split the line by spaces and read the 3rd and 6th columns
    char *fields[LINE_LENGTH];
    line[strcspn(line, "\r\n")] = 0;
    int nFields = split_line(line, fields, ",");
    if (seqIndex >= nFields || nameIndex >= nFields) {
        fprintf(stderr, "Error: Invalid line - two few fields %d \n",nFields);
        exit(EXIT_FAILURE);
    }
    //copy the name and sequence into the feature arrays
    char *tmpName = fields[nameIndex];
    strcpy(myfeatures->feature_names[count], fields[nameIndex]);
    if (myfeatures->feature_ids && idIndex >= 0 && idIndex < nFields) {
        strcpy(myfeatures->feature_ids[count], fields[idIndex]);
        if (count + 1 < myfeatures->number_of_features) {
            myfeatures->feature_ids[count + 1] = myfeatures->feature_ids[count] + strlen(fields[idIndex]) + 1;
        }
    } else if (myfeatures->feature_ids) {
        strcpy(myfeatures->feature_ids[count], fields[nameIndex]);
        if (count + 1 < myfeatures->number_of_features) {
            myfeatures->feature_ids[count + 1] = myfeatures->feature_ids[count] + strlen(tmpName) + 1;
        }
    }
    if (myfeatures->feature_types && featureTypeIndex >= 0 && featureTypeIndex < nFields) {
        strcpy(myfeatures->feature_types[count], fields[featureTypeIndex]);
        if (count + 1 < myfeatures->number_of_features) {
            myfeatures->feature_types[count + 1] = myfeatures->feature_types[count] + strlen(fields[featureTypeIndex]) + 1;
        }
    }
    if (count + 1 < myfeatures->number_of_features) {
        myfeatures->feature_names[count + 1] = myfeatures->feature_names[count] + strlen(tmpName) + 1;
    }
    char *tmpSeq = fields[seqIndex];
    strcpy(myfeatures->feature_sequences[count], tmpSeq);
    myfeatures->feature_lengths[count] = strlen(tmpSeq);
    myfeatures->feature_code_lengths[count] = string2code(tmpSeq, strlen(tmpSeq), myfeatures->feature_codes[count]);
    int offset = -1;
    const char *prefix_anchor = NULL;
    const char *suffix_anchor = NULL;
    int prefix_anchor_len = 0;
    int suffix_anchor_len = 0;
    if (patternIndex >= 0 && patternIndex < nFields) {
        parse_feature_pattern(fields[patternIndex], &offset, &prefix_anchor, &prefix_anchor_len, &suffix_anchor, &suffix_anchor_len);
    }
    if (myfeatures->feature_offsets) {
        // feature_offsets is 0-based; index i corresponds to feature index (i+1)
        myfeatures->feature_offsets[count] = offset;
    }
    if (myfeatures->feature_anchors) {
        if (prefix_anchor_len > 0 && prefix_anchor) {
            memcpy(myfeatures->feature_anchors[count], prefix_anchor, (size_t)prefix_anchor_len);
        }
        myfeatures->feature_anchors[count][prefix_anchor_len] = '\0';
        myfeatures->feature_anchor_lengths[count] = (unsigned int)prefix_anchor_len;
        if (count + 1 < myfeatures->number_of_features) {
            myfeatures->feature_anchors[count + 1] = myfeatures->feature_anchors[count] + prefix_anchor_len + 1;
        }
    }
    if (myfeatures->feature_suffix_anchors) {
        if (suffix_anchor_len > 0 && suffix_anchor) {
            memcpy(myfeatures->feature_suffix_anchors[count], suffix_anchor, (size_t)suffix_anchor_len);
        }
        myfeatures->feature_suffix_anchors[count][suffix_anchor_len] = '\0';
        myfeatures->feature_suffix_anchor_lengths[count] = (unsigned int)suffix_anchor_len;
        if (count + 1 < myfeatures->number_of_features) {
            myfeatures->feature_suffix_anchors[count + 1] = myfeatures->feature_suffix_anchors[count] + suffix_anchor_len + 1;
        }
    }
    if (myfeatures->feature_lengths[count] == myfeatures->common_length) {
        const char *fseq = myfeatures->feature_sequences[count];
        const int flen = (int)myfeatures->feature_lengths[count];
        if (feature_code_hash_mode == SEQ_KEY_64) {
            uint64_t ikey = seq_encode_64_fixed(fseq, flen);
            seq_hash_put_64(&feature_code_hash, ikey, (uint32_t)(count + 1));
        } else {
            seq128_t ikey = seq_encode_128_fixed(fseq, flen);
            seq_hash_put_128(&feature_code_hash, ikey, (uint32_t)(count + 1));
        }
    }
    if (count + 1 < myfeatures->number_of_features) {
        myfeatures->feature_sequences[count + 1] = myfeatures->feature_sequences[count] + strlen(tmpSeq) + 1;
        myfeatures->feature_codes[count + 1] = myfeatures->feature_codes[count] + myfeatures->feature_code_lengths[count];
    }
}
feature_arrays* allocate_feature_arrays(int name_size, int seq_size, int code_size, int anchor_size, int suffix_anchor_size, int id_size, int type_size, int count, int maxFeatureLength) {
        feature_arrays *myfeatures = malloc(sizeof(feature_arrays));
        if (myfeatures == NULL) {
            fprintf(stderr, "Failed to allocate memory for feature arrays\n");
            exit(EXIT_FAILURE);
        }
        memset(myfeatures, 0, sizeof(feature_arrays));
        myfeatures->max_length = maxFeatureLength;
        myfeatures->feature_names_storage = malloc(name_size);
        myfeatures->feature_sequences_storage = malloc(seq_size);
        myfeatures->feature_codes_storage = malloc(code_size);
        myfeatures->feature_names = malloc(count * sizeof(char*));
        myfeatures->feature_ids = NULL;
        myfeatures->feature_ids_storage = NULL;
        myfeatures->feature_types = NULL;
        myfeatures->feature_types_storage = NULL;
        if (id_size > 0) {
            myfeatures->feature_ids_storage = malloc(id_size);
            myfeatures->feature_ids = malloc(count * sizeof(char*));
        }
        if (type_size > 0) {
            myfeatures->feature_types_storage = malloc(type_size);
            myfeatures->feature_types = malloc(count * sizeof(char*));
        }
        myfeatures->feature_lengths = malloc(count * sizeof(unsigned int));
        myfeatures->feature_code_lengths = malloc(count * sizeof(unsigned char));
        myfeatures->feature_codes = malloc(count * sizeof(unsigned char*));
        myfeatures->feature_sequences = malloc(count * sizeof(char*));
        myfeatures->feature_anchors_storage = NULL;
        myfeatures->feature_anchors = NULL;
        myfeatures->feature_anchor_lengths = NULL;
        myfeatures->feature_suffix_anchors_storage = NULL;
        myfeatures->feature_suffix_anchors = NULL;
        myfeatures->feature_suffix_anchor_lengths = NULL;
        if (anchor_size > 0) {
            myfeatures->feature_anchors_storage = malloc(anchor_size);
            myfeatures->feature_anchors = malloc(count * sizeof(char*));
            myfeatures->feature_anchor_lengths = malloc(count * sizeof(unsigned int));
        }
        if (suffix_anchor_size > 0) {
            myfeatures->feature_suffix_anchors_storage = malloc(suffix_anchor_size);
            myfeatures->feature_suffix_anchors = malloc(count * sizeof(char*));
            myfeatures->feature_suffix_anchor_lengths = malloc(count * sizeof(unsigned int));
        }
        myfeatures->feature_offsets = malloc(count * sizeof(int));
        myfeatures->number_of_features = count;
        myfeatures->mismatched_feature_indices = malloc(count * sizeof(int));

        // Check if any of the mallocs failed by checking for NULL pointers
        if (myfeatures->feature_names_storage == NULL || myfeatures->feature_sequences_storage == NULL || myfeatures->feature_codes_storage == NULL || myfeatures->feature_names == NULL || myfeatures->feature_lengths == NULL || myfeatures->feature_code_lengths == NULL || myfeatures->feature_codes == NULL || myfeatures->mismatched_feature_indices == NULL || myfeatures->feature_offsets == NULL || (id_size > 0 && (!myfeatures->feature_ids_storage || !myfeatures->feature_ids)) || (type_size > 0 && (!myfeatures->feature_types_storage || !myfeatures->feature_types)) || (anchor_size > 0 && (!myfeatures->feature_anchors_storage || !myfeatures->feature_anchors || !myfeatures->feature_anchor_lengths)) || (suffix_anchor_size > 0 && (!myfeatures->feature_suffix_anchors_storage || !myfeatures->feature_suffix_anchors || !myfeatures->feature_suffix_anchor_lengths))) {
            fprintf(stderr, "Failed to allocate memory for feature arrays\n");
            exit(EXIT_FAILURE);
        }
        memset(myfeatures->feature_names_storage, 0, name_size);
        memset(myfeatures->feature_sequences_storage, 0, seq_size);
        memset(myfeatures->feature_codes_storage, 0, code_size);
        for (int i = 0; i < count; i++) {
            myfeatures->feature_offsets[i] = -1;
        }
        if (count > 0) {
            myfeatures->feature_names[0] = myfeatures->feature_names_storage;
            myfeatures->feature_sequences[0] = myfeatures->feature_sequences_storage;
            myfeatures->feature_codes[0] = myfeatures->feature_codes_storage;
            if (myfeatures->feature_ids_storage) {
                memset(myfeatures->feature_ids_storage, 0, id_size);
                myfeatures->feature_ids[0] = myfeatures->feature_ids_storage;
            }
            if (myfeatures->feature_types_storage) {
                memset(myfeatures->feature_types_storage, 0, type_size);
                myfeatures->feature_types[0] = myfeatures->feature_types_storage;
            }
            if (myfeatures->feature_anchors_storage) {
                memset(myfeatures->feature_anchors_storage, 0, anchor_size);
                for (int i = 0; i < count; i++) {
                    myfeatures->feature_anchor_lengths[i] = 0;
                }
                myfeatures->feature_anchors[0] = myfeatures->feature_anchors_storage;
            }
            if (myfeatures->feature_suffix_anchors_storage) {
                memset(myfeatures->feature_suffix_anchors_storage, 0, suffix_anchor_size);
                for (int i = 0; i < count; i++) {
                    myfeatures->feature_suffix_anchor_lengths[i] = 0;
                }
                myfeatures->feature_suffix_anchors[0] = myfeatures->feature_suffix_anchors_storage;
            }
        }

        return myfeatures;
    }
void find_feature_ref_fields(char *line, int *nameIndex, int *seqIndex, int *patternIndex,
                             int *idIndex, int *featureTypeIndex) {
    char *fields[LINE_LENGTH];
    line[strcspn(line, "\r\n")] = 0;
    int nFields = split_line(line, fields, ",");
    if (nameIndex) *nameIndex = -1;
    if (seqIndex) *seqIndex = -1;
    if (patternIndex) *patternIndex = -1;
    if (idIndex) *idIndex = -1;
    if (featureTypeIndex) *featureTypeIndex = -1;
    if (nFields < 2) {
        fprintf(stderr, "Error: Invalid header in tags file - there must be at least a name and sequence field \n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < nFields; i++) {
        char normalized[LINE_LENGTH];
        strncpy(normalized, fields[i], sizeof(normalized) - 1);
        normalized[sizeof(normalized) - 1] = '\0';
        for (char *p = normalized; *p; ++p) {
            if (*p >= 'A' && *p <= 'Z') *p = (char)(*p + ('a' - 'A'));
        }
        if (strcmp(normalized, "sequence") == 0) {
            if (seqIndex) *seqIndex = i;
        } else if (strcmp(normalized, "name") == 0) {
            if (nameIndex) *nameIndex = i;
        } else if (strcmp(normalized, "id") == 0) {
            if (idIndex) *idIndex = i;
        } else if (strcmp(normalized, "feature_type") == 0 || strcmp(normalized, "type") == 0) {
            if (featureTypeIndex) *featureTypeIndex = i;
        } else if (patternIndex && strcmp(normalized, "pattern") == 0) {
            *patternIndex = i;
        }
    }
    if (!seqIndex || *seqIndex < 0 || !nameIndex || *nameIndex < 0) {
        fprintf(stderr, "Error: Invalid header in tags file - there must be at least a name and sequence field \n");
        exit(EXIT_FAILURE);
    }
}

void find_name_and_sequence_fields(char *line, int *nameIndex, int *seqIndex, int *patternIndex) {
    find_feature_ref_fields(line, nameIndex, seqIndex, patternIndex, NULL, NULL);
}
int put_fastq_files_string_into_collection(char *fastqFilesString, char **fastq_files, int *nFiles, char *concatenated_fastq) {
    if (!fastqFilesString) {
        return 0;
    }
    char *this_file = concatenated_fastq;
    *nFiles = 1;
    strcpy(this_file, fastqFilesString);
    while (*this_file) {
        if (*this_file == ',') {
            (*nFiles)++;
            *this_file = '\0';
        }
        this_file++;
    }
    this_file = concatenated_fastq;
    for (int i = 0; i < *nFiles; i++) {
        fastq_files[i] = this_file;
        this_file += strlen(this_file) + 1;
    }
    return *nFiles;
}
void check_filecounts(fastq_files_collection *fastq_files) {
    if (!fastq_files->nbarcode_files) {
        fprintf(stderr, "Error: No barcode fastq files\n");
        exit(EXIT_FAILURE);
    }
    if (!fastq_files->nforward_files && !fastq_files->nreverse_files) {
        fprintf(stderr, "Error: No forward or reverse fastq files\n");
        exit(EXIT_FAILURE);
    }
    if (fastq_files->nforward_files && fastq_files->nforward_files != fastq_files->nbarcode_files) {
        fprintf(stderr, "Error: Unequal number of barcode and forward fastq files\n");
        exit(EXIT_FAILURE);
    }
    if (fastq_files->nreverse_files && fastq_files->nreverse_files != fastq_files->nbarcode_files) {
        fprintf(stderr, "Error: Unequal number of barcode and reverse fastq files\n");
        exit(EXIT_FAILURE);
    }
}
char* extract_sample_name(char *filename, char *pattern) {
    const char *samplename=get_basename(filename);
    char *pattern_position = strstr(samplename, pattern);
    if (!pattern_position) {
        fprintf(stderr, "Error: FASTQ file name %s does not contain the sample pattern %s\n", samplename, pattern);
        exit(EXIT_FAILURE);
    }
    *pattern_position = '\0';
    //find last underscore
    char *underscore_position = strrchr(samplename, '_');
    if (!underscore_position) {
        return (char*) samplename;
    }
    if (*(underscore_position+1) == 'L' && isdigit(*(underscore_position+2))){
        *underscore_position = '\0';
    }
    return (char*) samplename;
}
int count_character(char *string, char character) {
    int count = 0;
    for (int i = 0; string[i]; i++) {
        if (string[i] == character) {
            count++;
        }
    }
    return count;
}
int compare_filenames(const void *a, const void *b) {
    return strcmp(*(const char **)a, *(const char **)b);
}
int count_files_with_pattern(const char *directory_path, const char *pattern) {
    DIR *dir = opendir(directory_path);
    if (!dir) {
        perror("Unable to open directory");
        return -1;  // Return -1 to indicate an error
    }

    struct dirent *entry;
    struct stat file_stat;
    char filepath[FILENAME_LENGTH];
    int count = 0;

    // Iterate through the directory entries
    while ((entry = readdir(dir)) != NULL) {
        // Skip the "." and ".." entries
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue;
        }

        // Build the full file path
        snprintf(filepath, FILENAME_LENGTH, "%s/%s", directory_path, entry->d_name);

        // Check if the entry is a regular file (not a directory)
        if (stat(filepath, &file_stat) == 0 && S_ISREG(file_stat.st_mode)) {
            // Check if the filename contains the pattern
            if (strstr(entry->d_name, pattern) != NULL) {
                count++;  // Increment the count if a match is found
            }
        }
    }

    closedir(dir);
    return count;  // Return the count of matching files
}
char **find_files_with_pattern(const char *directory_path, const char *pattern, int *num_files_found) {
    DIR *dir = opendir(directory_path);
    if (!dir) {
        perror("Unable to open directory");
        return NULL;
    }

    struct dirent *entry;
    struct stat file_stat;
    char filepath[FILENAME_LENGTH];
    char **filepaths = NULL;  // Dynamically allocated array to store filepaths
    *num_files_found = 0;

    while ((entry = readdir(dir)) != NULL) {
        // Skip the "." and ".." entries
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue;
        }

        // Build the full file path
        strcpy(filepath, directory_path);
        if (filepath[strlen(filepath) - 1] != '/') {
            strcat(filepath, "/");
        }
        strcpy(filepath+strlen(filepath), entry->d_name);

        // Check if the entry is a regular file (not a directory)
        if (stat(filepath, &file_stat) == 0 && S_ISREG(file_stat.st_mode)) {
            // Check if the filename contains the pattern
            if (strstr(entry->d_name, pattern) != NULL) {
                // Store the matching filepath in the array
                filepaths = realloc(filepaths, (*num_files_found + 1) * sizeof(char *));
                if (!filepaths) {
                    perror("Memory allocation error");
                    exit(EXIT_FAILURE);
                }

                filepaths[*num_files_found] = strdup(filepath);  // Save a copy of the full file path
                if (!filepaths[*num_files_found]) {
                    perror("Memory allocation error");
                    exit(EXIT_FAILURE);
                }

                (*num_files_found)++;
            }
        }
    }

    closedir(dir);
    if (*num_files_found == 0) {
        free(filepaths);
        return NULL;  // Return NULL if no files were found
    }

    // Sort the filepaths alphabetically
    qsort(filepaths, *num_files_found, sizeof(char *), compare_filenames);

    return filepaths;  // Return the array of filepaths
}
void organize_fastq_files_by_directory(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern) {
    fastq_files->barcode_fastq = 0;
    fastq_files->forward_fastq = 0;
    fastq_files->reverse_fastq = 0;
    fastq_files->nbarcode_files = 0;
    fastq_files->nforward_files = 0;
    fastq_files->nreverse_files = 0;
    if (positional_arg_count) {
        //count the files in the first directory
        int barcode_file_exist=count_files_with_pattern(argv[optind], barcode_pattern);
        int forward_file_exist=count_files_with_pattern(argv[optind], forward_pattern);
        int reverse_file_exist=count_files_with_pattern(argv[optind], reverse_pattern);
        int total_barcode_files_found=0;
        for (int i=0;i < positional_arg_count;i++){
            total_barcode_files_found+=count_files_with_pattern(argv[optind+i], barcode_pattern);
        }
        if (!barcode_file_exist) {
            fprintf(stderr, "Error: No barcode fastq files found in directory %s\n", argv[optind]);
            exit(EXIT_FAILURE);
        }       
        if (!forward_file_exist && !reverse_file_exist) {
            fprintf(stderr, "Error: No forward or reverse fastq files found in directory %s\n", argv[optind]);
            exit(EXIT_FAILURE);
        }
        if (forward_file_exist) {
            fastq_files->forward_fastq = calloc(total_barcode_files_found, sizeof(char *));
            if (fastq_files->forward_fastq == NULL) {
                perror("Failed to allocate memory for forward fastq files");
                exit(EXIT_FAILURE);
            }       
        }
        if (reverse_file_exist) {
            fastq_files->reverse_fastq = calloc(total_barcode_files_found, sizeof(char *));
            if (fastq_files->reverse_fastq == NULL) {
                perror("Failed to allocate memory for reverse fastq files");
                exit(EXIT_FAILURE);
            }
        }
        fastq_files->barcode_fastq = calloc(total_barcode_files_found, sizeof(char *));
        fastq_files->sample_sizes=calloc(positional_arg_count,sizeof(int));
        fastq_files->sample_names=malloc(positional_arg_count*sizeof(char*));
        fastq_files->sample_offsets=malloc(positional_arg_count*sizeof(int));
        fastq_files->nsamples=positional_arg_count;
        fastq_files->sorted_index=malloc(positional_arg_count*sizeof(int));
        //check that the memory allocation was successful
        if (fastq_files->barcode_fastq == NULL || (forward_file_exist && fastq_files->forward_fastq == NULL) || (reverse_file_exist && fastq_files->reverse_fastq == NULL) || !fastq_files->sample_sizes || !fastq_files->sample_names || !fastq_files->sample_offsets || !fastq_files->sorted_index) {
            perror("Failed to allocate memory for fastq files");
            exit(EXIT_FAILURE);
        }
        total_barcode_files_found=0;
        for (int i = 0; i < positional_arg_count; i++) {
            int num_barcode_files_found = 0;
            int num_forward_files_found = 0;
            int num_reverse_files_found = 0;
            char *directory = strdup(argv[optind + i]); //this gets modified later
            char **sample_barcode_fastq=0, **sample_forward_fastq=0, **sample_reverse_fastq=0;

            sample_barcode_fastq=find_files_with_pattern(directory,barcode_pattern, &num_barcode_files_found);
            if (forward_file_exist) {
                sample_forward_fastq=find_files_with_pattern(directory,forward_pattern, &num_forward_files_found);
            }
            if (reverse_file_exist) {
                sample_reverse_fastq=find_files_with_pattern(directory,reverse_pattern, &num_reverse_files_found);
            }
            if (!num_barcode_files_found) {
                fprintf(stderr, "Error: No barcode fastq files found in directory %s\n", directory);
                exit(EXIT_FAILURE);
            }
            if (!num_forward_files_found && !num_reverse_files_found) {
                fprintf(stderr, "Error: No forward or reverse fastq files found in directory %s\n", directory);
                exit(EXIT_FAILURE);
            }
            if (num_forward_files_found && num_forward_files_found != num_barcode_files_found) {
                fprintf(stderr, "Error: Unequal number of barcode and forward fastq files in directory %s\n", directory);
                exit(EXIT_FAILURE);
            }
            if (num_reverse_files_found && num_reverse_files_found != num_barcode_files_found) {
                fprintf(stderr, "Error: Unequal number of barcode and reverse fastq files in directory %s\n", directory);
                exit(EXIT_FAILURE);
            }
            fastq_files->sample_sizes[i]=num_barcode_files_found;
            fastq_files->sample_offsets[i]=total_barcode_files_found;
            total_barcode_files_found+=num_barcode_files_found;
            for(int j=0;j<num_barcode_files_found;j++){
                fastq_files->barcode_fastq[fastq_files->nbarcode_files++]=sample_barcode_fastq[j];
                if (num_forward_files_found){
                    fastq_files->forward_fastq[fastq_files->nforward_files++]=sample_forward_fastq[j];
                }
                if (num_reverse_files_found){
                    fastq_files->reverse_fastq[fastq_files->nreverse_files++]=sample_reverse_fastq[j];
                }
            }
            
            fastq_files->sample_names[i] = strdup(get_basename(directory));
            fprintf(stderr, "directory %s Sample name %s\n", directory,fastq_files->sample_names[i]);
            free(directory);
            // Free the temporary arrays (individual strings were transferred, not copied)
            free(sample_barcode_fastq);
            if (sample_forward_fastq) free(sample_forward_fastq);
            if (sample_reverse_fastq) free(sample_reverse_fastq);
        }
    }    
    sort_samples_by_size(fastq_files, fastq_files->sorted_index);
    int max_sample_size=0;
    for (int i=0; i<fastq_files->nsamples; i++){
        if (fastq_files->sample_sizes[i]>max_sample_size){
            max_sample_size=fastq_files->sample_sizes[i];
        }
    }
    fastq_files->max_sample_size=max_sample_size;
}

void organize_fastq_files_by_type(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern, int sample_flag) {
    fastq_files->barcode_fastq = 0;
    fastq_files->forward_fastq = 0;
    fastq_files->reverse_fastq = 0;
    fastq_files->nbarcode_files = 0;
    fastq_files->nforward_files = 0;
    fastq_files->nreverse_files = 0;
    if (positional_arg_count) {
        //allocate memory for the fastq files for now

        //count the size of the positional arguments
        size_t block_size = 0;
        for (int i = 0; i < positional_arg_count; i++) {
            block_size += strlen(argv[optind + i]) + 1;
        }
        fastq_files->concatenated_files = malloc(block_size);
        //check if the memory allocation was successful
        if (fastq_files->concatenated_files == NULL) {
            perror("Failed to allocate memory for concatenated fastq files");
            exit(EXIT_FAILURE);
        }

        //copy the positional arguments to the concatenated fastq files
        memset(fastq_files->concatenated_files, 0, block_size);
        char *this_file = fastq_files->concatenated_files;
        //count the file types and allocate as necessary
        for (int i = 0; i < positional_arg_count; i++) {
            strcpy(this_file, argv[optind + i]);
            char *barcode_pattern_position = strstr(this_file, barcode_pattern);
            char *forward_pattern_position = strstr(this_file, forward_pattern);
            char *reverse_pattern_position = strstr(this_file, reverse_pattern);
            if (!barcode_pattern_position && !forward_pattern_position && !reverse_pattern_position) {
                fprintf(stderr, "Error: FASTQ file name %s does not contain the barcode, forward and reverse patterns\n", this_file);
                exit(EXIT_FAILURE);
            }
            if (barcode_pattern_position && !forward_pattern_position && !reverse_pattern_position) {
                if (!fastq_files->barcode_fastq){
                    fastq_files->barcode_fastq = calloc(positional_arg_count, sizeof(char *));
                    if (fastq_files->barcode_fastq == NULL) {
                        perror("Failed to allocate memory for barcode fastq files");
                        exit(EXIT_FAILURE);
                    }
                }
                fastq_files->barcode_fastq[fastq_files->nbarcode_files++] = this_file;
            } else if (!barcode_pattern_position && forward_pattern_position && !reverse_pattern_position) {
                if (!fastq_files->forward_fastq){
                    fastq_files->forward_fastq = calloc(positional_arg_count, sizeof(char *));
                    if (fastq_files->forward_fastq == NULL) {
                        perror("Failed to allocate memory for forward fastq files");
                        exit(EXIT_FAILURE);
                    }
                }
                fastq_files->forward_fastq[fastq_files->nforward_files++] = this_file;
            } else if (!barcode_pattern_position && !forward_pattern_position && reverse_pattern_position) {
                if (!fastq_files->reverse_fastq){
                    fastq_files->reverse_fastq = calloc(positional_arg_count, sizeof(char *));
                    if (fastq_files->reverse_fastq == NULL) {
                        perror("Failed to allocate memory for reverse fastq files");
                        exit(EXIT_FAILURE);
                    }
                }
                fastq_files->reverse_fastq[fastq_files->nreverse_files++] = this_file;
            } else {
                fprintf(stderr, "Error: FASTQ file name %s contains more than one of the barcode, forward and reverse patterns\n", this_file);
                exit(EXIT_FAILURE);
            }
            this_file += strlen(argv[optind + i]) + 1;
        }
    } else {
        //find number of files in each type
        int barcode_count = 0;
        int forward_count = 0;
        int reverse_count = 0;
        if (barcodeFastqFilesString) {
            barcode_count = count_character(barcodeFastqFilesString, ',') + 1;
            fastq_files->barcode_fastq = calloc(barcode_count, sizeof(char *));
            if (fastq_files->barcode_fastq == NULL) {
                perror("Failed to allocate memory for barcode fastq files");
                exit(EXIT_FAILURE);
            }
        }
        if (forwardFastqFilesString) {
            forward_count = count_character(forwardFastqFilesString, ',') + 1;
            fastq_files->forward_fastq = calloc(forward_count, sizeof(char *));
            if (fastq_files->forward_fastq == NULL) {
                perror("Failed to allocate memory for forward fastq files");
                exit(EXIT_FAILURE);
            }
        }
        if (reverseFastqFilesString) {
            reverse_count = count_character(reverseFastqFilesString, ',') + 1;
            fastq_files->reverse_fastq = calloc(reverse_count, sizeof(char *));
            if (fastq_files->reverse_fastq == NULL) {
                perror("Failed to allocate memory for reverse fastq files");
                exit(EXIT_FAILURE);
            }
        }
        //find concatenated string length
        size_t fwd_len = forwardFastqFilesString ? strlen(forwardFastqFilesString) : 0;
        size_t rev_len = reverseFastqFilesString ? strlen(reverseFastqFilesString) : 0;
        size_t bar_len = barcodeFastqFilesString ? strlen(barcodeFastqFilesString) : 0;
        int concatenated_length = bar_len + fwd_len + rev_len + 3;
        fastq_files->concatenated_files = malloc(concatenated_length);
        if (fastq_files->concatenated_files == NULL) {
            perror("Failed to allocate memory for concatenated fastq files");
            exit(EXIT_FAILURE);
        }
        char *this_file = fastq_files->concatenated_files;
        if (barcodeFastqFilesString) {
            put_fastq_files_string_into_collection(barcodeFastqFilesString, fastq_files->barcode_fastq, &fastq_files->nbarcode_files, this_file);
        } else {
            perror("Must have some barcode fastq files");
            exit(EXIT_FAILURE);
        }
        put_fastq_files_string_into_collection(forwardFastqFilesString, fastq_files->forward_fastq, &fastq_files->nforward_files, this_file);
        put_fastq_files_string_into_collection(reverseFastqFilesString, fastq_files->reverse_fastq, &fastq_files->nreverse_files, this_file);
    }
    check_filecounts(fastq_files);
    fastq_files->nsamples=0;
    int name_length=0;
    for (int i = 0; i < fastq_files->nbarcode_files; i++) {
        name_length+=strlen(fastq_files->barcode_fastq[i])+1;
    }
    fprintf(stderr, "Name length %d\n", name_length);
    fastq_files->concatenated_sample_names=malloc(name_length+1);
    fastq_files->sample_sizes=malloc(fastq_files->nbarcode_files*sizeof(int));
    fastq_files->sample_names=malloc(fastq_files->nbarcode_files*sizeof(char*));
    fastq_files->sample_offsets=malloc(fastq_files->nbarcode_files*sizeof(int));
    //check if the memory allocation was successful
    if (fastq_files->concatenated_sample_names == NULL || fastq_files->sample_sizes == NULL) {
        perror("Failed to allocate memory for sample names");
        exit(EXIT_FAILURE);
    }
    memset(fastq_files->concatenated_sample_names, 0, strlen(fastq_files->concatenated_files)+1);
    char name_copy[name_length+1];
    char *last_sample=fastq_files->concatenated_sample_names;
    char *sample_name=0;
    fastq_files->nsamples=0;
    fastq_files->sample_offsets[0]=0;
    if (sample_flag){
        //assume that files are sorted by name
        for (int i = 0; i < fastq_files->nbarcode_files; i++) {
            strcpy(name_copy, fastq_files->barcode_fastq[i]); //copy is nccecessary because extract_sample_name modifies the string
            sample_name=extract_sample_name(name_copy, barcode_pattern);
            fprintf(stderr, "Sample name %s\n", sample_name);
            fprintf(stderr, "last sample name %s\n", last_sample);
            if (!fastq_files->nsamples){
                strcpy(last_sample, sample_name);
                fastq_files->sample_sizes[0]=1;
                fastq_files->sample_names[0]=last_sample;
                fastq_files->nsamples=1;
            }
            else if (strcmp(last_sample, sample_name)){
                last_sample+=strlen(last_sample)+1;
                fastq_files->sample_offsets[fastq_files->nsamples]=i;
                fastq_files->sample_names[fastq_files->nsamples]=last_sample;
                strcpy(last_sample, sample_name);
                fastq_files->sample_sizes[fastq_files->nsamples]=1;
                fastq_files->nsamples++;
            }
            else{
                fastq_files->sample_sizes[fastq_files->nsamples-1]++;
            }
        }
    }
    else{
        fastq_files->nsamples=1;
        fastq_files->sample_sizes[0]=fastq_files->nbarcode_files;
    }
    //now allocate the sample buffers
    fastq_files->sorted_index=malloc(fastq_files->nsamples*sizeof(int));
    if (!fastq_files->sorted_index){
        perror("Failed to allocate memory for sorted index");
        exit(EXIT_FAILURE);
    }
    sort_samples_by_size(fastq_files, fastq_files->sorted_index);

    int max_sample_size=0;
    for (int i=0; i<fastq_files->nsamples; i++){
        if (fastq_files->sample_sizes[i]>max_sample_size){
            max_sample_size=fastq_files->sample_sizes[i];
        }
    }
    fastq_files->max_sample_size=max_sample_size;
    //check that memory allocations are non null and free them
}

void read_barcodes_into_hash(const char *filename, khash_t(strptr)* hash) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Failed to open barcode file %s\n", filename);
        exit(EXIT_FAILURE);
    }
    char line[LINE_LENGTH];
    while (fgets(line, LINE_LENGTH, file) != NULL) {
        line[strcspn(line, "\r\n")] = 0;
        char *line_copy = strdup(line);
        int ret;
        khint_t k = kh_put(strptr, hash, line_copy, &ret);
        kh_val(hash, k) = (void*)1; // Store as pointer to indicate presence
    }
    fclose(file);
}

void free_strptr_hash(khash_t(strptr)* hash) {
    if (!hash) return;
    khint_t k;
    for (k = kh_begin(hash); k != kh_end(hash); ++k) {
        if (kh_exist(hash, k)) {
            free((char*)kh_key(hash, k));  // Free the strdup'd key
        }
    }
    kh_destroy(strptr, hash);
}

int pf_file_exists(const char *filename){
    struct stat buffer;
    return (stat(filename, &buffer) == 0);
}
const char* get_basename(const char* path) {
    size_t len = strlen(path);
    while (len > 0 && path[len-1] == '/') len--;       // strip trailing slashes
    if (len == 0) return path;
    const char *start = path;
    for (size_t i = 0; i < len; ++i) if (path[i] == '/') start = path + i + 1;
    return start;
}

/* ------------------------------------------------------------------
 * Sample-size helper functions shared by assignBarcodes and demux_fastq
 * ------------------------------------------------------------------ */

size_t get_file_size(char *filepath)
{
    struct stat st;
    if (stat(filepath, &st) == 0)
        return (size_t)st.st_size;
    perror(filepath);
    return 0;
}

/* tiny selection-sort with context to avoid pulling qsort_r (BSD/GNU diff) */
static int compare_file_sizes_ctx(const void *a, const void *b, void *ctx)
{
    const size_t *sizes = (const size_t*)ctx;
    int ia = *(const int*)a;
    int ib = *(const int*)b;
    if (sizes[ia] < sizes[ib]) return 1;   /* we want descending order */
    if (sizes[ia] > sizes[ib]) return -1;
    return 0;
}

static void qsort_with_ctx(void *base, size_t n, size_t size,
                           int (*cmp)(const void*, const void*, void*),
                           void *ctx)
{
    char *p = (char*)base;
    for (size_t i = 0; i + 1 < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (cmp(p + i*size, p + j*size, ctx) > 0) {
                char tmp[size];
                memcpy(tmp,        p + i*size, size);
                memcpy(p + i*size, p + j*size, size);
                memcpy(p + j*size, tmp,        size);
            }
        }
    }
}

void sort_samples_by_size(fastq_files_collection *fastq_files, int *sample_order)
{
    if (!fastq_files || fastq_files->nsamples <= 0) return;

    size_t *sizes = (size_t*)malloc((size_t)fastq_files->nsamples * sizeof(size_t));
    if (!sizes) { perror("malloc sizes"); exit(EXIT_FAILURE); }

    int idx = 0;
    for (int s = 0; s < fastq_files->nsamples; ++s) {
        size_t total = 0;
        for (int j = 0; j < fastq_files->sample_sizes[s]; ++j) {
            total += get_file_size(fastq_files->barcode_fastq[idx]);
            if (fastq_files->forward_fastq)
                total += get_file_size(fastq_files->forward_fastq[idx]);
            if (fastq_files->reverse_fastq)
                total += get_file_size(fastq_files->reverse_fastq[idx]);
            idx++;
        }
        sizes[s] = total;
        sample_order[s] = s;
    }

    qsort_with_ctx(sample_order, (size_t)fastq_files->nsamples, sizeof(int),
                   compare_file_sizes_ctx, sizes);
    free(sizes);
}

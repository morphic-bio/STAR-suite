#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../include/io.h"

static int failures = 0;

#define CHECK(condition, message) do { \
    if (!(condition)) { \
        fprintf(stderr, "FAIL: %s\n", (message)); \
        failures++; \
    } else { \
        fprintf(stderr, "PASS: %s\n", (message)); \
    } \
} while (0)

static void touch_file(const char *directory, const char *name) {
    size_t length = strlen(directory) + strlen(name) + 2;
    char *path = malloc(length);
    if (!path) exit(EXIT_FAILURE);
    snprintf(path, length, "%s/%s", directory, name);
    FILE *handle = fopen(path, "w");
    if (!handle) {
        perror(path);
        exit(EXIT_FAILURE);
    }
    fclose(handle);
    free(path);
}

static void free_paths(char **paths, int count) {
    if (!paths) return;
    for (int i = 0; i < count; i++) free(paths[i]);
    free(paths);
}

static void remove_fixture(const char *directory, const char **names, int count) {
    for (int i = 0; i < count; i++) {
        size_t length = strlen(directory) + strlen(names[i]) + 2;
        char *path = malloc(length);
        if (!path) exit(EXIT_FAILURE);
        snprintf(path, length, "%s/%s", directory, names[i]);
        unlink(path);
        free(path);
    }
    rmdir(directory);
}

static void test_provider_r3_prefix(void) {
    char directory[] = "/tmp/pf-fastq-r3-prefix-XXXXXX";
    const char *names[] = {
        "CP_R3_gRNA_A1_S3_L001_R1_001.fastq.gz",
        "CP_R3_gRNA_A1_S3_L001_R2_001.fastq.gz",
        "CP_R3_gRNA_A1_S3_L001_I1_001.fastq.gz"
    };
    char **r1 = NULL, **r2 = NULL, **r3 = NULL;
    int n1 = 0, n2 = 0, n3 = 0;

    CHECK(mkdtemp(directory) != NULL, "create provider-prefix fixture");
    for (int i = 0; i < 3; i++) touch_file(directory, names[i]);
    int status = find_paired_fastq_files(directory, "_R1_", "_R2_", "_R3_",
                                         &r1, &n1, &r2, &n2, &r3, &n3);
    CHECK(status == 0, "R3 provider prefix resolves successfully");
    CHECK(n1 == 1 && n2 == 1 && n3 == 0, "R3 provider prefix is not classified as a reverse read");
    CHECK(n1 == 1 && strcmp(get_basename(r1[0]), names[0]) == 0, "correct R1 selected with R3 prefix");
    CHECK(n2 == 1 && strcmp(get_basename(r2[0]), names[1]) == 0, "correct R2 mate selected with R3 prefix");

    free_paths(r1, n1);
    free_paths(r2, n2);
    free_paths(r3, n3);
    remove_fixture(directory, names, 3);
}

static void test_duplicate_r1_token_resolved_by_mate(void) {
    char directory[] = "/tmp/pf-fastq-r1-prefix-XXXXXX";
    const char *names[] = {
        "CP_R1_gRNA_A1_S3_L001_R1_001.fastq.gz",
        "CP_R1_gRNA_A1_S3_L001_R2_001.fastq.gz"
    };
    char **r1 = NULL, **r2 = NULL, **r3 = NULL;
    int n1 = 0, n2 = 0, n3 = 0;

    CHECK(mkdtemp(directory) != NULL, "create duplicate-R1 fixture");
    for (int i = 0; i < 2; i++) touch_file(directory, names[i]);
    int status = find_paired_fastq_files(directory, "_R1_", "_R2_", "_R3_",
                                         &r1, &n1, &r2, &n2, &r3, &n3);
    CHECK(status == 0, "duplicate R1 token resolves by existing mate");
    CHECK(n1 == 1 && n2 == 1 && n3 == 0, "only the read-position R1 token is used");
    CHECK(n2 == 1 && strcmp(get_basename(r2[0]), names[1]) == 0, "duplicate-token R2 mate is correct");

    free_paths(r1, n1);
    free_paths(r2, n2);
    free_paths(r3, n3);
    remove_fixture(directory, names, 2);
}

static void test_true_three_read_layout(void) {
    char directory[] = "/tmp/pf-fastq-three-read-XXXXXX";
    const char *names[] = {
        "sample_L001_R1_001.fastq.gz",
        "sample_L001_R2_001.fastq.gz",
        "sample_L001_R3_001.fastq.gz"
    };
    char **r1 = NULL, **r2 = NULL, **r3 = NULL;
    int n1 = 0, n2 = 0, n3 = 0;

    CHECK(mkdtemp(directory) != NULL, "create three-read fixture");
    for (int i = 0; i < 3; i++) touch_file(directory, names[i]);
    int status = find_paired_fastq_files(directory, "_R1_", "_R2_", "_R3_",
                                         &r1, &n1, &r2, &n2, &r3, &n3);
    CHECK(status == 0, "true R1/R2/R3 layout resolves successfully");
    CHECK(n1 == 1 && n2 == 1 && n3 == 1, "true R3 mate is retained");

    free_paths(r1, n1);
    free_paths(r2, n2);
    free_paths(r3, n3);
    remove_fixture(directory, names, 3);
}

static void test_ambiguous_positions_fail(void) {
    char directory[] = "/tmp/pf-fastq-ambiguous-XXXXXX";
    const char *names[] = {
        "CP_R1_gRNA_L001_R1_001.fastq.gz",
        "CP_R2_gRNA_L001_R1_001.fastq.gz",
        "CP_R1_gRNA_L001_R2_001.fastq.gz"
    };
    char **r1 = NULL, **r2 = NULL, **r3 = NULL;
    int n1 = 0, n2 = 0, n3 = 0;

    CHECK(mkdtemp(directory) != NULL, "create ambiguous-position fixture");
    for (int i = 0; i < 3; i++) touch_file(directory, names[i]);
    int status = find_paired_fastq_files(directory, "_R1_", "_R2_", "_R3_",
                                         &r1, &n1, &r2, &n2, &r3, &n3);
    CHECK(status == -1, "multiple mate-producing R1 positions fail explicitly");

    free_paths(r1, n1);
    free_paths(r2, n2);
    free_paths(r3, n3);
    remove_fixture(directory, names, 3);
}

int main(void) {
    test_provider_r3_prefix();
    test_duplicate_r1_token_resolved_by_mate();
    test_true_three_read_layout();
    test_ambiguous_positions_fail();
    if (failures) {
        fprintf(stderr, "%d FASTQ pairing test(s) failed\n", failures);
        return EXIT_FAILURE;
    }
    fprintf(stderr, "All FASTQ pairing tests passed\n");
    return EXIT_SUCCESS;
}

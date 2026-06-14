#include "../include/barcode_match.h"
#include "../include/globals.h"
#include "../include/io.h"
#include "../include/prototypes.h"

#include <fcntl.h>
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

static int write_duplicate_feature_csv(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) {
        return -1;
    }
    fprintf(f, "name,sequence,pattern\n");
    fprintf(f, "GuideA,ACGTACGT,NNNN(BC)\n");
    fprintf(f, "GuideB,TGCATGCA,NNNN(BC)\n");
    fprintf(f, "GuideA,CCCCCCCC,AAAA(BC)\n");
    fprintf(f, "GuideC,GATCGATC,NNNN(BC)\n");
    fprintf(f, "GuideB,GGGGGGGG,TTTT(BC)\n");
    fclose(f);
    return 0;
}

static char *read_file(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) {
        return NULL;
    }
    if (fseek(f, 0, SEEK_END) != 0) {
        fclose(f);
        return NULL;
    }
    long len = ftell(f);
    if (len < 0 || fseek(f, 0, SEEK_SET) != 0) {
        fclose(f);
        return NULL;
    }
    char *buf = calloc((size_t)len + 1, 1);
    if (!buf) {
        fclose(f);
        return NULL;
    }
    size_t nread = fread(buf, 1, (size_t)len, f);
    buf[nread] = '\0';
    fclose(f);
    return buf;
}

int main(void) {
    barcode_match_init();
    initialize_complement();
    initialize_unit_sizes();

    char tmp_template[] = "/tmp/pf_feature_dups_XXXXXX";
    char *tmp_dir = mkdtemp(tmp_template);
    TEST_ASSERT(tmp_dir != NULL, "mkdtemp");

    char csv_path[512];
    char stderr_path[512];
    snprintf(csv_path, sizeof(csv_path), "%s/features.csv", tmp_dir);
    snprintf(stderr_path, sizeof(stderr_path), "%s/stderr.txt", tmp_dir);
    TEST_ASSERT(write_duplicate_feature_csv(csv_path) == 0, "write feature CSV");

    int saved_stderr = dup(STDERR_FILENO);
    TEST_ASSERT(saved_stderr >= 0, "dup stderr");
    int stderr_fd = open(stderr_path, O_CREAT | O_TRUNC | O_WRONLY, 0600);
    TEST_ASSERT(stderr_fd >= 0, "open stderr capture");
    TEST_ASSERT(dup2(stderr_fd, STDERR_FILENO) >= 0, "redirect stderr");

    feature_arrays *features = read_features_file(csv_path);

    fflush(stderr);
    TEST_ASSERT(dup2(saved_stderr, STDERR_FILENO) >= 0, "restore stderr");
    close(saved_stderr);
    close(stderr_fd);

    TEST_ASSERT(features != NULL, "read features");
    TEST_ASSERT(features->number_of_features == 3, "duplicate feature names ignored");
    TEST_ASSERT(strcmp(features->feature_names[0], "GuideA") == 0, "GuideA kept first");
    TEST_ASSERT(strcmp(features->feature_sequences[0], "ACGTACGT") == 0, "GuideA first sequence kept");
    TEST_ASSERT(strcmp(features->feature_names[1], "GuideB") == 0, "GuideB kept first");
    TEST_ASSERT(strcmp(features->feature_sequences[1], "TGCATGCA") == 0, "GuideB first sequence kept");
    TEST_ASSERT(strcmp(features->feature_names[2], "GuideC") == 0, "GuideC retained");

    char *stderr_text = read_file(stderr_path);
    TEST_ASSERT(stderr_text != NULL, "read stderr capture");
    TEST_ASSERT(strstr(stderr_text, "duplicate feature name 'GuideA'") != NULL, "GuideA warning emitted");
    TEST_ASSERT(strstr(stderr_text, "duplicate feature name 'GuideB'") != NULL, "GuideB warning emitted");
    TEST_ASSERT(strstr(stderr_text, "ignoring later definition") != NULL, "warning says later definition ignored");

    free(stderr_text);
    free_feature_arrays(features);
    feature_code_hash.h64 = NULL;

    unlink(csv_path);
    unlink(stderr_path);
    rmdir(tmp_dir);
    printf("Feature duplicate-name test passed.\n");
    return 0;
}

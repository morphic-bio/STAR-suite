/* Interleaved contexts retain cloned settings and their own feature references. */
#include "../include/pf_api.h"
#include "../include/globals.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void write_text(const char *path, const char *text) {
    FILE *file = fopen(path, "w");
    assert(file && fputs(text, file) >= 0);
    assert(fclose(file) == 0);
}

int main(void) {
    char directory[] = "/tmp/pf-context-isolation-XXXXXX";
    assert(mkdtemp(directory));
    char fa[256], fb[256], wl[256];
    snprintf(fa, sizeof(fa), "%s/a.csv", directory);
    snprintf(fb, sizeof(fb), "%s/b.csv", directory);
    snprintf(wl, sizeof(wl), "%s/whitelist.txt", directory);
    write_text(fa, "name,sequence,pattern\nA1,ACGTACGT,NNNN(BC)\nA2,TGCATGCA,NNNN(BC)\n");
    write_text(fb, "name,sequence,pattern\nB1,CCCCGGGG,NNNN(BC)\n");
    write_text(wl, "AAACCCAAGTAACCCC\n");

    pf_config *config = pf_config_create();
    assert(config);
    pf_config_set_umi_length(config, 10);
    pf_config_set_translate_nxt(config, 0);
    pf_context *a = pf_init(config);
    assert(a);
    pf_config_set_umi_length(config, 12);
    pf_config_set_translate_nxt(config, 1);
    pf_context *b = pf_init(config);
    assert(b);
    pf_config_destroy(config);

    assert(pf_load_feature_ref(a, fa) == PF_OK);
    assert(umi_length == 10 && translate_NXT == 0);
    assert(pf_load_feature_ref(b, fb) == PF_OK);
    assert(umi_length == 12 && translate_NXT == 1);
    assert(pf_get_num_features(a) == 2 && pf_get_num_features(b) == 1);
    assert(strcmp(pf_get_feature_name(a, 0), "A1") == 0);
    assert(strcmp(pf_get_feature_name(b, 0), "B1") == 0);

    for (int repeat = 0; repeat < 3; ++repeat) {
        assert(pf_load_whitelist(a, wl) == PF_OK);
        assert(umi_length == 10 && translate_NXT == 0 && number_of_features == 2);
        assert(pf_load_whitelist(b, wl) == PF_OK);
        assert(umi_length == 12 && translate_NXT == 1 && number_of_features == 1);
    }
    /* Destroy the inactive context, then reuse the surviving one. */
    pf_destroy(a);
    assert(pf_load_whitelist(b, wl) == PF_OK);
    assert(strcmp(pf_get_feature_name(b, 0), "B1") == 0);
    assert(strcmp(pf_get_feature_sequence(b, 0), "CCCCGGGG") == 0);
    assert(umi_length == 12 && translate_NXT == 1 && number_of_features == 1);
    pf_destroy(b);
    assert(unlink(fa) == 0 && unlink(fb) == 0 && unlink(wl) == 0);
    assert(rmdir(directory) == 0);
    puts("Context isolation checks passed.");
    return 0;
}

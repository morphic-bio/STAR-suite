/**
 * Regression test for prehash ambiguity bookkeeping.
 *
 * Reproduces: a third same-distance collider arriving after a key is
 * already marked ambiguous must still clear no_ambiguity[] for the
 * third feature.  Before the fix, the early-return in
 * feature_hash_insert_best skipped the no_ambiguity update.
 *
 * Setup (4-base features for minimal combinatorics):
 *   A = "CAAA"  (index 0, 1-based 1)
 *   B = "ACAA"  (index 1, 1-based 2)
 *   C = "AACA"  (index 2, 1-based 3)
 *
 * All three generate the 1-mismatch variant "AAAA":
 *   A at pos 0: C->A  =>  AAAA  (d=1)
 *   B at pos 1: C->A  =>  AAAA  (d=1)
 *   C at pos 2: C->A  =>  AAAA  (d=1)
 *
 * Expected: no_ambiguity_le1[0..2] == 0 for all three features.
 * Bug: no_ambiguity_le1[2] was stale (remained 1) because the third
 *      insert hit the "already ambiguous" early-return without updating.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../include/pf_api.h"

static void write_file(const char *path, const char *content) {
    FILE *fp = fopen(path, "w");
    if (!fp) { perror(path); exit(1); }
    fputs(content, fp);
    fclose(fp);
}

int main(void) {
    char tmp_template[] = "/tmp/pf_prehash_ambig_XXXXXX";
    char *tmp_dir = mkdtemp(tmp_template);
    if (!tmp_dir) { perror("mkdtemp"); return 1; }

    char feature_csv[512], whitelist[512];
    snprintf(feature_csv, sizeof(feature_csv), "%s/features.csv", tmp_dir);
    snprintf(whitelist, sizeof(whitelist), "%s/whitelist.txt", tmp_dir);

    write_file(feature_csv,
        "id,name,sequence,feature_type\n"
        "f1,FeatA,CAAA,CRISPR Guide Capture\n"
        "f2,FeatB,ACAA,CRISPR Guide Capture\n"
        "f3,FeatC,AACA,CRISPR Guide Capture\n");

    write_file(whitelist, "AAAA\n");

    pf_config *cfg = pf_config_create();
    pf_config_set_max_hamming_distance(cfg, 1);
    pf_config_set_barcode_length(cfg, 4);
    pf_config_set_umi_length(cfg, 0);

    pf_context *ctx = pf_init(cfg);
    pf_config_destroy(cfg);
    if (!ctx) {
        fprintf(stderr, "FAIL: pf_init returned NULL\n");
        return 1;
    }

    pf_error err = pf_load_feature_ref(ctx, feature_csv);
    if (err != PF_OK) {
        fprintf(stderr, "FAIL: pf_load_feature_ref: %s\n", pf_get_error(ctx));
        pf_destroy(ctx);
        return 1;
    }

    int nf = pf_get_num_features(ctx);
    if (nf != 3) {
        fprintf(stderr, "FAIL: expected 3 features, got %d\n", nf);
        pf_destroy(ctx);
        return 1;
    }

    int pass = 1;

    for (int i = 0; i < 3; i++) {
        int val = pf_get_feature_no_ambiguity(ctx, 1, i);
        if (val == -1) {
            fprintf(stderr, "FAIL: le1 prehash not built (feature %d returned -1)\n", i);
            pass = 0;
            continue;
        }
        if (val != 0) {
            fprintf(stderr,
                "FAIL: no_ambiguity_le1[%d] = %d, expected 0 "
                "(feature %s participates in ambiguous key \"AAAA\" at d=1)\n",
                i, val, pf_get_feature_name(ctx, i));
            pass = 0;
        } else {
            fprintf(stderr,
                "  OK: no_ambiguity_le1[%d] = 0 (feature %s correctly marked)\n",
                i, pf_get_feature_name(ctx, i));
        }
    }

    pf_destroy(ctx);

    unlink(feature_csv);
    unlink(whitelist);
    rmdir(tmp_dir);

    if (pass) {
        printf("PASS: prehash ambiguity bookkeeping for third collider\n");
        return 0;
    }
    printf("FAIL: prehash ambiguity bookkeeping regression\n");
    return 1;
}

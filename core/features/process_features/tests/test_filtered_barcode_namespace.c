/* Filtered barcode lookup is exact after explicit ingress normalization. */
#include "../include/common.h"
#include "../include/io.h"
#include "../include/pf_api.h"
#include "../include/prototypes.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main(void) {
    const char *tru = "AAACCCAAGTAACCCC";
    const char *nxt = "AAACCCATCTAACCCC";
    khash_t(strptr) *hash = kh_init(strptr);
    assert(hash);
    int added;
    khint_t key = kh_put(strptr, hash, strdup(tru), &added);
    assert(added > 0);
    kh_val(hash, key) = NULL;
    assert(filtered_barcode_hash_contains(hash, tru));
    assert(!filtered_barcode_hash_contains(hash, nxt));
    assert(pf_normalize_hash_namespace(hash) == 1);
    assert(!filtered_barcode_hash_contains(hash, tru));
    assert(filtered_barcode_hash_contains(hash, nxt));
    assert(pf_normalize_hash_namespace(hash) == 1);
    assert(filtered_barcode_hash_contains(hash, tru));
    assert(!filtered_barcode_hash_contains(hash, nxt));
    free_strptr_hash(hash);

    char path[] = "/tmp/pf-filtered-namespace-XXXXXX";
    int fd = mkstemp(path);
    assert(fd >= 0);
    FILE *file = fdopen(fd, "w");
    assert(file);
    assert(fprintf(file, "%s\n", tru) > 0);
    assert(fclose(file) == 0);
    const pf_namespace_t namespaces[] = {PF_NS_UNKNOWN, PF_NS_NXT, PF_NS_TRU, PF_NS_UNION};
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            pf_config *config = pf_config_create();
            assert(config);
            pf_config_set_source_namespace(config, namespaces[i]);
            pf_config_set_target_namespace(config, namespaces[j]);
            pf_context *ctx = pf_init(config);
            pf_config_destroy(config);
            assert(ctx);
            int both_known = (i == 1 || i == 2) && (j == 1 || j == 2);
            assert(pf_load_filtered_barcodes(ctx, path) == (both_known ? PF_OK : PF_ERR_NAMESPACE));
            pf_destroy(ctx);
        }
    }
    pf_config *config = pf_config_create();
    assert(config);
    pf_config_set_allow_union_whitelist(config, 1);
    pf_context *ctx = pf_init(config);
    pf_config_destroy(config);
    assert(ctx && pf_load_filtered_barcodes(ctx, path) == PF_OK);
    pf_destroy(ctx);
    assert(unlink(path) == 0);
    puts("Filtered barcode namespace checks passed.");
    return 0;
}

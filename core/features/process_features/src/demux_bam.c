#include "../include/common.h"
#include "../include/utils.h"
#include "../include/globals.h"
#include "../include/io.h"
#include "../include/barcode_match.h"
#include "../include/prototypes.h"   /* for code2string() */
#include "../include/memory.h"
#include <htslib/sam.h>
#include <endian.h>
#include <libgen.h>

static memory_pool_collection *global_pools = NULL;

// Memory duplication helper (replaces g_memdup/g_memdup2)
#define MEMDUP(ptr, size) ({ \
    void *_dst = malloc(size); \
    if (_dst) memcpy(_dst, ptr, size); \
    _dst; \
})

#define PROBE_LEN 8

/*
 * demux_bam.c – Read STARsolo BAM, deduplicate (CB,UB,gene) and output MatrixMarket.
 * Features = sample barcode indices (probe variants). Dedup ignores sample index.
 */

#define DEFAULT_CB_TAG "CB"
#define DEFAULT_UB_TAG "UB"
#define DEFAULT_GENE_TAG "GX"
#define FALLBACK_GENE_TAG "GE"

/* ---------------- Intern tables ---------------- */

typedef struct {
    khash_t(strptr) *str_to_id;   /* char*  -> uint32_t */
    vec_ptr_t  *id_to_str;   /* index -> char*   */
} intern_table;

static inline uint32_t intern_get(intern_table *t, const char *s) {
    khint_t k = kh_get(strptr, t->str_to_id, s);
    if (k != kh_end(t->str_to_id)) return (uint32_t)(uintptr_t)kh_val(t->str_to_id, k);
    uint32_t new_id = vec_ptr_size(t->id_to_str);
    char *s_copy = strdup(s);
    int ret;
    khint_t kh = kh_put(strptr, t->str_to_id, s_copy, &ret);
    kh_val(t->str_to_id, kh) = (void*)(uintptr_t)new_id;
    vec_ptr_add(t->id_to_str, s_copy);
    return new_id;
}

static inline const char* id_to_str(intern_table *t, uint32_t id) {
    return (const char*)vec_ptr_get(t->id_to_str, id);
}

static intern_table* intern_table_new(void) {
    intern_table *t = calloc(1, sizeof(intern_table));
    t->str_to_id = kh_init(strptr);
    t->id_to_str = vec_ptr_init(free);
    return t;
}

/* ---------------- Sample probe features ---------------- */
static feature_arrays *probe_fa = NULL;
static int direct_probe = -1; /* -1 auto; 0 hash; 1 direct */
static int probe_offset = 68;  /* default */

static inline char bam_code_to_base(uint8_t c){
    switch(c){
        case 1: return 'A';
        case 2: return 'C';
        case 4: return 'G';
        case 8: return 'T';
        default: return 'N';
    }
}

static inline char complement_base(char b)
{
    switch (b) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default : return 'N';
    }
}

static int extract_kmer_from_bam(const bam1_t *b, int offset,
                                      char out8[PROBE_LEN + 1])
{
    int n = b->core.l_qseq;
    if (offset < 0 || offset + PROBE_LEN > n) return 0;

    const uint8_t *seq = bam_get_seq(b);

    if (b->core.flag & BAM_FREVERSE) {              /* reverse-strand read */
        int start = n - offset - PROBE_LEN;         /* 5′ on reference   */
        if (start < 0) return 0;
        for (int i = 0; i < PROBE_LEN; ++i) {
            /* walk right→left through the read, left→right into out8 */
            char base = bam_code_to_base(bam_seqi(seq,
                                                  start + PROBE_LEN - 1 - i));
            if (base != 'A' && base != 'C' && base != 'G' && base != 'T')
                return 0;
            out8[i] = complement_base(base);
        }
    } else {                                        /* forward-strand read */
        for (int i = 0; i < PROBE_LEN; ++i) {
            char base = bam_code_to_base(bam_seqi(seq, offset + i));
            if (base != 'A' && base != 'C' && base != 'G' && base != 'T')
                return 0;
            out8[i] = base;
        }
    }

    out8[PROBE_LEN] = '\0';
    return 1;
}


/* TripKey dedup helpers removed */

/* ---------------- Shard data ---------------- */

typedef struct {
    khash_t(u32ptr) *cb_counts;     /* key: cb_id (uint32_t)  val: uint16_t[nprobes] */
    uint64_t     total_reads;
    uint64_t     usable_reads;
    uint64_t    *probe_tot;     /* per-probe total usable reads */
} shard_data;

static shard_data* shard_data_new(int nprobes){
    shard_data *s = calloc(1, sizeof(shard_data));
    s->cb_counts = kh_init(u32ptr);
    s->probe_tot = calloc(nprobes, sizeof(uint64_t));
    return s;
}


/* ---- helper: qsort comparator for barcode strings ----------------- */
static int cmp_str (const void *a, const void *b)
{
    return strcmp (*(char* const*)a, *(char* const*)b);
}
/* ------------------------------------------------------------------- */

/* ---------------- CLI config ---------------- */

typedef struct {
    const char *bam_path;
    const char *outdir;
    const char *sample_map_path;
    const char *cb_tag;
    const char *ub_tag;
    const char *gene_tag;
    int hts_threads;
    int consumer_threads; /* shards */
    int min_mapq;
    int primary_only;
    int skip_dup;
    long long max_records;
    int debug;
    const char *sample_probes; /* path */
    int probe_offset;
    int direct_probe;
    int save_read_to_cb;
    int search_nearby; /* if non-zero, enable fallback search for probes */
    int count_intergene;        /* keep GX starting with '-'     */
    const char *cbub_path;
    int use_cbub;
    const char *whitelist_path;
    const char *bin_out_path;   /* binary sidecar output path */
    int use_bin_out;             /* boolean: binary output mode active */
} cfg_t;

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s --bam FILE --outdir DIR [options]\n"
        "Options:\n"
        "  --cell_sample_map FILE   Optional CB→sample mapping TSV (unused in probe mode)\n"
        "  --sample_probes FILE     Probe 8-mer table (variant,unused,BCname)\n"
        "  --probe_offset N         0-based probe offset (default 68)\n"
        "  --cb_tag STR             Default CB\n"
        "  --ub_tag STR             Default UB\n"
        "  --gene_tag STR           Default GX (fallback GE)\n"
        "  --hts_threads N          BGZF threads (default 2) or -S N\n"
        "  --threads N              Consumer shards (default 1) or -t N\n"
        "  --min_mapq N             Minimum MAPQ (default 0)\n"
        "  --no_primary_filter      Keep secondary/supplementary\n"
        "  --keep_dup               Keep duplicate reads\n"
        "  --max_records N          Process first N records only\n"
        "  -v                       Debug\n"
        "  --search_nearby           Try offsets ±1,±2 for probe search\n"
        "  --count_intergene         Include reads whose GX tag starts with '-'\n"
        "  --CBUB_file FILE         Optional STAR cb/ub tag stream (Aligned.out.cb_ub.bin)\n"
        "  --whitelist FILE         Cell barcode whitelist (required for CBUB mode)\n"
        "  --bin_out FILE           Binary sidecar output (suppresses MatrixMarket TSVs)\n",
        prog);
}

static feature_arrays* load_probe_variants_to_features(const char *path) {
    /* Three-column TSV: VARIANT_8MER  CANONICAL_8MER  BARCODE_ID */
    FILE *fp = fopen(path, "r");
    if (!fp) { perror("sample_probes"); exit(EXIT_FAILURE); }

    /* first pass – collect unique barcode IDs and their canonical sequence */
    khash_t(strptr) *bc2canon = kh_init(strptr);
    char line[LINE_LENGTH];
    int uniq_cnt = 0, name_size = 0, seq_size = 0, code_size = 0, variant_cnt = 0;
    while (fgets(line, sizeof(line), fp)) {
        char *save=NULL;
        char *variant = strtok_r(line, "\t \r\n", &save);
        char *canonical = strtok_r(NULL, "\t \r\n", &save);
        char *bc_id   = strtok_r(NULL, "\t \r\n", &save);
        if (!variant || !bc_id) continue;
        if ((int)strlen(variant) != PROBE_LEN) continue;
        if (!check_sequence(variant, PROBE_LEN)) continue;
        if (!canonical || strlen(canonical)!=PROBE_LEN || !check_sequence(canonical, PROBE_LEN))
            canonical = variant; /* fallback */
        khint_t k = kh_get(strptr, bc2canon, bc_id);
        if (k == kh_end(bc2canon)) {
            char *bc_id_copy = strdup(bc_id);
            char *canon_copy = strdup(canonical);
            int ret;
            khint_t kh = kh_put(strptr, bc2canon, bc_id_copy, &ret);
            kh_val(bc2canon, kh) = canon_copy;
            uniq_cnt++;
            name_size += (int)strlen(bc_id) + 1;
            seq_size  += PROBE_LEN + 1;
            code_size += (PROBE_LEN+3)/4;
        }
        if (canonical != variant) variant_cnt++;    
    }

    feature_arrays *fa = allocate_feature_arrays(name_size, seq_size, code_size, uniq_cnt, PROBE_LEN);
    fa->common_length = PROBE_LEN;
    //allocate a array of blocks of size (PROBE_LEN+3)/4
    uint8_t *variant_codes = malloc(variant_cnt * (PROBE_LEN+3)/4);
    memset(variant_codes, 0, variant_cnt * (PROBE_LEN+3)/4);

    /* ---------- collect & sort unique BC-IDs ---------- */
    char **bc_list = malloc(uniq_cnt * sizeof(char*));
    {
        khint_t k;
        int i = 0;
        for (k = kh_begin(bc2canon); k != kh_end(bc2canon); ++k) {
            if (kh_exist(bc2canon, k))
                bc_list[i++] = (char*)kh_key(bc2canon, k);          /* key already dup'ed */
        }
    }

    /* alphabetic sort */
    qsort (bc_list, uniq_cnt, sizeof(char*), cmp_str);

    /* second pass – build arrays */
    khash_t(strptr) *bc2idx = kh_init(strptr);
    int idx = 0;
    fa->feature_names[0]     = fa->feature_names_storage;
    fa->feature_sequences[0] = fa->feature_sequences_storage;
    fa->feature_codes[0]     = fa->feature_codes_storage;
    for (idx = 0; idx < uniq_cnt; idx++) {
        const char *bc_id = bc_list[idx];
        khint_t kcanon = kh_get(strptr, bc2canon, bc_id);
        const char *canon = (kcanon != kh_end(bc2canon)) ? (const char*)kh_val(bc2canon, kcanon) : NULL;
        /* store name */
        strcpy(fa->feature_names[idx], bc_id);
        if (idx+1 < uniq_cnt)
            fa->feature_names[idx+1] = fa->feature_names[idx] + strlen(bc_id) + 1;
        /* store canonical sequence */
        strcpy(fa->feature_sequences[idx], canon);
        fa->feature_lengths[idx] = PROBE_LEN;
        /* cast away const because string2code expects char* */
        fa->feature_code_lengths[idx] =
            string2code((char *)canon, PROBE_LEN, fa->feature_codes[idx]);
        if (idx+1 < uniq_cnt) {
            fa->feature_sequences[idx+1] = fa->feature_sequences[idx] + PROBE_LEN + 1;
            fa->feature_codes[idx+1]     = fa->feature_codes[idx] + fa->feature_code_lengths[idx];
        }
        /* insert canonical code into lookup */
        if (fa->feature_lengths[idx] == PROBE_LEN) {
            var_key_t k = {.ptr = fa->feature_codes[idx], .len = fa->feature_code_lengths[idx]};
            int ret;
            khint_t kh = kh_put(codeu32, feature_code_hash, k, &ret);
            kh_val(feature_code_hash, kh) = idx+1;
        }
        int ret;
        khint_t kh = kh_put(strptr, bc2idx, bc_id, &ret);
        kh_val(bc2idx, kh) = (void*)(uintptr_t)(idx + 1);
    }
    free(bc_list);

    /* third pass – map every variant 8-mer to the same index */
    fseek(fp, 0, SEEK_SET);
    int variant_idx = 0;
    while (fgets(line, sizeof(line), fp)) {
        char *save=NULL;
        char *variant = strtok_r(line, "\t \r\n", &save);
        char *canonical = strtok_r(NULL, "\t \r\n", &save);
        char *bc_id   = strtok_r(NULL, "\t \r\n", &save);
        if (!variant || !bc_id) continue;
        if (strcmp(variant, canonical) == 0) continue;
        if ((int)strlen(variant) != PROBE_LEN) continue;
        if (!check_sequence(variant, PROBE_LEN)) continue;

        khint_t kbc = kh_get(strptr, bc2idx, bc_id);
        if (kbc == kh_end(bc2idx)) continue;              /* should not happen */
        uint32_t bc_index = (uint32_t)(uintptr_t)kh_val(bc2idx, kbc) - 1;   /* back to 0-based */
        uint8_t *codebuf = variant_codes + variant_idx * (PROBE_LEN+3)/4;
        int codelen = string2code(variant, PROBE_LEN, codebuf);
        var_key_t k = {.ptr = codebuf, .len = (uint16_t)codelen};
        int ret;
        khint_t kh = kh_put(codeu32, feature_code_hash, k, &ret);
        kh_val(feature_code_hash, kh) = bc_index+1;
        variant_idx++;
    }
    //print the size of the feature_code_hash
    fprintf(stderr, "DEBUG: variant_idx=%d\n", variant_idx);
    fprintf(stderr, "DEBUG: feature_code_hash size=%u\n", kh_size(feature_code_hash));   
    fclose(fp);
    // Free bc2canon keys and values (both are strings)
    khint_t k;
    for (k = kh_begin(bc2canon); k != kh_end(bc2canon); ++k) {
        if (kh_exist(bc2canon, k)) {
            free((char*)kh_key(bc2canon, k));
            free((char*)kh_val(bc2canon, k));
        }
    }
    kh_destroy(strptr, bc2canon);
    kh_destroy(strptr, bc2idx);

    fprintf(stderr, "Loaded %d unique sample barcode IDs (len=%d)\n", uniq_cnt, PROBE_LEN);
    return fa;
}

static void parse_cli(cfg_t *cfg, int argc, char **argv) {
    static struct option long_opts[] = {
        {"bam", required_argument, 0, 0},
        {"outdir", required_argument,0,0},
        {"cell_sample_map", required_argument,0,0},
        {"sample_probes", required_argument,0,0},
        {"probe_offset", required_argument,0,0},
        {"cb_tag", required_argument,0,0},
        {"ub_tag", required_argument,0,0},
        {"gene_tag", required_argument,0,0},
        {"hts_threads", required_argument,0,0},
        {"threads", required_argument,0,0},
        {"min_mapq", required_argument,0,0},
        {"no_primary_filter", no_argument,0,0},
        {"keep_dup", no_argument,0,0},
        {"max_records", required_argument,0,0},
        {"save_read_to_cb", no_argument,0,0},
        {"search_nearby", no_argument, 0, 0},
        {"count_intergene", no_argument, 0, 0},
        {"CBUB_file", required_argument, 0, 0},
        {"whitelist", required_argument, 0, 0},
        {"bin_out", required_argument, 0, 0},
        {"debug", no_argument,0,'v'},
        {0,0,0,0}
    };
    int opt_idx=0;
    while (1) {
        int c = getopt_long(argc, argv, "vS:t:", long_opts, &opt_idx);
        if (c==-1) break;
        if (c=='v') { cfg->debug=1; continue; }
        if (c=='S') { cfg->hts_threads = atoi(optarg); continue; }
        if (c=='t') { cfg->consumer_threads = atoi(optarg); continue; }
        if (c!=0) continue;
        const char *opt_name = long_opts[opt_idx].name;
        if (!strcmp(opt_name, "bam")) cfg->bam_path = optarg;
        else if (!strcmp(opt_name, "outdir")) cfg->outdir = optarg;
        else if (!strcmp(opt_name, "cell_sample_map")) cfg->sample_map_path = optarg;
        else if (!strcmp(opt_name, "sample_probes")) cfg->sample_probes = optarg;
        else if (!strcmp(opt_name, "probe_offset")) cfg->probe_offset = atoi(optarg);
        else if (!strcmp(opt_name, "cb_tag")) cfg->cb_tag = optarg;
        else if (!strcmp(opt_name, "ub_tag")) cfg->ub_tag = optarg;
        else if (!strcmp(opt_name, "gene_tag")) cfg->gene_tag = optarg;
        else if (!strcmp(opt_name, "hts_threads")) cfg->hts_threads = atoi(optarg);
        else if (!strcmp(opt_name, "threads")) cfg->consumer_threads = atoi(optarg);
        else if (!strcmp(opt_name, "min_mapq")) cfg->min_mapq = atoi(optarg);
        else if (!strcmp(opt_name, "no_primary_filter")) cfg->primary_only = 0;
        else if (!strcmp(opt_name, "keep_dup")) cfg->skip_dup = 0;
        else if (!strcmp(opt_name, "max_records")) cfg->max_records = atoll(optarg);
        else if (!strcmp(opt_name, "save_read_to_cb")) cfg->save_read_to_cb = 1;
        else if (!strcmp(opt_name, "search_nearby")) cfg->search_nearby = 1;
        else if (!strcmp(opt_name, "count_intergene")) cfg->count_intergene = 1;
        else if (!strcmp(opt_name, "CBUB_file")) {
            cfg->cbub_path = optarg;
            cfg->use_cbub = 1;
        }
        else if (!strcmp(opt_name, "whitelist")) cfg->whitelist_path = optarg;
        else if (!strcmp(opt_name, "bin_out")) {
            cfg->bin_out_path = optarg;
            cfg->use_bin_out = 1;
        }
        else {
            fprintf(stderr, "Error: unsupported option --%s\n", opt_name);
            usage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }
    if (!cfg->bam_path || !cfg->outdir) { usage(argv[0]); exit(EXIT_FAILURE);}    
    if (cfg->use_cbub && access(cfg->cbub_path, R_OK) != 0) {
        perror("CBUB_file");
        exit(EXIT_FAILURE);
    }
    if (cfg->use_cbub && !cfg->whitelist_path) {
        fprintf(stderr, "Error: --whitelist is required when using --CBUB_file\n");
        exit(EXIT_FAILURE);
    }
    if (cfg->whitelist_path && access(cfg->whitelist_path, R_OK) != 0) {
        perror("whitelist");
        exit(EXIT_FAILURE);
    }
    if (!cfg->cb_tag) cfg->cb_tag = DEFAULT_CB_TAG;
    if (!cfg->ub_tag) cfg->ub_tag = DEFAULT_UB_TAG;
    if (!cfg->gene_tag) cfg->gene_tag = DEFAULT_GENE_TAG;
    if (cfg->hts_threads<=0) cfg->hts_threads=2;
    if (cfg->consumer_threads<=0) cfg->consumer_threads=1;
    if (cfg->probe_offset<=0) cfg->probe_offset = probe_offset;
}

/* ---------------- Whitelist Management ---------------- */

typedef struct {
    char **barcodes;
    size_t count;
} whitelist_t;

static whitelist_t* load_whitelist(const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        perror("load_whitelist");
        exit(EXIT_FAILURE);
    }
    
    whitelist_t *wl = calloc(1, sizeof(whitelist_t));
    char line[256];
    size_t capacity = 1000000; /* Start with 1M capacity */
    wl->barcodes = malloc(capacity * sizeof(char*));
    
    while (fgets(line, sizeof(line), fp)) {
        /* Remove newline */
        char *nl = strchr(line, '\n');
        if (nl) *nl = '\0';
        
        /* Skip empty lines */
        if (strlen(line) == 0) continue;
        
        /* Expand capacity if needed */
        if (wl->count >= capacity) {
            capacity *= 2;
            wl->barcodes = realloc(wl->barcodes, capacity * sizeof(char*));
        }
        
        wl->barcodes[wl->count] = strdup(line);
        wl->count++;
    }
    
    fclose(fp);
    fprintf(stderr, "Loaded %zu barcodes from whitelist\n", wl->count);
    return wl;
}

static void free_whitelist(whitelist_t *wl) {
    if (wl) {
        for (size_t i = 0; i < wl->count; i++) {
            free(wl->barcodes[i]);
        }
        free(wl->barcodes);
        free(wl);
    }
}

/* ---------------- CBUB Reader Abstraction ---------------- */

typedef struct {
    FILE *fp;
    uint64_t status_bits;
    uint64_t cb_bits;
    uint64_t umi_bits;
    uint64_t record_count;
    uint64_t record_idx;
    size_t record_bytes;
    unsigned char *raw;
    size_t raw_cap;
    size_t umi_len;
} cbub_reader;

static cbub_reader* cbub_reader_open(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        perror("cbub_reader_open");
        exit(EXIT_FAILURE);
    }
    
    cbub_reader *r = calloc(1, sizeof(cbub_reader));
    r->fp = fp;
    
    /* Read header: 4 uint64_t values */
    uint64_t header[4];
    if (fread(header, sizeof(uint64_t), 4, fp) != 4) {
        fprintf(stderr, "Error reading CBUB header\n");
        exit(EXIT_FAILURE);
    }
    
    r->status_bits = header[0];
    r->cb_bits = header[1];
    r->umi_bits = header[2];
    r->record_count = header[3];
    
    /* Validate header */
    if (r->status_bits != 1) {
        fprintf(stderr, "Invalid CBUB file: status_bits must be 1, got %lu\n", 
                (unsigned long)r->status_bits);
        exit(EXIT_FAILURE);
    }
    if (r->umi_bits == 0 || r->umi_bits % 2 != 0 || r->umi_bits > 64) {
        fprintf(stderr, "Invalid CBUB file: umi_bits must be even, non-zero, and <= 64, got %lu\n", 
                (unsigned long)r->umi_bits);
        exit(EXIT_FAILURE);
    }
    if (r->cb_bits < 1 || r->cb_bits > 32) {
        fprintf(stderr, "Invalid CBUB file: cb_bits must be >= 1 and <= 32, got %lu\n", 
                (unsigned long)r->cb_bits);
        exit(EXIT_FAILURE);
    }
    if (r->record_count == 0) {
        fprintf(stderr, "Invalid CBUB file: record_count must be > 0, got %lu\n", 
                (unsigned long)r->record_count);
        exit(EXIT_FAILURE);
    }
    
    /* Compute record size and UMI length */
    r->record_bytes = (r->status_bits + r->cb_bits + r->umi_bits + 7) / 8;
    r->umi_len = r->umi_bits / 2;
    
    /* Allocate raw buffer */
    r->raw_cap = r->record_bytes;
    r->raw = malloc(r->raw_cap);
    
    /* Optional: sanity-check file length */
    long current_pos = ftell(fp);
    fseek(fp, 0, SEEK_END);
    long file_size = ftell(fp);
    long expected_size = 32 + r->record_count * r->record_bytes;
    if (file_size != expected_size) {
        fprintf(stderr, "Warning: CBUB file size mismatch. Expected %ld, got %ld\n",
                expected_size, file_size);
    }
    fseek(fp, current_pos, SEEK_SET);
    
    return r;
}

static int cbub_reader_next(cbub_reader *r, int *has_tags, uint32_t *cb_idx, char *umi_out) {
    if (fread(r->raw, 1, r->record_bytes, r->fp) != r->record_bytes) {
        return 0; /* EOF or error */
    }
    
    /* Decode bitstream LSB-first */
    uint64_t accumulator = 0;
    int bits_available = 0;
    int byte_idx = 0;
    
    /* Helper function to get next N bits */
    #define get_bits(n) ({ \
        while (bits_available < (n) && byte_idx < (int)r->record_bytes) { \
            accumulator |= ((uint64_t)r->raw[byte_idx]) << bits_available; \
            bits_available += 8; \
            byte_idx++; \
        } \
        uint64_t result = accumulator & ((1ULL << (n)) - 1); \
        accumulator >>= (n); \
        bits_available -= (n); \
        result; \
    })
    
    /* Extract status (1 bit) */
    uint64_t status = get_bits(1);
    
    /* Extract CB index (cb_bits) */
    uint64_t cb_index = get_bits(r->cb_bits);
    
    /* Extract UMI (umi_bits, decode 2 bits per base) - read in reverse order */
    for (int i = r->umi_len - 1; i >= 0; i--) {
        uint64_t code = get_bits(2);
        switch (code) {
            case 0: umi_out[i] = 'A'; break;
            case 1: umi_out[i] = 'C'; break;
            case 2: umi_out[i] = 'G'; break;
            case 3: umi_out[i] = 'T'; break;
        }
    }
    umi_out[r->umi_len] = '\0';
    
    /* Set outputs */
    *has_tags = (status != 0 && cb_index != 0);
    *cb_idx = (uint32_t)cb_index;
    
    r->record_idx++;
    return 1;
    
    #undef get_bits
}

static void cbub_reader_close(cbub_reader *r) {
    if (r) {
        if (r->fp) fclose(r->fp);
        if (r->raw) free(r->raw);
        free(r);
    }
}

/* Helper function to escape JSON string (simple: escape quotes and backslashes) */
static void json_escape_qname(const char *qname, char *out, size_t out_size) {
    if (!qname) {
        out[0] = '\0';
        return;
    }
    size_t i = 0;
    size_t j = 0;
    while (qname[i] != '\0' && j < out_size - 1) {
        if (qname[i] == '"' || qname[i] == '\\') {
            if (j < out_size - 2) {
                out[j++] = '\\';
                out[j++] = qname[i++];
            } else {
                break;
            }
        } else {
            out[j++] = qname[i++];
        }
    }
    out[j] = '\0';
}

/* ---------------- Worker function (single-thread) ---------------- */

static void process_bam_single(cfg_t *cfg,
                               intern_table *cb_tab,
                               intern_table *gene_tab)
{
    samFile *in = sam_open(cfg->bam_path, "r");
    if (!in) { perror("sam_open"); exit(EXIT_FAILURE);}    
    bam_hdr_t *hdr = sam_hdr_read(in);
    hts_set_threads(in, cfg->hts_threads);

    cbub_reader *cbub = NULL;
    whitelist_t *whitelist = NULL;
    if (cfg->use_cbub) {
        cbub = cbub_reader_open(cfg->cbub_path);
        whitelist = load_whitelist(cfg->whitelist_path);
    }

    /* Debug log setup for NH instrumentation */
    FILE *nh_staging_fp = NULL;
    FILE *nh_tap_fp = NULL;
    long long nh_staging_count = 0;
    const char *nh_staging_path = getenv("STAR_DEBUG_NH_STAGING");
    const char *nh_tap_path = getenv("STAR_DEBUG_NH_TAP");
    const long long NH_STAGING_MAX = 1000;
    
    if (nh_staging_path && nh_staging_path[0] != '\0') {
        char *path_copy = strdup(nh_staging_path);
        char *parent_dir = dirname(path_copy);
        if (strcmp(parent_dir, ".") != 0 && strcmp(parent_dir, "/") != 0) {
            if (mkdir_p(parent_dir)) {
                fprintf(stderr, "Failed to create parent directory for nh_staging: %s\n", parent_dir);
                free(path_copy);
                exit(EXIT_FAILURE);
            }
        }
        free(path_copy);
        
        nh_staging_fp = fopen(nh_staging_path, "w");
        if (!nh_staging_fp) {
            perror("STAR_DEBUG_NH_STAGING");
            exit(EXIT_FAILURE);
        }
    }
    
    if (nh_tap_path && nh_tap_path[0] != '\0') {
        char *path_copy = strdup(nh_tap_path);
        char *parent_dir = dirname(path_copy);
        if (strcmp(parent_dir, ".") != 0 && strcmp(parent_dir, "/") != 0) {
            if (mkdir_p(parent_dir)) {
                fprintf(stderr, "Failed to create parent directory for nh_tap: %s\n", parent_dir);
                free(path_copy);
                exit(EXIT_FAILURE);
            }
        }
        free(path_copy);
        
        nh_tap_fp = fopen(nh_tap_path, "w");
        if (!nh_tap_fp) {
            perror("STAR_DEBUG_NH_TAP");
            exit(EXIT_FAILURE);
        }
    }

    /* Binary output setup */
    FILE *bin_fp = NULL;
    uint64_t bytes_written = 0;
    if (cfg->use_bin_out) {
        /* Create parent directories for bin_out_path */
        char *path_copy = strdup(cfg->bin_out_path);
        char *parent_dir = dirname(path_copy);
        if (strcmp(parent_dir, ".") != 0 && strcmp(parent_dir, "/") != 0) {
            if (mkdir_p(parent_dir)) {
                fprintf(stderr, "Failed to create parent directory for bin_out: %s\n", parent_dir);
                free(path_copy);
                exit(EXIT_FAILURE);
            }
        }
        free(path_copy);
        
        bin_fp = fopen(cfg->bin_out_path, "wb");
        if (!bin_fp) {
            perror("bin_out");
            exit(EXIT_FAILURE);
        }
        /* Write 9-byte placeholder header: 8-byte uint64 length + 1-byte probe count */
        uint64_t placeholder = 0;
        uint8_t probe_count_placeholder = 0;
        if (fwrite(&placeholder, sizeof(uint64_t), 1, bin_fp) != 1 ||
            fwrite(&probe_count_placeholder, sizeof(uint8_t), 1, bin_fp) != 1) {
            perror("bin_out header write");
            exit(EXIT_FAILURE);
        }
        bytes_written = 9;
    }

    int nprobes = probe_fa ? probe_fa->number_of_features : 0;
    shard_data *shd = shard_data_new(nprobes);
    /* seen_reads guard removed: we rely on primary-only filter and TripKey dedup */
    /* (no separate hash table for read IDs anymore) */
    khash_t(strptr) *read_map = NULL;
    if (cfg->save_read_to_cb)
        read_map = kh_init(strptr);

    bam1_t *b = bam_init1();
    long long recs=0;
    while (sam_read1(in, hdr, b) >= 0) {
        if (cfg->max_records && recs>=cfg->max_records) break;
        recs++;
        
        /* Read CBUB data for this record FIRST to keep streams synchronized */
        const char *cb = NULL;
        const char *ub = NULL;
        int has_tags = 1;
        char umi_buf[65]; /* enough for umi_len<=64, overwritten each read */
        uint32_t cb_idx = 0;
        const char *gene = NULL;
        
        if (cfg->use_cbub) {
            if (!cbub_reader_next(cbub, &has_tags, &cb_idx, umi_buf)) {
                fprintf(stderr, "CBUB stream ended early at record %llu\n", (unsigned long long)cbub->record_idx);
                exit(EXIT_FAILURE);
            }
        }
        
        /* Now apply BAM-level filters */
        const char *qname = bam_get_qname(b);
        uint16_t flag = b->core.flag;
        if (cfg->primary_only && (flag & (BAM_FSECONDARY|BAM_FSUPPLEMENTARY))) {
            continue; 
        }
        if (cfg->skip_dup && (flag & BAM_FDUP)) { 
            continue; 
        }
        if (b->core.qual < cfg->min_mapq) { 
            continue; 
        }

        /* Read NH tag (number of hits) */
        int32_t meta_nh = 1;  /* Default to 1 if NH tag missing */
        const uint8_t *nh_tag = bam_aux_get(b, "NH");
        if (nh_tag) {
            meta_nh = bam_aux2i(nh_tag);
            if (meta_nh < 1) meta_nh = 1;  /* Sanity check */
        }

        /* Binary output: write byte for each record that passes BAM-level filters */
        uint8_t output_byte = 0;

        /* Tag acquisition - branch on CBUB vs BAM tags */
        int valid_cb = 1;  /* Track if CB is in whitelist for later filtering */
        if (cfg->use_cbub) {
            if (!has_tags) {
                if (cfg->use_bin_out) {
                    if (fwrite(&output_byte, 1, 1, bin_fp) != 1) {
                        perror("bin_out write"); exit(EXIT_FAILURE);
                    }
                    bytes_written++;
                }
                continue; /* mirror skip for missing tags */
            }
            /* Check if CB index is valid, but don't skip yet - still do probe lookup */
            if (cb_idx == 0 || cb_idx > whitelist->count) {
                valid_cb = 0;
            }
            /* Convert CB index to actual barcode string (1-based indexing) */
            if (valid_cb) {
                cb = whitelist->barcodes[cb_idx - 1];
            } else {
                /* Use a placeholder for out-of-whitelist CBs - won't be counted */
                cb = "-";
            }
            ub = umi_buf;
            
            /* Still need gene tag from BAM */
            const uint8_t *gx_tag = bam_aux_get(b, cfg->gene_tag);
            if (!gx_tag) gx_tag = bam_aux_get(b, FALLBACK_GENE_TAG);
            if (!gx_tag) {
                if (cfg->use_bin_out) {
                    if (fwrite(&output_byte, 1, 1, bin_fp) != 1) {
                        perror("bin_out write"); exit(EXIT_FAILURE);
                    }
                    bytes_written++;
                }
                continue;
            }
            gene = bam_aux2Z(gx_tag);
            if (!gene) {
                if (cfg->use_bin_out) {
                    if (fwrite(&output_byte, 1, 1, bin_fp) != 1) {
                        perror("bin_out write"); exit(EXIT_FAILURE);
                    }
                    bytes_written++;
                }
                continue;
            }
        } else {
            /* existing bam_aux_get logic */
            const uint8_t *cb_tag = bam_aux_get(b, cfg->cb_tag);
            const uint8_t *ub_tag = bam_aux_get(b, cfg->ub_tag);
            const uint8_t *gx_tag = bam_aux_get(b, cfg->gene_tag);
            if (!gx_tag) gx_tag = bam_aux_get(b, FALLBACK_GENE_TAG);
            if (!cb_tag || !ub_tag || !gx_tag) {
                if (cfg->use_bin_out) {
                    if (fwrite(&output_byte, 1, 1, bin_fp) != 1) {
                        perror("bin_out write"); exit(EXIT_FAILURE);
                    }
                    bytes_written++;
                }
                continue; 
            }
            cb = bam_aux2Z(cb_tag);
            ub = bam_aux2Z(ub_tag);
            gene = bam_aux2Z(gx_tag);
            if (!cb || !ub || !gene) {
                if (cfg->use_bin_out) {
                    if (fwrite(&output_byte, 1, 1, bin_fp) != 1) {
                        perror("bin_out write"); exit(EXIT_FAILURE);
                    }
                    bytes_written++;
                }
                continue;
            }
            if (strcmp(cb, "-") == 0 || strcmp(ub, "-") == 0) {
                valid_cb = 0; /* Mark as invalid for counting, but continue to probe lookup */
            }
        }
        
        /* Check intergenic - if not counting, mark invalid but still do probe lookup */
        int count_this_record = 1;
        if (gene[0] == '-' && !cfg->count_intergene) {
            count_this_record = 0;
        }
        /* Also don't count if CB is invalid (out of whitelist or "-") */
        if (!valid_cb) {
            count_this_record = 0;
        }

        /* packed UMI stored in 'out' but no longer used for dedup */

        /* sample probe lookup */
        int sample_idx = 0;
        if (probe_fa) {
            char kmer[PROBE_LEN + 1];
            /* try primary offset first */
            if (extract_kmer_from_bam(b, cfg->probe_offset, kmer)) {
                int idx = feature_lookup_kmer(kmer, PROBE_LEN, probe_fa, direct_probe);
                if (idx > 0 && idx <= probe_fa->number_of_features)
                    sample_idx = idx;
            }
            /* fallback search if enabled and not found */
            if (sample_idx == 0 && cfg->search_nearby) {
                int offs[4] = { cfg->probe_offset + 1,
                                cfg->probe_offset + 2,
                                cfg->probe_offset - 1,
                                cfg->probe_offset - 2 };
                for (int oi = 0; oi < 4 && sample_idx == 0; ++oi) {
                    if (!extract_kmer_from_bam(b, offs[oi], kmer))
                        continue; /* invalid offset or non-ACGT base */
                    int idx = feature_lookup_kmer(kmer, PROBE_LEN, probe_fa, direct_probe);
                    if (idx > 0 && idx <= probe_fa->number_of_features)
                        sample_idx = idx;   /* early exit once found */
                }
            }
        }

        /* Staging log: log meta-level info right after tag extraction and probe lookup (where meta.nh is set) */
        if (nh_staging_fp && nh_staging_count < NH_STAGING_MAX) {
            char qname_escaped[512];
            json_escape_qname(qname, qname_escaped, sizeof(qname_escaped));
            uint32_t cb_idx_plus1 = (cfg->use_cbub && valid_cb && cb_idx > 0) ? cb_idx : 0;
            fprintf(nh_staging_fp, "{\"event\":\"nh_staging\",\"index\":%lld,\"qname\":\"%s\",\"meta_nh\":%d,\"sampleIdx\":%d,\"cbIdxPlus1\":%u}\n",
                    nh_staging_count, qname_escaped, meta_nh, sample_idx, cb_idx_plus1);
            nh_staging_count++;
        }

        /* Write binary output byte - always write for all reads that reach here */
        if (cfg->use_bin_out) {
            output_byte = (sample_idx > 0) ? (uint8_t)sample_idx : 0;
            if (fwrite(&output_byte, 1, 1, bin_fp) != 1) {
                perror("bin_out write");
                exit(EXIT_FAILURE);
            }
            bytes_written++;
        }

        /* Only count records with valid CB in whitelist and valid tags */
        if (count_this_record && sample_idx > 0) {
            uint32_t cb_id = intern_get(cb_tab, cb);
            
            /* Tap log: log final cand.nh value in CRKeyAggregator::onRecord equivalent */
            if (nh_tap_fp) {
                char qname_escaped[512];
                json_escape_qname(qname, qname_escaped, sizeof(qname_escaped));
                uint32_t cb_idx_plus1 = (cfg->use_cbub && valid_cb && cb_idx > 0) ? cb_idx : 0;
                fprintf(nh_tap_fp, "{\"event\":\"nh_tap\",\"qname\":\"%s\",\"nh\":%d,\"sampleIdx\":%d,\"cbIdxPlus1\":%u}\n",
                        qname_escaped, meta_nh, sample_idx, cb_idx_plus1);
            }
            
            /* get / create per-CB counts array */
            khint_t kcb = kh_get(u32ptr, shd->cb_counts, cb_id);
            uint32_t *arr;
            if (kcb == kh_end(shd->cb_counts)){
                /* allocate from global pools if available, else calloc */
                #ifdef CB_COUNTS_BLOCK_SIZE
                if (global_pools && global_pools->cb_counts_pool){
                    arr = allocate_memory_from_pool(global_pools->cb_counts_pool);
                } else {
                    arr = calloc(nprobes, sizeof(uint32_t));
                }
                #else
                arr = calloc(nprobes, sizeof(uint32_t));
                #endif
                int ret;
                khint_t kh = kh_put(u32ptr, shd->cb_counts, cb_id, &ret);
                kh_val(shd->cb_counts, kh) = arr;
            } else {
                arr = (uint32_t*)kh_val(shd->cb_counts, kcb);
            }
            if (arr[sample_idx-1] < UINT32_MAX)
                arr[sample_idx-1]++;

            shd->usable_reads++;
            shd->probe_tot[sample_idx-1]++;
            
            if (cfg->save_read_to_cb){
                char val_buf[512];
                snprintf(val_buf, sizeof(val_buf), "%s\t%s\t%s", cb, ub, gene);
                char *qname_copy = strdup(qname);
                char *val_copy = strdup(val_buf);
                int ret;
                khint_t k = kh_put(strptr, read_map, qname_copy, &ret);
                kh_val(read_map, k) = val_copy;
            }
        }
        shd->total_reads++;
    }

    /* TripKey dedupulation disabled (Stage 3). Stats simplified */
    fprintf(stderr, "records=%lu usable=%lu\n",
            (unsigned long)shd->total_reads,
            (unsigned long)shd->usable_reads);

    if (cfg->use_cbub) {
        if (cbub->record_idx != cbub->record_count) {
            fprintf(stderr, "Warning: CBUB stream has %llu leftover records\n",
                    (unsigned long long)(cbub->record_count - cbub->record_idx));
        }
        cbub_reader_close(cbub);
        free_whitelist(whitelist);
    }

    /* Close debug log files */
    if (nh_staging_fp) {
        fclose(nh_staging_fp);
        fprintf(stderr, "NH staging log: wrote %lld entries to %s\n", nh_staging_count, nh_staging_path);
    }
    if (nh_tap_fp) {
        fclose(nh_tap_fp);
        fprintf(stderr, "NH tap log: wrote entries to %s\n", nh_tap_path);
    }

    /* Finalize binary output header */
    if (cfg->use_bin_out) {
        uint8_t probe_count = (uint8_t)nprobes;
        /* Seek to beginning and write actual header */
        if (fseek(bin_fp, 0, SEEK_SET) != 0) {
            perror("bin_out seek");
            exit(EXIT_FAILURE);
        }
        /* Convert to little-endian and write 8-byte total length including header */
        uint64_t length_le = htole64(bytes_written);
        if (fwrite(&length_le, sizeof(uint64_t), 1, bin_fp) != 1) {
            perror("bin_out header finalize");
            exit(EXIT_FAILURE);
        }
        /* Write probe count */
        if (fwrite(&probe_count, sizeof(uint8_t), 1, bin_fp) != 1) {
            perror("bin_out header finalize");
            exit(EXIT_FAILURE);
        }
        fclose(bin_fp);
        fprintf(stderr, "Binary output: wrote %lu bytes (%lu records, %u probes)\n",
                (unsigned long)bytes_written, (unsigned long)(bytes_written - 9), probe_count);
    }

    /* Write outputs (skip MatrixMarket if binary mode) */
    if (!cfg->use_bin_out) {
        if (mkdir_p(cfg->outdir)) { fprintf(stderr,"Failed to create outdir\n"); exit(EXIT_FAILURE);}    
        char path[FILENAME_LENGTH];

        /* ---- Write barcodes.tsv & build CB→column map ---- */
        khash_t(u32u32) *cb_to_col = kh_init(u32u32);
        snprintf(path,sizeof(path), "%s/barcodes.tsv", cfg->outdir);
        FILE *fbc = fopen(path,"w");
        uint32_t col=0;
        khint_t k;
        for (k = kh_begin(shd->cb_counts); k != kh_end(shd->cb_counts); ++k) {
            if (!kh_exist(shd->cb_counts, k)) continue;
            uint32_t cb_id = kh_key(shd->cb_counts, k);
            fprintf(fbc, "%s\n", id_to_str(cb_tab, cb_id));
            int ret;
            khint_t kh = kh_put(u32u32, cb_to_col, cb_id, &ret);
            kh_val(cb_to_col, kh) = ++col;
        }
        fclose(fbc);

        /* features.tsv = probe names */
        snprintf(path,sizeof(path), "%s/features.tsv", cfg->outdir);
        FILE *ffe = fopen(path,"w");
        for (int i=0;i<nprobes;i++) fprintf(ffe, "%s\n", probe_fa->feature_names[i]);
        fclose(ffe);

        /* ---- Compute NNZ ---- */
        uint64_t nnz=0;
        for (k = kh_begin(shd->cb_counts); k != kh_end(shd->cb_counts); ++k) {
            if (!kh_exist(shd->cb_counts, k)) continue;
            uint32_t *arr = (uint32_t*)kh_val(shd->cb_counts, k);
            for (int p=0;p<nprobes;p++) if (arr[p]) nnz++;
        }

        /* matrix.mtx */
        snprintf(path,sizeof(path), "%s/matrix.mtx", cfg->outdir);
        FILE *fmtx = fopen(path,"w");
        fprintf(fmtx, "%%MatrixMarket matrix coordinate integer general\n");
        fprintf(fmtx, "%d %u %llu\n", nprobes, col, (unsigned long long)nnz);

        for (k = kh_begin(shd->cb_counts); k != kh_end(shd->cb_counts); ++k) {
            if (!kh_exist(shd->cb_counts, k)) continue;
            uint32_t cb_id = kh_key(shd->cb_counts, k);
            khint_t kcol = kh_get(u32u32, cb_to_col, cb_id);
            uint32_t col_idx = (kcol != kh_end(cb_to_col)) ? kh_val(cb_to_col, kcol) : 0;
            uint32_t *arr = (uint32_t*)kh_val(shd->cb_counts, k);
            for (int p=0;p<nprobes;p++)
                if (arr[p]) fprintf(fmtx, "%d %u %u\n", p+1, col_idx, arr[p]);
        }
        fclose(fmtx);
        kh_destroy(u32u32, cb_to_col);

        /* stats.txt */
        snprintf(path,sizeof(path), "%s/stats.txt", cfg->outdir);
        FILE *fst = fopen(path,"w");
        fprintf(fst,"Total_records\t%lu\n", (unsigned long)shd->total_reads);
        fprintf(fst,"Usable_records\t%lu\n", (unsigned long)shd->usable_reads);
        fprintf(fst,"BC_ID\tTotal_reads\n");
        for(int i=0;i<nprobes;i++){
            fprintf(fst,"%s\t%lu\n", probe_fa->feature_names[i],
                    (unsigned long)shd->probe_tot[i]);
        }
        fclose(fst);
    }

    /* Clean */
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(in);

    /* probe_single removed */

    /* write read_to_cb (works in both binary and MatrixMarket modes) */
    if (cfg->save_read_to_cb){
        if (mkdir_p(cfg->outdir)) { fprintf(stderr,"Failed to create outdir\n"); exit(EXIT_FAILURE);}
        char path[FILENAME_LENGTH];
        snprintf(path, sizeof(path), "%s/read_to_cb_umi_gene.txt", cfg->outdir);
        FILE *fp = fopen(path, "w");
        if (!fp) { perror("read_to_cb"); exit(EXIT_FAILURE);} 
        khint_t k;
        for (k = kh_begin(read_map); k != kh_end(read_map); ++k) {
            if (kh_exist(read_map, k))
                fprintf(fp, "%s\t%s\n", (char*)kh_key(read_map, k), (char*)kh_val(read_map, k));
        }
        fclose(fp);
        // Free keys and values
        for (k = kh_begin(read_map); k != kh_end(read_map); ++k) {
            if (kh_exist(read_map, k)) {
                free((char*)kh_key(read_map, k));
                free((char*)kh_val(read_map, k));
            }
        }
        kh_destroy(strptr, read_map);
    }
}

int main(int argc, char **argv) {
    cfg_t cfg = {0};
    cfg.primary_only = 1; cfg.skip_dup = 1; cfg.min_mapq = 0; cfg.hts_threads=2; cfg.consumer_threads=1;
    cfg.probe_offset = 68; cfg.direct_probe = -1;
    parse_cli(&cfg, argc, argv);

    /* init matching tables */
    barcode_match_init();
    feature_code_hash = kh_init(codeu32);

    if (cfg.sample_probes) {
        probe_fa     = load_probe_variants_to_features(cfg.sample_probes);
        probe_offset = cfg.probe_offset;
        direct_probe = 0;          /* disable direct search for now */
    }

    /* Validate probe count for binary output mode */
    if (cfg.use_bin_out) {
        if (!probe_fa) {
            fprintf(stderr, "Error: --bin_out requires --sample_probes to be specified\n");
            exit(EXIT_FAILURE);
        }
        if (probe_fa->number_of_features > 255) {
            fprintf(stderr, "Error: --bin_out requires probe count <= 255, got %d\n", 
                    probe_fa->number_of_features);
            exit(EXIT_FAILURE);
        }
    }

    /* ---- Stage 4: set up memory pools for CB×probe counts ---- */
    dynamic_struct_sizes.cb_probe_counts = (probe_fa?probe_fa->number_of_features:1) * sizeof(uint32_t);
    global_pools = initialize_memory_pool_collection();

    intern_table *cb_tab = intern_table_new();
    intern_table *gene_tab = intern_table_new();

    process_bam_single(&cfg, cb_tab, gene_tab);

    if (global_pools) free_memory_pool_collection(global_pools);
    return 0;
}

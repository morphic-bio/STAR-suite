#include "../include/common.h"
#include "../include/io.h"
#include "../include/utils.h"
#include "../include/barcode_match.h"

#define PROBE_LEN 8
static pthread_mutex_t sinks_map_mutex;
static pthread_mutex_t sample_map_mutex;
typedef struct {
  gzFile r1, r2, r3;
  pthread_mutex_t lock;    // NEW: per-sink write lock
} demux_sink;

typedef struct {
  feature_arrays *features;        // loaded probe variants (names = BC ids)
  GHashTable *libprobe_to_sample;  // key: "LIBID\tBCxxx"; val: sample name
  GHashTable *sample_to_sink;      // key: sample name; val: demux_sink*
  char outdir[FILENAME_LENGTH];
  char undetermined[64];
  char unspecified[64];
  int probe_read_idx;  // 0=R1(barcode),1=R2(forward),2=R3(reverse)
  int probe_offset;    // 0-based
  long long max_records; // 0 = no limit
  int direct_search;   // -1=auto (default), 0=hash, 1=direct compare
  int threads;         // number of consumer threads (default 1)
} demux_cfg;

// Use shared implementations from barcode_match.c

// Simple size-based ordering without stat() to avoid heavy deps
#if 0   /* duplicate – now provided in io.c */
void sort_samples_by_size(fastq_files_collection *fastq_files, int *sample_order) {
  for (int i = 0; i < fastq_files->nsamples; ++i) sample_order[i] = i;
  for (int i = 0; i < fastq_files->nsamples; ++i) {
    for (int j = i + 1; j < fastq_files->nsamples; ++j) {
      if (fastq_files->sample_sizes[sample_order[j]] > fastq_files->sample_sizes[sample_order[i]]) {
        int t = sample_order[i]; sample_order[i] = sample_order[j]; sample_order[j] = t;
      }
    }
  }
}
#endif


static inline void chomp(char *s) {
  size_t n = strlen(s);
  while (n && (s[n-1] == '\n' || s[n-1] == '\r')) s[--n] = 0;
}

static inline int read_fastq_record(gzFile f, char *l1, char *l2, char *l3, char *l4, size_t buflen) {
  if (!f) return 0;
  if (!gzgets(f, l1, buflen)) return 0;
  if (!gzgets(f, l2, buflen)) return -1;
  if (!gzgets(f, l3, buflen)) return -1;
  if (!gzgets(f, l4, buflen)) return -1;
  return 1;
}

static gboolean extract_kmer_at_offset(const char *seq, int offset, char out8[PROBE_LEN+1]) {
  const size_t n = strlen(seq);
  if ((int)n < offset + PROBE_LEN) return FALSE;
  for (int i = 0; i < PROBE_LEN; ++i) {
    char c = seq[offset + i];
    if (c >= 'a' && c <= 'z') c -= 32;
    if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return FALSE;
    out8[i] = c;
  }
  out8[PROBE_LEN] = 0;
  return TRUE;
}

static feature_arrays* load_probe_variants_to_features(const char *path) {
  FILE *fp = fopen(path, "r");
  if (!fp) { perror("probe_barcodes"); exit(EXIT_FAILURE); }

  // First pass: count lines and size requirements
  int count = 0;
  int name_size = 0, seq_size = 0, code_size = 0;
  int maxFeatureLength = PROBE_LEN;
  char line[LINE_LENGTH];
  while (fgets(line, sizeof(line), fp)) {
    char *save = NULL;
    char *v = strtok_r(line, "\t \r\n", &save);   // variant
    (void)strtok_r(NULL, "\t \r\n", &save);        // skip original
    char *bc = strtok_r(NULL, "\t \r\n", &save);   // BCxxx
    if (!v || !bc) continue;
    if ((int)strlen(v) != PROBE_LEN) continue;
    if (!check_sequence(v, PROBE_LEN)) continue;
    count++;
    name_size += (int)strlen(bc) + 1;
    seq_size  += PROBE_LEN + 1;
    code_size += (PROBE_LEN + 3) / 4;
  }

  feature_arrays *fa = allocate_feature_arrays(name_size, seq_size, code_size, count, maxFeatureLength);
  fa->common_length = PROBE_LEN;

  // Second pass: fill arrays and populate feature_code_hash
  fseek(fp, 0, SEEK_SET);
  int idx = 0;
  fa->feature_names[0] = fa->feature_names_storage;
  fa->feature_sequences[0] = fa->feature_sequences_storage;
  fa->feature_codes[0] = fa->feature_codes_storage;

  while (fgets(line, sizeof(line), fp)) {
    char *save = NULL;
    char *v = strtok_r(line, "\t \r\n", &save);
    (void)strtok_r(NULL, "\t \r\n", &save);
    char *bc = strtok_r(NULL, "\t \r\n", &save);
    if (!v || !bc) continue;
    if ((int)strlen(v) != PROBE_LEN) continue;
    if (!check_sequence(v, PROBE_LEN)) continue;

    // Copy name and sequence
    strcpy(fa->feature_names[idx], bc);
    if (idx + 1 < fa->number_of_features) {
      fa->feature_names[idx + 1] = fa->feature_names[idx] + strlen(bc) + 1;
    }
    strcpy(fa->feature_sequences[idx], v);
    fa->feature_lengths[idx] = PROBE_LEN;

    // Code and insert into feature_code_hash
    fa->feature_code_lengths[idx] = string2code(v, PROBE_LEN, fa->feature_codes[idx]);
    if (fa->feature_lengths[idx] == fa->common_length) {
      GBytes *key = g_bytes_new_static(fa->feature_codes[idx], fa->feature_code_lengths[idx]);
      g_hash_table_insert(feature_code_hash, key, GUINT_TO_POINTER(idx + 1));
    }

    if (idx + 1 < fa->number_of_features) {
      fa->feature_sequences[idx + 1] = fa->feature_sequences[idx] + PROBE_LEN + 1;
      fa->feature_codes[idx + 1] = fa->feature_codes[idx] + fa->feature_code_lengths[idx];
    }
    idx++;
  }
  fclose(fp);
  fprintf(stderr, "Loaded %d probe variants into feature arrays (len=%d)\n", idx, fa->common_length);
  return fa;
}

static gboolean parse_library_id(const char *basename, char *out, size_t outsz) {
  const char *p = strstr(basename, "SC");
  if (!p) return FALSE;
  size_t i = 2;
  while (isdigit((unsigned char)p[i])) i++;
  if (i == 2) return FALSE;
  size_t len = i;
  if (len >= outsz) return FALSE;
  memcpy(out, p, len);
  out[len] = 0;
  return TRUE;
}

static G_GNUC_UNUSED GHashTable* load_probe_variants(const char *path) {
  FILE *fp = fopen(path, "r");
  if (!fp) { perror("probe_barcodes"); exit(EXIT_FAILURE); }
  GHashTable *h = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
  char line[LINE_LENGTH];
  while (fgets(line, sizeof(line), fp)) {
    chomp(line);
    if (!line[0]) continue;
    char *a = strtok(line, "\t ");
    strtok(NULL, "\t "); // skip 2nd col
    char *c = strtok(NULL, "\t ");
    if (!a || !c) continue;
    if (strlen(a) != PROBE_LEN) continue;
    g_hash_table_replace(h, g_strdup(a), g_strdup(c)); // variant -> BCxxx
  }
  fclose(fp);
  fprintf(stderr, "Loaded %u probe variants\n", g_hash_table_size(h));
  return h;
}

static GHashTable* load_sample_map(const char *path) {
  FILE *fp = fopen(path, "r");
  if (!fp) { perror("sample_map"); exit(EXIT_FAILURE); }
  GHashTable *h = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
  char line[LINE_LENGTH];
  int first = 1;
  while (fgets(line, sizeof(line), fp)) {
    chomp(line);
    if (!line[0]) continue;
    if (first) { first = 0; continue; } // skip header
    char *lib = strtok(line, "\t ");
    char *bc  = strtok(NULL, "\t ");
    char *sm  = strtok(NULL, "\t ");
    if (!lib || !bc || !sm) continue;
    char key[NAME_LENGTH*2];
    snprintf(key, sizeof(key), "%s\t%s", lib, bc);
    g_hash_table_replace(h, g_strdup(key), g_strdup(sm));
  }
  fclose(fp);
  fprintf(stderr, "Loaded %u lib+probe → sample mappings\n", g_hash_table_size(h));
  return h;
}

static void free_sink(gpointer data) {
  demux_sink *s = (demux_sink*)data;
  if (!s) return;
  if (s->r1) gzclose(s->r1);
  if (s->r2) gzclose(s->r2);
  if (s->r3) gzclose(s->r3);
  pthread_mutex_destroy(&s->lock);    // NEW
  free(s);
}

static demux_sink* get_or_open_sink(demux_cfg *cfg, const char *sample,
                                    const char *base_r1, const char *base_r2, const char *base_r3) {
  pthread_mutex_lock(&sinks_map_mutex);            // NEW: protect map
  demux_sink *s = (demux_sink*)g_hash_table_lookup(cfg->sample_to_sink, sample);
  if (s) { pthread_mutex_unlock(&sinks_map_mutex); return s; }

  char dir[FILENAME_LENGTH*2];
  snprintf(dir, sizeof(dir), "%s/%s", cfg->outdir, sample);
  if (mkdir_p(dir)) { fprintf(stderr, "Failed to create %s\n", dir); exit(EXIT_FAILURE); }

  s = (demux_sink*)calloc(1, sizeof(demux_sink));
  pthread_mutex_init(&s->lock, NULL);             // NEW

  if (base_r1) { char p[FILENAME_LENGTH*3]; snprintf(p, sizeof(p), "%s/%s", dir, base_r1);
                 s->r1 = gzopen(p, "wb"); if (!s->r1) { perror("gzopen r1"); exit(EXIT_FAILURE); } }
  if (base_r2) { char p[FILENAME_LENGTH*3]; snprintf(p, sizeof(p), "%s/%s", dir, base_r2);
                 s->r2 = gzopen(p, "wb"); if (!s->r2) { perror("gzopen r2"); exit(EXIT_FAILURE); } }
  if (base_r3) { char p[FILENAME_LENGTH*3]; snprintf(p, sizeof(p), "%s/%s", dir, base_r3);
                 s->r3 = gzopen(p, "wb"); if (!s->r3) { perror("gzopen r3"); exit(EXIT_FAILURE); } }

  g_hash_table_replace(cfg->sample_to_sink, g_strdup(sample), s);
  pthread_mutex_unlock(&sinks_map_mutex);         // NEW
  return s;
}

static void write_record(demux_sink *s,
                         const char *r1_1, const char *r1_2, const char *r1_3, const char *r1_4,
                         const char *r2_1, const char *r2_2, const char *r2_3, const char *r2_4,
                         const char *r3_1, const char *r3_2, const char *r3_3, const char *r3_4) {
  if (s->r1) { gzputs(s->r1, r1_1); gzputs(s->r1, r1_2); gzputs(s->r1, r1_3); gzputs(s->r1, r1_4); }
  if (s->r2) { gzputs(s->r2, r2_1); gzputs(s->r2, r2_2); gzputs(s->r2, r2_3); gzputs(s->r2, r2_4); }
  if (s->r3) { gzputs(s->r3, r3_1); gzputs(s->r3, r3_2); gzputs(s->r3, r3_3); gzputs(s->r3, r3_4); }
}

static G_GNUC_UNUSED const char* resolve_sample(demux_cfg *cfg, const char *libid, const char *probe_id) {
  if (!probe_id || !libid) return cfg->undetermined;
  char key[NAME_LENGTH*2];
  snprintf(key, sizeof(key), "%s\t%s", libid, probe_id);
  char *s = (char*)g_hash_table_lookup(cfg->libprobe_to_sample, key);
  return s ? s : cfg->undetermined;
}

static int path_is_directory(const char *p) {
  struct stat st;
  if (stat(p, &st) != 0) return 0;
  return S_ISDIR(st.st_mode) ? 1 : 0;
}

static void usage(const char *prog) {
  fprintf(stderr,
    "Usage: %s --probe_barcodes <file> --sample_map <file> --outdir <dir> [input FASTQs or dirs]\n"
    "Options:\n"
    "  --barcode_fastqs STR     Comma-separated R1 paths\n"
    "  --forward_fastqs STR     Comma-separated R2 paths\n"
    "  --reverse_fastqs STR     Comma-separated R3 paths\n"
    "  --barcode_fastq_pattern STR   Default _R1_\n"
    "  --forward_fastq_pattern STR   Default _R2_\n"
    "  --reverse_fastq_pattern STR   Default _R3_\n"
    "  --probe_read {R1|R2|R3}  Default R2\n"
    "  --probe_offset INT       0-based offset into probe read (required)\n"
    "  --max_records N          Limit for testing (default 0 = all)\n"
    "  --direct_search          Force direct 64-bit compare "
    "(default auto: direct unless features > 128)\n"
    "  --threads N              Number of consumer threads (default 1)\n"
    "  -v                       Debug\n",
    prog);
}

// Multithreaded FASTQ block reader for demux_fastq
typedef struct {
  gzFile r1, r2, r3;
  int have_r1, have_r2, have_r3;

  char   **buffer;
  char    *buffer_storage;
  size_t   read_buffer_lines;   // total slots in ring (must be multiple of lines_per_block)
  size_t   produce_index;
  size_t   consume_index;
  size_t   filled;              // number of occupied slots (lines)
  pthread_mutex_t mutex;
  pthread_cond_t  can_produce;
  pthread_cond_t  can_consume;
  int done;                     // producer done
  int error;                    // producer error
} dmx_reader_set;



static inline int dmx_lines_per_block(const dmx_reader_set *s) {
  return 4 * ((s->have_r1 ? 1 : 0) + (s->have_r2 ? 1 : 0) + (s->have_r3 ? 1 : 0));
}

static dmx_reader_set* dmx_alloc_reader_set(gzFile r1, gzFile r2, gzFile r3, size_t blocks) {
  dmx_reader_set *s = (dmx_reader_set*)calloc(1, sizeof(dmx_reader_set));
  s->r1 = r1; s->r2 = r2; s->r3 = r3;
  s->have_r1 = (r1 != NULL);
  s->have_r2 = (r2 != NULL);
  s->have_r3 = (r3 != NULL);
  const int lpblock = dmx_lines_per_block(s);
  if (lpblock <= 0) { free(s); return NULL; }
  s->read_buffer_lines = lpblock * blocks;
  s->buffer_storage = (char*)malloc(s->read_buffer_lines * LINE_LENGTH);
  s->buffer = (char**)malloc(s->read_buffer_lines * sizeof(char*));
  for (size_t i = 0; i < s->read_buffer_lines; ++i)
    s->buffer[i] = s->buffer_storage + i * LINE_LENGTH;

  pthread_mutex_init(&s->mutex, NULL);
  pthread_cond_init(&s->can_produce, NULL);
  pthread_cond_init(&s->can_consume, NULL);
  return s;
}

static void dmx_free_reader_set(dmx_reader_set *s) {
  if (!s) return;
  pthread_mutex_destroy(&s->mutex);
  pthread_cond_destroy(&s->can_produce);
  pthread_cond_destroy(&s->can_consume);
  free(s->buffer);
  free(s->buffer_storage);
  free(s);
}

static void* dmx_producer(void *arg) {
  dmx_reader_set *s = (dmx_reader_set*)arg;
  const int lpblock = dmx_lines_per_block(s);

  char r1_1[LINE_LENGTH], r1_2[LINE_LENGTH], r1_3[LINE_LENGTH], r1_4[LINE_LENGTH];
  char r2_1[LINE_LENGTH], r2_2[LINE_LENGTH], r2_3[LINE_LENGTH], r2_4[LINE_LENGTH];
  char r3_1[LINE_LENGTH], r3_2[LINE_LENGTH], r3_3[LINE_LENGTH], r3_4[LINE_LENGTH];

  for (;;) {
    int ok1 = 1, ok2 = 1, ok3 = 1;

    if (s->have_r1) ok1 = read_fastq_record(s->r1, r1_1, r1_2, r1_3, r1_4, sizeof(r1_1));
    if (s->have_r2) ok2 = read_fastq_record(s->r2, r2_1, r2_2, r2_3, r2_4, sizeof(r2_1));
    if (s->have_r3) ok3 = read_fastq_record(s->r3, r3_1, r3_2, r3_3, r3_4, sizeof(r3_1));

    // EOF handling: all present must hit EOF together
    int present = (s->have_r1?1:0) + (s->have_r2?1:0) + (s->have_r3?1:0);
    int n_ok = (s->have_r1 && ok1 > 0) + (s->have_r2 && ok2 > 0) + (s->have_r3 && ok3 > 0);
    int n_eof = (s->have_r1 && ok1 == 0) + (s->have_r2 && ok2 == 0) + (s->have_r3 && ok3 == 0);

    if ((s->have_r1 && ok1 < 0) || (s->have_r2 && ok2 < 0) || (s->have_r3 && ok3 < 0)) {
      fprintf(stderr, "Error: FASTQ read error in producer thread\n");
      s->error = 1;
      break;
    }
    if (n_eof > 0) {
      if (n_eof != present) {
        fprintf(stderr, "FASTQ length mismatch across R1/R2/R3\n");
        s->error = 1;
      }
      break; // EOF
    }
    if (n_ok != present) {
      fprintf(stderr, "FASTQ length mismatch across R1/R2/R3\n");
      s->error = 1;
      break;
    }

    // Enqueue a full block
    pthread_mutex_lock(&s->mutex);
    while (s->filled > s->read_buffer_lines - (size_t)lpblock)
      pthread_cond_wait(&s->can_produce, &s->mutex);

    size_t p = s->produce_index;
    if (s->have_r1) {
      strcpy(s->buffer[p], r1_1); p = (p+1) % s->read_buffer_lines;
      strcpy(s->buffer[p], r1_2); p = (p+1) % s->read_buffer_lines;
      strcpy(s->buffer[p], r1_3); p = (p+1) % s->read_buffer_lines;
      strcpy(s->buffer[p], r1_4); p = (p+1) % s->read_buffer_lines;
    }
    if (s->have_r2) {
      strcpy(s->buffer[p], r2_1); p = (p+1) % s->read_buffer_lines;
      strcpy(s->buffer[p], r2_2); p = (p+1) % s->read_buffer_lines;
      strcpy(s->buffer[p], r2_3); p = (p+1) % s->read_buffer_lines;
      strcpy(s->buffer[p], r2_4); p = (p+1) % s->read_buffer_lines;
    }
    if (s->have_r3) {
      strcpy(s->buffer[p], r3_1); p = (p+1) % s->read_buffer_lines;
      strcpy(s->buffer[p], r3_2); p = (p+1) % s->read_buffer_lines;
      strcpy(s->buffer[p], r3_3); p = (p+1) % s->read_buffer_lines;
      strcpy(s->buffer[p], r3_4); p = (p+1) % s->read_buffer_lines;
    }

    s->produce_index = (s->produce_index + (size_t)lpblock) % s->read_buffer_lines;
    s->filled += (size_t)lpblock;
    pthread_cond_signal(&s->can_consume);
    pthread_mutex_unlock(&s->mutex);
  }

  if (s->r1) gzclose(s->r1);
  if (s->r2) gzclose(s->r2);
  if (s->r3) gzclose(s->r3);

  pthread_mutex_lock(&s->mutex);
  s->done = 1;
  pthread_cond_broadcast(&s->can_consume);
  pthread_mutex_unlock(&s->mutex);

  return NULL;
}

typedef struct {
  dmx_reader_set *dmx;
  demux_cfg *cfg;
  const char *base_r1, *base_r2, *base_r3;
  char libid[64];
  int lpblock;
  long long total, matched, unspecified_cnt, undetermined_cnt;
} dmx_consumer_args;

static void* dmx_consumer(void *arg) {
  dmx_consumer_args *a = (dmx_consumer_args*)arg;
  dmx_reader_set *dmx = a->dmx;
  demux_cfg *cfg = a->cfg;
  const int lp = a->lpblock;

  // local copies of FASTQ lines
  char r1_1[LINE_LENGTH], r1_2[LINE_LENGTH], r1_3[LINE_LENGTH], r1_4[LINE_LENGTH];
  char r2_1[LINE_LENGTH], r2_2[LINE_LENGTH], r2_3[LINE_LENGTH], r2_4[LINE_LENGTH];
  char r3_1[LINE_LENGTH], r3_2[LINE_LENGTH], r3_3[LINE_LENGTH], r3_4[LINE_LENGTH];

  for (;;) {
    // get a full block
    pthread_mutex_lock(&dmx->mutex);
    while (dmx->filled < (size_t)lp && !dmx->done)
      pthread_cond_wait(&dmx->can_consume, &dmx->mutex);
    if (dmx->filled < (size_t)lp) { pthread_mutex_unlock(&dmx->mutex); break; }

    size_t c = dmx->consume_index;
    if (dmx->have_r1) { strcpy(r1_1, dmx->buffer[c]); c=(c+1)%dmx->read_buffer_lines;
                        strcpy(r1_2, dmx->buffer[c]); c=(c+1)%dmx->read_buffer_lines;
                        strcpy(r1_3, dmx->buffer[c]); c=(c+1)%dmx->read_buffer_lines;
                        strcpy(r1_4, dmx->buffer[c]); c=(c+1)%dmx->read_buffer_lines; }
    if (dmx->have_r2) { strcpy(r2_1, dmx->buffer[c]); c=(c+1)%dmx->read_buffer_lines;
                        strcpy(r2_2, dmx->buffer[c]); c=(c+1)%dmx->read_buffer_lines;
                        strcpy(r2_3, dmx->buffer[c]); c=(c+1)%dmx->read_buffer_lines;
                        strcpy(r2_4, dmx->buffer[c]); c=(c+1)%dmx->read_buffer_lines; }
    if (dmx->have_r3) { strcpy(r3_1, dmx->buffer[c]); c=(c+1)%dmx->read_buffer_lines;
                        strcpy(r3_2, dmx->buffer[c]); c=(c+1)%dmx->read_buffer_lines;
                        strcpy(r3_3, dmx->buffer[c]); c=(c+1)%dmx->read_buffer_lines;
                        strcpy(r3_4, dmx->buffer[c]); c=(c+1)%dmx->read_buffer_lines; }

    dmx->consume_index = (dmx->consume_index + (size_t)lp) % dmx->read_buffer_lines;
    dmx->filled -= (size_t)lp;
    pthread_cond_signal(&dmx->can_produce);
    pthread_mutex_unlock(&dmx->mutex);

    if (dmx->error) break;

    // pick probe sequence based on probe_read selection
    const char *probe_seq = (cfg->probe_read_idx == 0 ? r1_2 :
                             cfg->probe_read_idx == 1 ? r2_2 : r3_2);
    if (!probe_seq) break;
    chomp((char*)probe_seq);

    char kmer[PROBE_LEN+1];
    const char *probe_id = NULL;
    if (extract_kmer_at_offset(probe_seq, cfg->probe_offset, kmer)) {
      int index1 = feature_lookup_kmer(kmer, PROBE_LEN, cfg->features, cfg->direct_search);
      if (index1 > 0 && index1 <= cfg->features->number_of_features) {
        probe_id = cfg->features->feature_names[index1 - 1];
      }
    }

    // Resolve sample group (lock only the sample_map for lookup)
    char sample_group[FILENAME_LENGTH];
    if (!probe_id) {
      snprintf(sample_group, sizeof(sample_group), "%s/%s_unk", cfg->undetermined, a->libid);
      a->undetermined_cnt++;
    } else {
      char key[NAME_LENGTH*2];
      snprintf(key, sizeof(key), "%s\t%s", a->libid, probe_id);
      pthread_mutex_lock(&sample_map_mutex);
      const char *mapped_sample = (const char*)g_hash_table_lookup(cfg->libprobe_to_sample, key);
      pthread_mutex_unlock(&sample_map_mutex);
      if (mapped_sample) { snprintf(sample_group, sizeof(sample_group), "%s", mapped_sample); a->matched++; }
      else { snprintf(sample_group, sizeof(sample_group), "%s/%s_%s", cfg->unspecified, a->libid, probe_id); a->unspecified_cnt++; }
    }

    // Get sink (map protected inside), then lock only that sink to write
    demux_sink *sink = get_or_open_sink(cfg, sample_group, a->base_r1, a->base_r2, a->base_r3);
    pthread_mutex_lock(&sink->lock);
    write_record(sink,
                 dmx->have_r1 ? r1_1 : "", dmx->have_r1 ? r1_2 : "", dmx->have_r1 ? r1_3 : "", dmx->have_r1 ? r1_4 : "",
                 dmx->have_r2 ? r2_1 : "", dmx->have_r2 ? r2_2 : "", dmx->have_r2 ? r2_3 : "", dmx->have_r2 ? r2_4 : "",
                 dmx->have_r3 ? r3_1 : "", dmx->have_r3 ? r3_2 : "", dmx->have_r3 ? r3_3 : "", dmx->have_r3 ? r3_4 : "");
    pthread_mutex_unlock(&sink->lock);

    a->total++;
    if (cfg->max_records && a->total >= cfg->max_records) break;
  }
  return NULL;
}

int main(int argc, char **argv) {
  demux_cfg cfg;
  memset(&cfg, 0, sizeof(cfg));
  strcpy(cfg.outdir, ".");
  strcpy(cfg.undetermined, "undetermined");
  strcpy(cfg.unspecified, "unspecified");
  cfg.probe_read_idx = 1; // default R2
  cfg.probe_offset = -1;
  cfg.direct_search = -1; // auto by default
  cfg.threads = 1;        // default 1 consumer

  char barcode_pattern[LINE_LENGTH] = "_R1_";
  char forward_pattern[LINE_LENGTH] = "_R2_";
  char reverse_pattern[LINE_LENGTH] = "_R3_";
  char *barcodeFastqFilesString = NULL;
  char *forwardFastqFilesString = NULL;
  char *reverseFastqFilesString = NULL;
  const char *probe_barcodes_path = NULL;
  const char *sample_map_path = NULL;
  int debug_flag = 0; (void)debug_flag;

  static struct option long_options[] = {
    {"probe_barcodes", required_argument, 0, 'p'},
    {"sample_map", required_argument, 0, 's'},
    {"outdir", required_argument, 0, 'o'},
    {"probe_read", required_argument, 0, 'r'},
    {"probe_offset", required_argument, 0, 'O'},
    {"max_records", required_argument, 0, 'M'},
    {"barcode_fastqs", required_argument, 0, 100},
    {"forward_fastqs", required_argument, 0, 101},
    {"reverse_fastqs", required_argument, 0, 102},
    {"barcode_fastq_pattern", required_argument, 0, 200},
    {"forward_fastq_pattern", required_argument, 0, 201},
    {"reverse_fastq_pattern", required_argument, 0, 202},
    {"direct_search", no_argument, 0, 300},
    {"threads", required_argument, 0, 't'},
    {"debug", no_argument, 0, 'v'},
    {0,0,0,0}
  };

  int opt, idx = 0;
  while ((opt = getopt_long(argc, argv, "p:s:o:r:O:M:t:v", long_options, &idx)) != -1) {
    switch (opt) {
      case 'p': probe_barcodes_path = optarg; break;
      case 's': sample_map_path = optarg; break;
      case 'o': g_strlcpy(cfg.outdir, optarg, sizeof(cfg.outdir)); break;
      case 'r':
        if      (!strcmp(optarg, "R1")) cfg.probe_read_idx = 0;
        else if (!strcmp(optarg, "R2")) cfg.probe_read_idx = 1;
        else if (!strcmp(optarg, "R3")) cfg.probe_read_idx = 2;
        else { fprintf(stderr, "Invalid --probe_read %s\n", optarg); return 1; }
        break;
      case 'O': cfg.probe_offset = atoi(optarg); break;
      case 'M': cfg.max_records = atoll(optarg); break;
      case 't': cfg.threads = atoi(optarg); if (cfg.threads < 1) cfg.threads = 1; break;
      case 'v': debug_flag = 1; break;
      case 100: barcodeFastqFilesString = strdup(optarg); break;
      case 101: forwardFastqFilesString = strdup(optarg); break;
      case 102: reverseFastqFilesString = strdup(optarg); break;
      case 200: g_strlcpy(barcode_pattern,  optarg, sizeof(barcode_pattern));  break;
      case 201: g_strlcpy(forward_pattern,  optarg, sizeof(forward_pattern));  break;
      case 202: g_strlcpy(reverse_pattern,  optarg, sizeof(reverse_pattern));  break;
      case 300: cfg.direct_search = 1; break;
      default: usage(argv[0]); return 1;
    }
  }
  if (!probe_barcodes_path || !sample_map_path || cfg.probe_offset < 0) {
    usage(argv[0]); return 1;
  }
  if (cfg.outdir[strlen(cfg.outdir)-1] == '/') cfg.outdir[strlen(cfg.outdir)-1] = 0;
  if (mkdir_p(cfg.outdir)) { fprintf(stderr, "Failed to create outdir %s\n", cfg.outdir); return 1; }

  // Initialise core matching tables
  barcode_match_init();
  initialize_complement();
  feature_code_hash = g_hash_table_new_full(g_bytes_hash, g_bytes_equal, (GDestroyNotify)g_bytes_unref, NULL);

  cfg.features = load_probe_variants_to_features(probe_barcodes_path);
  cfg.libprobe_to_sample = load_sample_map(sample_map_path);
  cfg.sample_to_sink = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, free_sink);

  fastq_files_collection fastq_files;
  memset(&fastq_files, 0, sizeof(fastq_files));

  int positional_arg_count = argc - optind;
  if (positional_arg_count && path_is_directory(argv[optind])) {
    organize_fastq_files_by_directory(positional_arg_count, argc, argv, optind,
      barcodeFastqFilesString, forwardFastqFilesString, reverseFastqFilesString,
      &fastq_files, barcode_pattern, forward_pattern, reverse_pattern);
  } else {
    organize_fastq_files_by_type(positional_arg_count, argc, argv, optind,
      barcodeFastqFilesString, forwardFastqFilesString, reverseFastqFilesString,
      &fastq_files, barcode_pattern, forward_pattern, reverse_pattern, 1);
  }
  check_filecounts(&fastq_files);

  for (int j = 0; j < fastq_files.nbarcode_files; ++j) {
    const char *r1_path = fastq_files.barcode_fastq ? fastq_files.barcode_fastq[j] : NULL;
    const char *r2_path = (fastq_files.nforward_files ? fastq_files.forward_fastq[j] : NULL);
    const char *r3_path = (fastq_files.nreverse_files ? fastq_files.reverse_fastq[j] : NULL);

    const char *probe_path = (cfg.probe_read_idx == 0 ? r1_path :
                              cfg.probe_read_idx == 1 ? r2_path : r3_path);
    if (!probe_path) { fprintf(stderr, "Probe read missing for set %d\n", j); return 1; }

    const char *base_r1 = r1_path ? get_basename(r1_path) : NULL;
    const char *base_r2 = r2_path ? get_basename(r2_path) : NULL;
    const char *base_r3 = r3_path ? get_basename(r3_path) : NULL;

    char libid[64] = {0};
    if (!parse_library_id(get_basename(probe_path), libid, sizeof(libid))) {
      fprintf(stderr, "Failed to parse library_id from %s\n", get_basename(probe_path));
      return 1;
    }

    gzFile r1 = r1_path ? gzopen(r1_path, "rb") : NULL;
    gzFile r2 = r2_path ? gzopen(r2_path, "rb") : NULL;
    gzFile r3 = r3_path ? gzopen(r3_path, "rb") : NULL;
    if ((r1_path && !r1) || (r2_path && !r2) || (r3_path && !r3)) { perror("gzopen input"); return 1; }

    // Start producer
    dmx_reader_set *dmx = dmx_alloc_reader_set(r1, r2, r3, /*blocks*/ 256);
    if (!dmx) { fprintf(stderr, "Failed to init demux reader\n"); return 1; }
    pthread_t prod;
    if (pthread_create(&prod, NULL, dmx_producer, dmx) != 0) { perror("pthread_create"); return 1; }
    const int lpblock = dmx_lines_per_block(dmx);

    long long total = 0, matched = 0, unspecified_cnt = 0, undetermined_cnt = 0;

    // Launch consumer threads
    const int ncons = (cfg.threads > 0 ? cfg.threads : 1);
    pthread_t *cons = (pthread_t*)malloc((size_t)ncons * sizeof(pthread_t));
    dmx_consumer_args *cargs = (dmx_consumer_args*)calloc((size_t)ncons, sizeof(dmx_consumer_args));
    for (int i = 0; i < ncons; ++i) {
      cargs[i].dmx = dmx;
      cargs[i].cfg = &cfg;
      cargs[i].base_r1 = base_r1;
      cargs[i].base_r2 = base_r2;
      cargs[i].base_r3 = base_r3;
      g_strlcpy(cargs[i].libid, libid, sizeof(cargs[i].libid));
      cargs[i].lpblock = lpblock;
      pthread_create(&cons[i], NULL, dmx_consumer, &cargs[i]);
    }

    // Wait for consumers to finish
    for (int i = 0; i < ncons; ++i) {
      pthread_join(cons[i], NULL);
      total          += cargs[i].total;
      matched        += cargs[i].matched;
      unspecified_cnt+= cargs[i].unspecified_cnt;
      undetermined_cnt+=cargs[i].undetermined_cnt;
    }
    free(cons);
    free(cargs);

    // Join producer and cleanup ring
    pthread_join(prod, NULL);
    dmx_free_reader_set(dmx);

    fprintf(stderr, "Set %d: lib=%s total=%lld matched=%lld unspecified=%lld undetermined=%lld\n",
            j, libid, total, matched, unspecified_cnt, undetermined_cnt);

    g_hash_table_remove_all(cfg.sample_to_sink); // closes sinks for this set
  }

  g_hash_table_destroy(cfg.sample_to_sink);
  free_feature_arrays(cfg.features);
  g_hash_table_destroy(cfg.libprobe_to_sample);
  g_hash_table_destroy(feature_code_hash);
  free_fastq_files_collection(&fastq_files);
  if (barcodeFastqFilesString) free(barcodeFastqFilesString);
  if (forwardFastqFilesString) free(forwardFastqFilesString);
  if (reverseFastqFilesString) free(reverseFastqFilesString);
  return 0;
}

#ifndef PF_BGZF_INPUT_H
#define PF_BGZF_INPUT_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum { PF_BGZF_AUTO = 0, PF_BGZF_OFF = 1, PF_BGZF_RANGE = 2 } pf_bgzf_mode;
typedef struct pf_bgzf_input pf_bgzf_input;
typedef struct {
    const char *name, *sequence, *quality;
    size_t name_length, sequence_length, quality_length;
    uint64_t ordinal;
} pf_bgzf_record;
typedef struct {
    void *context;
    uint64_t (*acquire)(void *);
    void (*release)(void *, uint64_t, uint64_t, uint64_t, uint64_t);
} pf_bgzf_permits;

/* Returns 1 for regular BGZF, 0 for a legacy input, -1 for an I/O error.
 * Non-regular paths are never opened by detection (e.g. FIFOs). */
int pf_bgzf_detect(const char *path, char *error, size_t error_size);
/* Stream order is barcode, forward (if present), reverse (if present).
 * worker_threads is a TOTAL across this lane's streams; zero inflates on the
 * calling producer. No workers or buffers depend on the STAR executable. */
pf_bgzf_input *pf_bgzf_open(const char *const *paths, int streams,
                           unsigned worker_threads, int check_crc,
                           const pf_bgzf_permits *permits,
                           char *error, size_t error_size);
/* Returns 1 for a tuple, 0 at joint EOF, -1 on error. Views remain valid until
 * the next call or close. Names use the first token, ignoring terminal /1,/2,/3.
 * Every stream must have the same name and logical record ordinal. */
int pf_bgzf_next(pf_bgzf_input *input, pf_bgzf_record *records,
                 char *error, size_t error_size);
void pf_bgzf_close(pf_bgzf_input *input);
typedef int (*pf_bgzf_batch_consumer)(void *, unsigned, const pf_bgzf_record *, size_t, int);
/* Bounded leased batches, shared assignment workers, stable lane-local ordinals.
 * Callback owns a worker ID exclusively and must return held permits before it
 * returns. No view survives the callback. max_reads applies per lane. */
int pf_bgzf_process_batches(const char *const *paths, unsigned lanes, int streams,
    unsigned workers, unsigned inflater_threads, int crc, uint64_t max_reads,
    const pf_bgzf_permits *permits, pf_bgzf_batch_consumer consume, void *context,
    char *error, size_t error_size);

#ifdef __cplusplus
}
#endif
#endif

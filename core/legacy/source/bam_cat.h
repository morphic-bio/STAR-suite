#ifndef CODE_bam_cat
#define CODE_bam_cat

#if defined(WITH_CHROMAP) && WITH_CHROMAP
#include <htslib/sam.h>
#else
#include "htslib/htslib/sam.h"
#endif

int bam_cat(int nfn, char * const *fn, const bam_hdr_t *h, const char* outbam);

#endif

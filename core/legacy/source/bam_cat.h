#ifndef CODE_bam_cat
#define CODE_bam_cat

#if defined(STAR_EXTERNAL_HTSLIB) && STAR_EXTERNAL_HTSLIB
#include <htslib/sam.h>
#else
#include "htslib/htslib/sam.h"
#endif

int bam_cat(int nfn, char * const *fn, const bam_hdr_t *h, const char* outbam);

#endif

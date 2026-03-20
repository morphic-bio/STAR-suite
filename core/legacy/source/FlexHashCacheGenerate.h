#ifndef H_FlexHashCacheGenerate
#define H_FlexHashCacheGenerate

class Parameters;
class Genome;
class Transcriptome;

namespace libem {
class Transcriptome;
}

void runFlexHashCacheGenerate(Parameters& P, Genome& genome, Transcriptome* transcriptomeMain,
                              const libem::Transcriptome* libemTr);

#endif

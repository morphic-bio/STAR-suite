#ifndef H_ZGZXTags
#define H_ZGZXTags

#include "IncludeDefine.h"
#include "ReadAnnotations.h"
#include "Transcriptome.h"

class ZGZXTags {
public:
    // Format ZG tag (gene set) from ReadAnnotFeature
    static string formatZGTag(const ReadAnnotFeature &annFeat, const Transcriptome &transcriptome);
    
    // Format ZX tag (overlap status) from ReadAnnotFeature
    static string formatZXTag(const ReadAnnotFeature &annFeat);
    
private:
    // Convert overlap type enum to string
    static const char* overlapTypeToString(uint32 ovType);
};

#endif

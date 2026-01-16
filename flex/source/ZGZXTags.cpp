#include "ZGZXTags.h"
#include <sstream>
#include <algorithm>
#include <unordered_set>

string ZGZXTags::formatZGTag(const ReadAnnotFeature &annFeat, const Transcriptome &transcriptome) {
    if (annFeat.fSet.empty()) {
        return "-";
    }

    const vector<string> &canonicalIDs = transcriptome.geIDCanonical.empty()
        ? transcriptome.geID
        : transcriptome.geIDCanonical;

    std::unordered_set<string> seen;
    seen.reserve(annFeat.fSet.size());

    std::ostringstream oss;
    bool first = true;

    for (uint32 geneIdx : annFeat.fSet) {
        string geneId;
        if (geneIdx < canonicalIDs.size()) {
            geneId = canonicalIDs[geneIdx];
            if (geneId.empty() && geneIdx < transcriptome.geID.size()) {
                geneId = transcriptome.geID[geneIdx];
            }
        } else {
            geneId = "UNKNOWN_GENE_" + to_string(geneIdx);
        }

        if (geneId.empty()) {
            continue;
        }

        if (!seen.insert(geneId).second) {
            continue; // avoid emitting duplicates when multiple indices map to same gene
        }

        if (!first) {
            oss << ",";
        }
        oss << geneId;
        first = false;
    }

    if (first) {
        return "-";
    }

    return oss.str();
}

string ZGZXTags::formatZXTag(const ReadAnnotFeature &annFeat) {
    return overlapTypeToString(annFeat.ovType);
}

const char* ZGZXTags::overlapTypeToString(uint32 ovType) {
    switch (ovType) {
        case ReadAnnotFeature::overlapTypes::none:
            return "none";
        case ReadAnnotFeature::overlapTypes::exonic:
            return "exonic";
        case ReadAnnotFeature::overlapTypes::exonicAS:
            return "exonic";  // Treat antisense exonic as exonic for simplicity
        case ReadAnnotFeature::overlapTypes::exonic50p:
            return "exonic";  // Treat 50% exonic as exonic
        case ReadAnnotFeature::overlapTypes::exonic50pAS:
            return "exonic";  // Treat antisense 50% exonic as exonic
        case ReadAnnotFeature::overlapTypes::intronic:
            return "intronic";
        case ReadAnnotFeature::overlapTypes::intronicAS:
            return "intronic";  // Treat antisense intronic as intronic
        case ReadAnnotFeature::overlapTypes::intergenic:
            return "none";  // Treat intergenic as none
        default:
            return "spanning";  // Default for unknown types
    }
}

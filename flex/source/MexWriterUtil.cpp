#include "MexWriterUtil.h"

namespace MexWriterUtil {

int writeMexFromDedup(
    const std::string& outputPrefix,
    const std::vector<std::string>& barcodes,
    const std::vector<std::string>& geneIds,
    const std::vector<MexWriter::Triplet>& triplets)
{
    return writeMexFromDedup(outputPrefix, barcodes, geneIds, triplets, 1u);
}

int writeMexFromDedup(
    const std::string& outputPrefix,
    const std::vector<std::string>& barcodes,
    const std::vector<std::string>& geneIds,
    const std::vector<MexWriter::Triplet>& triplets,
    unsigned int matrixThreads)
{
    // Convert gene IDs to Feature structs with default settings
    std::vector<MexWriter::Feature> features;
    features.reserve(geneIds.size());
    
    for (const auto& geneId : geneIds) {
        features.emplace_back(geneId, geneId, "Gene Expression");
    }
    
    // Delegate to MexWriter
    return MexWriter::writeMex(
        outputPrefix, barcodes, features, triplets, -1, matrixThreads);
}

int writeMexFromDedup(
    const std::string& outputPrefix,
    const std::vector<std::string>& barcodes,
    const std::vector<MexWriter::Feature>& features,
    const std::vector<MexWriter::Triplet>& triplets)
{
    return writeMexFromDedup(outputPrefix, barcodes, features, triplets, 1u);
}

int writeMexFromDedup(
    const std::string& outputPrefix,
    const std::vector<std::string>& barcodes,
    const std::vector<MexWriter::Feature>& features,
    const std::vector<MexWriter::Triplet>& triplets,
    unsigned int matrixThreads)
{
    // Direct delegation to MexWriter
    return MexWriter::writeMex(
        outputPrefix, barcodes, features, triplets, -1, matrixThreads);
}

} // namespace MexWriterUtil

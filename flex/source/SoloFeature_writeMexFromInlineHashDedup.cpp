#include "SoloFeature.h"
#include "MexWriterUtil.h"
#include "Parameters.h"
#include "TimeFunctions.h"
#include "SampleDetector.h"
#include <sstream>
#include <iomanip>
#include <map>
#include <fstream>
#include <algorithm>
#include <stdexcept>

/**
 * @brief Build composite MEX matrix data from inline hash dedup counts
 *
 * Populates composite barcodes, gene list, sparse triplets, and SampleMatrixData.
 */
SoloFeature::InlineMatrixBundle SoloFeature::buildInlineMatrixFromHash(
    const std::unordered_map<uint64_t, std::vector<std::pair<uint32_t, uint32_t>>>& cbTagGeneCounts)
{
    InlineMatrixBundle bundle;

    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) 
                    << " ... Building in-memory MEX from inline-hash dedup data" << endl;

    if (cbTagGeneCounts.empty()) {
        P.inOut->logMain << "WARNING: No (CB,TAG) data to build MEX in memory" << endl;
        return bundle;
    }

    // Step 1: Build sorted list of (cbIdx, tagIdx) pairs and create composite barcodes
    std::vector<uint64_t> cbTagKeys;
    cbTagKeys.reserve(cbTagGeneCounts.size());
    for (const auto &entry : cbTagGeneCounts) {
        cbTagKeys.push_back(entry.first);
    }
    std::sort(cbTagKeys.begin(), cbTagKeys.end(), [](uint64_t a, uint64_t b) {
        uint32_t cbA = a >> 8, cbB = b >> 8;
        uint8_t tagA = a & 0xFF, tagB = b & 0xFF;
        if (cbA != cbB) return cbA < cbB;
        return tagA < tagB;
    });

    std::vector<std::string> compositeBarcodes;
    std::map<uint64_t, uint32_t> cbTagToColIdx;
    compositeBarcodes.reserve(cbTagKeys.size());
    uint32_t colIdx = 0;

    for (uint64_t cbTagKey : cbTagKeys) {
        uint32_t cbIdx = cbTagKey >> 8;
        uint8_t tagIdx = cbTagKey & 0xFF;
        if (tagIdx == 0) {
            continue;
        }

        if (cbIdx >= pSolo.cbWLstr.size()) {
            std::ostringstream errMsg;
            errMsg << "ERROR: Cell barcode index " << cbIdx 
                   << " out of range (whitelist size: " << pSolo.cbWLstr.size() << ")";
            P.inOut->logMain << errMsg.str() << endl;
            throw std::runtime_error(errMsg.str());
        }

        if (tagIdx >= gCanonicalTags.size() || gCanonicalTags[tagIdx].empty()) {
            std::ostringstream errMsg;
            errMsg << "ERROR: Missing canonical tag sequence for tagIdx=" << (int)tagIdx;
            P.inOut->logMain << errMsg.str() << endl;
            throw std::runtime_error(errMsg.str());
        }

        std::string composite = pSolo.cbWLstr[cbIdx] + gCanonicalTags[tagIdx];
        compositeBarcodes.push_back(composite);
        cbTagToColIdx[cbTagKey] = colIdx;
        colIdx++;
    }

    P.inOut->logMain << "  Unique (CB,TAG) combinations: " << compositeBarcodes.size() << endl;

    // Step 2: Build gene list (ordered by probe index)
    std::vector<std::string> geneIds;
    uint32_t maxProbeIdx = 0;
    for (const auto& entry : cbTagGeneCounts) {
        for (const auto& gc : entry.second) {
            if (gc.first > maxProbeIdx) maxProbeIdx = gc.first;
        }
    }
    bool usedProbeList = false;
    if (!P.pSolo.probeListPath.empty() && P.pSolo.probeListPath != "-") {
        std::ifstream probeFile(P.pSolo.probeListPath);
        if (probeFile.is_open()) {
            std::string line;
            while (std::getline(probeFile, line)) {
                if (!line.empty() && line[0] == '#') continue;
                size_t beg = line.find_first_not_of(" \t\r\n");
                size_t end = line.find_last_not_of(" \t\r\n");
                if (beg == std::string::npos) continue;
                geneIds.push_back(line.substr(beg, end - beg + 1));
            }
            probeFile.close();
            if (maxProbeIdx <= geneIds.size()) {
                usedProbeList = true;
            } else {
                P.inOut->logMain << "WARNING: probe list too small for max probe index " << maxProbeIdx
                                 << " (size " << geneIds.size() << "), falling back to transcriptome features" << endl;
                geneIds.clear();
            }
        } else {
            P.inOut->logMain << "WARNING: could not open probe list at " << P.pSolo.probeListPath
                             << ", falling back to transcriptome features" << endl;
        }
    }
    if (!usedProbeList) {
        std::ostringstream err;
        err << "ERROR: Probe list unavailable or too small for max probe index " << maxProbeIdx
            << " (probe list path: " << P.pSolo.probeListPath << ")";
        P.inOut->logMain << err.str() << endl;
        throw std::runtime_error(err.str());
    }

    // Step 3: Build sparse triplets
    bundle.triplets.reserve(compositeBarcodes.size() * 4);
    for (const auto& entry : cbTagGeneCounts) {
        auto colIt = cbTagToColIdx.find(entry.first);
        if (colIt == cbTagToColIdx.end()) {
            continue;
        }
        uint32_t cellIdx = colIt->second;
        for (const auto& geneCount : entry.second) {
            uint32_t probeIdx = geneCount.first;
            if (probeIdx == 0 || probeIdx > geneIds.size()) {
                continue;
            }
            bundle.triplets.push_back({cellIdx, probeIdx - 1, geneCount.second});
        }
    }

    if (bundle.triplets.empty()) {
        P.inOut->logMain << "WARNING: No non-zero entries after dedup (zero triplets)" << endl;
    }

    std::sort(bundle.triplets.begin(), bundle.triplets.end(), [](const MexWriter::Triplet& a, const MexWriter::Triplet& b) {
        return a.cell_idx < b.cell_idx;
    });

    // Step 4: Populate matrix data
    SampleMatrixData &matrix = bundle.matrixData;
    matrix.nCells = compositeBarcodes.size();
    matrix.nGenes = geneIds.size();
    matrix.countMatStride = 3;
    matrix.barcodes = compositeBarcodes;
    matrix.features = geneIds;
    matrix.nUMIperCB.assign(matrix.nCells, 0);
    matrix.nGenePerCB.assign(matrix.nCells, 0);
    matrix.countCellGeneUMIindex.assign(matrix.nCells + 1, 0);

    uint32_t offset = 0;
    size_t tripletIdx = 0;
    for (uint32_t cell = 0; cell < matrix.nCells; ++cell) {
        matrix.countCellGeneUMIindex[cell] = offset;
        while (tripletIdx < bundle.triplets.size() && bundle.triplets[tripletIdx].cell_idx == cell) {
            matrix.countCellGeneUMI.push_back(bundle.triplets[tripletIdx].gene_idx);
            matrix.countCellGeneUMI.push_back(bundle.triplets[tripletIdx].count);
            matrix.countCellGeneUMI.push_back(0);
            matrix.nUMIperCB[cell] += bundle.triplets[tripletIdx].count;
            matrix.nGenePerCB[cell] += 1;
            offset += matrix.countMatStride;
            ++tripletIdx;
        }
    }
    matrix.countCellGeneUMIindex[matrix.nCells] = offset;

    P.inOut->logMain << "  Genes: " << matrix.nGenes << ", Entries: " << bundle.triplets.size() << endl;
    return bundle;
}

/**
 * @brief Write composite MEX from precalculated inline matrix bundle.
 */
void SoloFeature::writeMexFromInlineHashDedup(
    const std::string& outputPrefix,
    const InlineMatrixBundle& bundle)
{
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) 
                     << " ... Writing MEX from inline-hash dedup data" << endl;

    if (bundle.matrixData.nCells == 0 || bundle.matrixData.features.empty()) {
        P.inOut->logMain << "WARNING: Inline matrix bundle is empty, skipping MEX write" << endl;
        return;
    }

    int result = MexWriterUtil::writeMexFromDedup(
        outputPrefix,
        bundle.matrixData.barcodes,
        bundle.matrixData.features,
        bundle.triplets
    );

    if (result == 0) {
        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime) 
                        << " ... Finished writing MEX format:" << endl;
        std::string matrixPath = outputPrefix + "matrix.mtx";
        std::string barcodesPath = outputPrefix + "barcodes.tsv";
        std::string featuresPath = outputPrefix + "features.tsv";
        P.inOut->logMain << "  " << matrixPath << endl;
        P.inOut->logMain << "  " << barcodesPath << endl;
        P.inOut->logMain << "  " << featuresPath << endl;
        P.inOut->logMain << "  Cells (CB+TAG combos): " << bundle.matrixData.nCells << endl;
        P.inOut->logMain << "  Features: " << bundle.matrixData.nGenes << endl;
        P.inOut->logMain << "  Entries: " << bundle.triplets.size() << endl;
    } else {
        P.inOut->logMain << "ERROR: Failed to write MEX format" << endl;
    }
}

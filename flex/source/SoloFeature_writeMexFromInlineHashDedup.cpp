#include "SoloFeature.h"
#include "MexWriterUtil.h"
#include "Parameters.h"
#include "TimeFunctions.h"
#include <algorithm>
#include <chrono>
#include <fstream>

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

    const auto mexStart = std::chrono::steady_clock::now();
    int result = MexWriterUtil::writeMexFromDedup(
        outputPrefix,
        bundle.matrixData.barcodes,
        bundle.matrixData.features,
        bundle.triplets,
        static_cast<unsigned int>(std::max(1, P.runThreadN))
    );
    P.inOut->logMain << "Solo timing: raw MEX "
                     << std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - mexStart).count()
                     << " s, matrix_threads=" << std::max(1, P.runThreadN)
                     << endl << std::flush;

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

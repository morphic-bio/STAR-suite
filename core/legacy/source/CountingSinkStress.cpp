#include "CountingSinkStress.h"
#include "Parameters.h"
#include "SoloFeature.h"
#include "SoloReadInfoSink.h"
#include "SoloMemoryProfile.h"
#include "Transcriptome.h"
#include "TimeFunctions.h"
#include "systemFunctions.h"
#include <cstdlib>
#include <string>

namespace {

uint64_t envUint64(const char *name, uint64_t defaultVal)
{
    const char *env = std::getenv(name);
    if (env == nullptr || env[0] == '\0') {
        return defaultVal;
    }
    return std::strtoull(env, nullptr, 10);
}

uint32_t envUint32(const char *name, uint32_t defaultVal)
{
    return static_cast<uint32_t>(envUint64(name, defaultVal));
}

} // namespace

void runCountingSinkStress(Parameters &P)
{
    const uint64_t nRecords = envUint64("STAR_COUNTING_SINK_STRESS_NRECORDS", 50000000ULL);
    const uint32_t activeCBs = envUint32("STAR_COUNTING_SINK_STRESS_ACTIVE_CBS", 20000U);
    const uint32_t featureMod = envUint32("STAR_COUNTING_SINK_STRESS_FEATURE_MOD", 5000U);
    const bool runCollapse = (std::getenv("STAR_COUNTING_SINK_STRESS_COLLAPSE") != nullptr);

    if (P.pSolo.cbWLsize == 0) {
        P.inOut->logMain << "EXITING: countingSinkStress requires --soloCBwhitelist (cbWLsize>0)\n" << std::flush;
        std::exit(EXIT_CODE_INPUT_FILES);
    }
    if (activeCBs == 0 || activeCBs > P.pSolo.cbWLsize) {
        P.inOut->logMain << "EXITING: STAR_COUNTING_SINK_STRESS_ACTIVE_CBS must be in 1..cbWLsize ("
                         << P.pSolo.cbWLsize << ")\n"
                         << std::flush;
        std::exit(EXIT_CODE_INPUT_FILES);
    }

    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ..... countingSinkStress: nRecords=" << nRecords
                     << " activeCBs=" << activeCBs
                     << " cbWLsize=" << P.pSolo.cbWLsize
                     << " featureMod=" << featureMod
                     << " collapse=" << (runCollapse ? "yes" : "no")
                     << std::endl;

    Transcriptome trans(P);
    SoloFeature *soloFeatArr[1] = {nullptr};
    soloFeatArr[0] = new SoloFeature(P, nullptr, trans, SoloFeatureTypes::GeneFull, nullptr, soloFeatArr);
    SoloFeature &feat = *soloFeatArr[0];
    feat.nReadsInput = nRecords;
    feat.nReadsMapped = nRecords;

    soloMemoryProfileCheckpoint(P.inOut->logMain, "countingSinkStress_begin");

    CountingSink sink;
    for (uint64_t i = 0; i < nRecords; ++i) {
        ReadInfoRecord rec{};
        rec.readId = i;
        rec.readIndex = static_cast<uint32_t>(i);
        rec.cbIdx = static_cast<uint32_t>(i % activeCBs);
        rec.umi = static_cast<uint32_t>((i >> 8) & 0xFFFFF);
        rec.status = 1;
        rec.featureId = static_cast<uint32_t>(i % featureMod);
        sink.onRecord(feat, rec);
    }

    soloMemoryProfileCheckpoint(P.inOut->logMain, "countingSinkStress_loader_done");
    sink.finalize(feat);

    if (runCollapse) {
        feat.nReadPerCB.resize(feat.nCB);
        feat.nReadPerCBmax = 0;
        for (uint32_t iCB = 0; iCB < feat.nCB; ++iCB) {
            feat.nReadPerCB[iCB] = feat.rCBn[iCB];
            if (feat.nReadPerCB[iCB] > feat.nReadPerCBmax) {
                feat.nReadPerCBmax = feat.nReadPerCB[iCB];
            }
        }
        feat.nUMIperCB.resize(feat.nCB);
        feat.nGenePerCB.resize(feat.nCB);
        feat.countMatStride = feat.pSolo.umiDedup.yes.N + 1;
        feat.countCellGeneUMI.resize(static_cast<size_t>(feat.nReadsMapped * feat.countMatStride / 5 + 16));
        feat.countCellGeneUMIindex.resize(feat.nCB + 1, 0);
        soloMemoryProfileCheckpoint(P.inOut->logMain, "countingSinkStress_collapse_begin");
        feat.collapseUMIall();
        soloMemoryProfileCheckpoint(P.inOut->logMain, "countingSinkStress_collapse_done");
    }

    soloMemoryProfileCheckpoint(P.inOut->logMain, "countingSinkStress_done");
    P.inOut->logMain << "countingSinkStress complete\n" << std::flush;

    delete soloFeatArr[0];
}

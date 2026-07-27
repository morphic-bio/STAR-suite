#include "FlexGdna.h"
#include "FlexHashScreen.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << "\n";
        std::exit(1);
    }
}

bool near(double lhs, double rhs, double tolerance = 1.0e-10)
{
    return std::fabs(lhs - rhs) <= tolerance;
}

std::string makeTempDirectory()
{
    char pattern[] = "/tmp/star_flex_gdna_test.XXXXXX";
    char* path = mkdtemp(pattern);
    require(path != nullptr, "mkdtemp");
    return path;
}

void testPackedValue()
{
    const uint32_t spliced = flexGdnaPackValue(7, FlexGdnaSpliced);
    const uint32_t unspliced = flexGdnaPackValue(11, FlexGdnaUnspliced);
    require(flexGdnaValueCount(spliced) == 7, "spliced count decode");
    require(flexGdnaValueRegion(spliced) == FlexGdnaSpliced,
            "spliced region decode");

    const uint32_t merged = flexGdnaMergeValue(spliced, unspliced);
    require(flexGdnaValueCount(merged) == 18, "merged count");
    require(flexGdnaValueRegion(merged) == FlexGdnaConflicting,
            "opposite regions conflict");

    const uint32_t unknown = flexGdnaPackValue(3, FlexGdnaUnknown);
    const uint32_t known = flexGdnaMergeValue(unknown, spliced);
    require(flexGdnaValueCount(known) == 10, "unknown merge count");
    require(flexGdnaValueRegion(known) == FlexGdnaSpliced,
            "unknown evidence is identity");

    const uint32_t saturated = flexGdnaMergeValue(
        flexGdnaPackValue(kFlexGdnaCountMask, FlexGdnaSpliced),
        flexGdnaPackValue(1, FlexGdnaSpliced));
    require(flexGdnaValueCount(saturated) == kFlexGdnaCountMask,
            "count saturation");
}

void writeFixture(const std::string& directory,
                  std::string& probeListPath,
                  std::string& probeCsvPath)
{
    probeListPath = directory + "/probe_list.txt";
    probeCsvPath = directory + "/filtered_probe_set.csv";

    std::ofstream probeList(probeListPath.c_str());
    std::ofstream probeCsv(probeCsvPath.c_str());
    require(probeList.is_open() && probeCsv.is_open(), "open metadata fixtures");

    probeCsv << "gene_id,probe_id,region,included\n";
    for (unsigned int gene = 1; gene <= 12; ++gene) {
        const std::string geneId = "GENE" + std::to_string(gene);
        probeList << geneId << "\n";
        probeCsv << geneId << "," << geneId << "_S,spliced,true\n";
        probeCsv << geneId << "," << geneId << "_U,unspliced,true\n";
    }
    probeCsv << "GENE1,EXCLUDED,unspliced,false\n";
    probeCsv << "GENE1,deprecated_probe,unspliced,true\n";
}

void testMetadataAndEstimator(const std::string& directory)
{
    std::string probeListPath;
    std::string probeCsvPath;
    writeFixture(directory, probeListPath, probeCsvPath);

    FlexGdnaProbeMetadata& metadata = FlexGdnaProbeMetadata::instance();
    std::string error;
    require(metadata.load(probeCsvPath, probeListPath, &error),
            "metadata load: " + error);
    require(metadata.ready(), "metadata ready");
    require(metadata.totalProbes() == 24, "included probe count");
    require(metadata.totalUnsplicedProbes() == 12,
            "included unspliced probe count");
    require(metadata.controlGeneCount() == 12, "control gene count");
    require(metadata.regionForProbeId("GENE2_S") == FlexGdnaSpliced,
            "spliced probe lookup");
    require(metadata.regionForProbeId("GENE2_U") == FlexGdnaUnspliced,
            "unspliced probe lookup");
    require(metadata.regionForProbeId("missing") == FlexGdnaUnknown,
            "missing probe lookup");
    require(FlexGdnaProbeMetadata::discoverProbeCsv(
                "auto", probeListPath, "/missing") == probeCsvPath,
            "probe CSV discovery");

    std::vector<FlexGdnaGeneMoleculeCounts> molecules(13);
    uint64_t classified = 0;
    for (unsigned int gene = 1; gene <= 12; ++gene) {
        molecules[gene].spliced = gene - 1;
        molecules[gene].unspliced = 1;
        classified += molecules[gene].spliced + molecules[gene].unspliced;
    }
    const FlexGdnaEstimate estimate =
        flexGdnaEstimate(metadata, molecules, classified, 4, 2, 5);
    require(estimate.valid && estimate.status == "ok", "valid estimate");
    require(estimate.controlGenes == 12, "estimate control genes");
    require(estimate.totalFilteredMolecules == 84,
            "gene-assigned filtered molecule audit");
    require(estimate.classifiedMolecules == 78, "classified molecule audit");
    require(estimate.unknownMolecules == 4, "unknown molecule audit");
    require(estimate.conflictingMolecules == 2, "conflicting molecule audit");
    require(estimate.unassignedMolecules == 5, "unassigned molecule audit");
    require(near(estimate.estimatedGdnaPerProbe, 1.0),
            "constant unspliced estimate");
    require(near(estimate.estimatedGdnaFraction, 12.0 / 84.0),
            "gDNA fraction");
    require(estimate.threshold == 1.0, "gDNA threshold");

    std::vector<FlexGdnaGeneMoleculeCounts> tooFew(13);
    metadata.reset();
    require(!flexGdnaEstimate(metadata, tooFew, 1, 0, 0).valid,
            "missing metadata rejection");
}

void testCacheEncoding(const std::string& directory)
{
    FlexHashScreenCache::Record record;
    const std::string sequence(50, 'A');
    require(FlexHashScreenCache::encodeProbeWindow(
                sequence.c_str(), 0, record.seqLo, record.seqHi),
            "probe window encoding");
    record.resolvedGeneIdx15 = 23;
    record.cacheClass = 1;
    record.sampleIdx = 5;
    record.probeRegion = FlexGdnaUnspliced;

    struct CacheHeaderRaw {
        char magic[8];
        uint16_t version;
        uint16_t kmerLength;
        uint32_t recordSize;
        uint64_t recordCount;
    };
    struct CacheRecordRaw {
        uint64_t seqLo;
        uint64_t seqHi;
        uint32_t resolvedGeneIdx15;
        uint8_t cacheClass;
        uint8_t negativeCode;
        uint16_t reserved;
    };

    for (int regionComplete = 0; regionComplete <= 1; ++regionComplete) {
        std::vector<FlexHashScreenCache::Record> records(1, record);
        const std::string path = directory
            + (regionComplete ? "/cache_v3.bin" : "/cache_v2.bin");
        std::string error;
        require(FlexHashScreenCache::writeHashCacheFile(
                    path, records, &error, regionComplete != 0),
                "cache write: " + error);
        std::ifstream in(path.c_str(), std::ios::binary);
        CacheHeaderRaw header {};
        CacheRecordRaw raw {};
        in.read(reinterpret_cast<char*>(&header), sizeof(header));
        in.read(reinterpret_cast<char*>(&raw), sizeof(raw));
        require(in.good(), "cache readback");
        require(header.version == (regionComplete ? 3 : 2),
                "cache version");
        require((raw.resolvedGeneIdx15 & 0x7FFFu) == 23,
                "cache gene encoding");
        const uint32_t encodedRegion = (raw.resolvedGeneIdx15 >> 30) & 0x3u;
        require(encodedRegion == (regionComplete ? FlexGdnaUnspliced
                                                  : FlexGdnaUnknown),
                "cache region encoding");
    }
}

} // namespace

int main()
{
    const std::string directory = makeTempDirectory();
    testPackedValue();
    testMetadataAndEstimator(directory);
    testCacheEncoding(directory);
    std::cout << "Flex gDNA diagnostic tests passed\n";
    return 0;
}

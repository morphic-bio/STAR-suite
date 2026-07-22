#include "SpatialFeatureSidecar.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace sfs = spatial_feature_sidecar;

int main(int argc, char **argv)
{
    try {
        if (argc != 2) {
            throw std::invalid_argument("usage: spatial_feature_sidecar_dump SIDECAR.bin");
        }
        sfs::Reader reader;
        std::string error;
        if (!reader.open(argv[1], error) || !reader.validateAll(error)) {
            throw std::runtime_error(error);
        }
        const sfs::Header &header = reader.header();
        std::cout << "#schema_version=" << header.schemaVersion << '\n'
                  << "#header_bytes=" << header.headerBytes << '\n'
                  << "#record_bytes=" << header.recordBytes << '\n'
                  << "#complete=" << (header.complete ? 1 : 0) << '\n'
                  << "#total_reads=" << header.totalReads << '\n'
                  << "#feature_count=" << header.featureCount << '\n'
                  << "#gene_dictionary_sha256=" << header.geneDictionarySha256 << '\n'
                  << "ordinal\tgene_index\tstatus_flags\treserved\n";
        for (std::uint64_t ordinal = 0; ordinal < header.totalReads; ++ordinal) {
            sfs::Record record;
            if (!reader.read(ordinal, record, error)) throw std::runtime_error(error);
            std::cout << ordinal << '\t' << record.geneIndex << '\t'
                      << record.statusFlags << '\t' << record.reserved << '\n';
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "spatial_feature_sidecar_dump: ERROR: " << error.what() << '\n';
        return 2;
    }
}

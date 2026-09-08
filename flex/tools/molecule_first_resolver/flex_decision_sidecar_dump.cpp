#include "FlexDecisionSidecar.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fds = flex_decision_sidecar;

int main(int argc, char **argv)
{
    try {
        if (argc != 2) {
            throw std::invalid_argument(
                "usage: flex_decision_sidecar_dump SIDECAR.bin");
        }
        fds::Reader reader;
        std::string error;
        if (!reader.open(argv[1], error) || !reader.validateAll(error)) {
            throw std::runtime_error(error);
        }
        std::cout
            << "global_ordinal\tlane\tlane_ordinal\tqname_fnv1a64"
               "\tcache_action\tcache_class\tmatched_cache_class"
               "\tsingle_n_attempted\tsingle_n_resolved\tgene_idx15"
               "\tnegative_code\thash_offset\tsample_checked\tsample_matched"
               "\tsample_rejected\tsample_token\talignment_handoff"
               "\talignment_ran\talignment_resolved\talignment_rejected"
               "\talignment_source\tno_align_dropped\tfinal_reason\n";
        for (std::uint64_t ordinal = 0; ordinal < reader.header().totalReads;
             ++ordinal) {
            fds::Record record;
            if (!reader.read(ordinal, record, error)) {
                throw std::runtime_error(error);
            }
            std::cout << ordinal << '\t';
            if (record.laneIndex == fds::kMissingLane) std::cout << '.';
            else std::cout << record.laneIndex;
            std::cout << '\t';
            if (record.laneOrdinal == fds::kMissingLaneOrdinal) std::cout << '.';
            else std::cout << record.laneOrdinal;
            std::cout << '\t' << std::hex << std::setw(16) << std::setfill('0')
                      << record.qnameHash << std::dec << std::setfill(' ')
                      << '\t' << fds::cacheActionName(record.cacheAction)
                      << '\t' << fds::cacheClassName(record.cacheClass)
                      << '\t' << fds::cacheClassName(record.matchedCacheClass)
                      << '\t' << ((record.statusFlags & fds::kSingleNAttempted) != 0)
                      << '\t' << ((record.statusFlags & fds::kSingleNResolved) != 0)
                      << '\t' << record.geneIdx15
                      << '\t' << static_cast<unsigned>(record.negativeCode)
                      << '\t' << static_cast<int>(record.hashOffset)
                      << '\t' << ((record.statusFlags & fds::kSampleChecked) != 0)
                      << '\t' << ((record.statusFlags & fds::kSampleMatched) != 0)
                      << '\t' << ((record.statusFlags & fds::kSampleRejected) != 0)
                      << '\t';
            if (record.sampleToken == 0xFF) std::cout << '.';
            else std::cout << static_cast<unsigned>(record.sampleToken);
            const char *alignmentSource = ".";
            if (record.statusFlags & fds::kAlignmentProbe) alignmentSource = "PROBE";
            else if (record.statusFlags & fds::kAlignmentGenomic) alignmentSource = "GENOMIC";
            std::cout << '\t' << ((record.statusFlags & fds::kAlignmentHandoff) != 0)
                      << '\t' << ((record.statusFlags & fds::kAlignmentRan) != 0)
                      << '\t' << ((record.statusFlags & fds::kAlignmentResolved) != 0)
                      << '\t' << ((record.statusFlags & fds::kAlignmentRejected) != 0)
                      << '\t' << alignmentSource
                      << '\t' << ((record.statusFlags & fds::kNoAlignDropped) != 0)
                      << '\t' << fds::finalReasonName(record.finalReason) << '\n';
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "flex_decision_sidecar_dump: ERROR: " << error.what() << '\n';
        return 1;
    }
}

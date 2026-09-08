#include "FlexDecisionSidecar.h"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <thread>
#include <vector>

#include <unistd.h>

namespace fds = flex_decision_sidecar;

namespace {

std::vector<unsigned char> readFile(const std::string &path)
{
    std::ifstream input(path.c_str(), std::ios::binary);
    assert(input.good());
    return std::vector<unsigned char>((std::istreambuf_iterator<char>(input)),
                                      std::istreambuf_iterator<char>());
}

fds::WriterConfig config(const std::string &path)
{
    fds::WriterConfig value;
    value.path = path;
    value.starSuiteVersion = "test-version";
    value.sourceRevision = "test-revision";
    value.cachePath = "/test/h01x2_cache.bin";
    return value;
}

void writeComplete(const std::string &path, bool threaded)
{
    fds::Writer writer;
    std::string error;
    assert(writer.open(config(path), error));
    auto writeOne = [&](std::uint64_t ordinal) {
        FlexHashScreenDecision decision;
        decision.action = ordinal == 0 ? FlexHashScreenDecision::Keep
                                       : FlexHashScreenDecision::Pass;
        decision.geneIdx15 = ordinal == 0 ? 17 : 0;
        decision.cacheClass = ordinal == 0 ? FlexHashCacheH0 : 0xFF;
        std::string localError;
        const std::string qname = "@read-" + std::to_string(ordinal) + "/1 extra";
        const bool triageOk = writer.recordTriage(
            ordinal, static_cast<std::uint32_t>(ordinal % 2), ordinal / 2,
            qname.data(), qname.size(), decision, true, true, 1,
            ordinal != 0, false, localError);
        if (!triageOk) std::cerr << "recordTriage: " << localError << '\n';
        assert(triageOk);
        if (ordinal != 0) {
            assert(writer.recordAlignment(
                ordinal, true, ordinal % 2 != 0,
                static_cast<std::uint16_t>(20 + ordinal),
                ordinal % 2 != 0 ? fds::kReasonAlignmentGenomic
                                 : fds::kReasonAlignmentProbe,
                localError));
        }
    };
    if (threaded) {
        std::vector<std::thread> threads;
        for (std::uint64_t ordinal = 0; ordinal < 8; ++ordinal) {
            threads.push_back(std::thread(writeOne, 7 - ordinal));
        }
        for (std::thread &thread : threads) thread.join();
    } else {
        for (std::uint64_t ordinal = 0; ordinal < 8; ++ordinal) writeOne(ordinal);
    }
    assert(writer.finalize(8, error));
}

} // namespace

int main()
{
    // Published FNV-1a-64 vector for the single byte "a".
    assert(fds::normalizedReadNameHash("a", 1) == UINT64_C(0xaf63dc4c8601ec8c));
    assert(fds::normalizedReadNameHash("@read/1 rest", 12)
           == fds::normalizedReadNameHash("read/2", 6));

    char directoryTemplate[] = "/tmp/star-flex-decision-sidecar-XXXXXX";
    const char *directory = ::mkdtemp(directoryTemplate);
    assert(directory != nullptr);
    const std::string serial = std::string(directory) + "/serial.bin";
    const std::string threaded = std::string(directory) + "/threaded.bin";
    writeComplete(serial, false);
    writeComplete(threaded, true);
    assert(readFile(serial) == readFile(threaded));

    fds::Reader reader;
    std::string error;
    assert(reader.open(serial, error));
    assert(reader.header().complete);
    assert(reader.header().schemaVersion == fds::kSchemaVersion);
    assert(reader.header().recordBytes == fds::kRecordBytes);
    assert(reader.header().totalReads == 8);
    assert(reader.header().recordsWritten == 8);
    assert(reader.header().cachePath == "/test/h01x2_cache.bin");
    assert(reader.validateAll(error));
    for (std::uint64_t ordinal = 0; ordinal < 8; ++ordinal) {
        fds::Record record;
        assert(reader.read(ordinal, record, error));
        assert(record.statusFlags & fds::kRecordPresent);
        assert(record.statusFlags & fds::kNameHashPresent);
        if (ordinal == 0) {
            assert(record.cacheClass == FlexHashCacheH0);
            assert(record.finalReason == fds::kReasonCacheKeep);
            assert(record.geneIdx15 == 17);
        } else {
            assert(record.statusFlags & fds::kAlignmentRan);
            assert(record.statusFlags & fds::kAlignmentResolved);
            assert(record.geneIdx15 == 0);
            assert(record.alignmentGeneIdx15 == 20 + ordinal);
        }
    }

    const std::string anchors = std::string(directory) + "/anchors.bin";
    {
        fds::Writer writer;
        assert(writer.open(config(anchors), error));

        FlexHashScreenDecision agreed;
        agreed.action = FlexHashScreenDecision::Pass;
        agreed.residualAnchorGeneIdx15 = 42;
        assert(writer.recordTriage(0, 0, 0, "anchor-agree", 12, agreed,
                                   true, true, 1, true, false, error));
        assert(writer.recordAlignment(0, true, false, 42,
                                      fds::kReasonAlignmentProbe, error));

        FlexHashScreenDecision disagreed;
        disagreed.action = FlexHashScreenDecision::Pass;
        disagreed.residualAnchorGeneIdx15 = 43;
        assert(writer.recordTriage(1, 0, 1, "anchor-disagree", 15, disagreed,
                                   true, true, 1, true, false, error));
        assert(writer.recordAlignment(1, false, false, 99,
                                      fds::kReasonAlignmentAnchorDisagree,
                                      error));

        FlexHashScreenDecision absent;
        absent.action = FlexHashScreenDecision::Deny;
        absent.negativeCode = FlexHashNegHalfNoAnchor;
        assert(writer.recordTriage(2, 0, 2, "anchor-absent", 13, absent,
                                   true, true, 1, false, false, error));

        FlexHashScreenDecision ambiguous;
        ambiguous.action = FlexHashScreenDecision::Deny;
        ambiguous.negativeCode = FlexHashNegHalfGeneAmbig;
        assert(writer.recordTriage(3, 0, 3, "anchor-ambiguous", 16, ambiguous,
                                   true, true, 1, false, false, error));

        FlexHashScreenDecision scoreFail;
        scoreFail.action = FlexHashScreenDecision::Deny;
        scoreFail.negativeCode = FlexHashNegHalfScoreFail;
        assert(writer.recordTriage(4, 0, 4, "probe-score-fail", 16, scoreFail,
                                   true, true, 1, false, false, error));

        FlexHashScreenDecision splitProbe;
        splitProbe.action = FlexHashScreenDecision::Deny;
        splitProbe.negativeCode = FlexHashNegHalfSplitProbe;
        assert(writer.recordTriage(5, 0, 5, "probe-split", 11, splitProbe,
                                   true, true, 1, false, false, error));

        assert(writer.finalize(6, error));
        fds::Reader audit;
        assert(audit.open(anchors, error));
        assert(audit.validateAll(error));
        fds::Record record;
        assert(audit.read(0, record, error));
        assert(record.statusFlags & fds::kResidualAnchorUnique);
        assert(record.statusFlags & fds::kAlignmentAnchorAgreed);
        assert(record.residualAnchorGeneIdx15 == 42);
        assert(record.alignmentGeneIdx15 == 42);
        assert(audit.read(1, record, error));
        assert(record.statusFlags & fds::kAlignmentAnchorDisagreed);
        assert(record.residualAnchorGeneIdx15 == 43);
        assert(record.alignmentGeneIdx15 == 99);
        assert(record.finalReason == fds::kReasonAlignmentAnchorDisagree);
        assert(audit.read(2, record, error));
        assert(record.statusFlags & fds::kResidualAnchorAbsent);
        assert(!(record.statusFlags & fds::kCacheTerminal));
        assert(record.cacheClass == 0xFF);
        assert(record.finalReason == fds::kReasonResidualNoAnchor);
        assert(audit.read(3, record, error));
        assert(record.statusFlags & fds::kResidualAnchorAmbiguous);
        assert(!(record.statusFlags & fds::kCacheTerminal));
        assert(record.cacheClass == 0xFF);
        assert(record.finalReason == fds::kReasonResidualAnchorAmbiguous);
        assert(audit.read(4, record, error));
        assert(record.statusFlags & fds::kProbeScoreFailed);
        assert(record.finalReason == fds::kReasonProbeScoreFail);
        assert(audit.read(5, record, error));
        assert(record.statusFlags & fds::kProbeSplit);
        assert(record.finalReason == fds::kReasonProbeSplit);
    }

    const std::string singleN = std::string(directory) + "/single-n.bin";
    {
        fds::Writer writer;
        assert(writer.open(config(singleN), error));
        FlexHashScreenDecision decision;
        decision.action = FlexHashScreenDecision::Keep;
        decision.geneIdx15 = 42;
        decision.cacheClass = FlexHashCacheH1;
        decision.singleN = true;
        decision.singleNCacheClass = FlexHashCacheH1X2;
        assert(writer.recordTriage(0, 3, 9, "read-n", 6, decision,
                                   true, true, 2, false, false, error));
        assert(writer.finalize(1, error));
        fds::Reader one;
        assert(one.open(singleN, error));
        fds::Record record;
        assert(one.read(0, record, error));
        assert(record.statusFlags & fds::kSingleNAttempted);
        assert(record.statusFlags & fds::kSingleNResolved);
        assert(record.cacheClass == FlexHashCacheH1);
        assert(record.matchedCacheClass == FlexHashCacheH1X2);
    }

    const std::string missing = std::string(directory) + "/missing.bin";
    {
        fds::Writer writer;
        assert(writer.open(config(missing), error));
        FlexHashScreenDecision decision;
        decision.action = FlexHashScreenDecision::Deny;
        assert(writer.recordTriage(1, 0, 1, "read", 4, decision,
                                   false, true, 0xFF, false, false, error));
        assert(!writer.finalize(2, error));
    }

    const std::string outOfRange = std::string(directory) + "/out-of-range.bin";
    {
        fds::Writer writer;
        assert(writer.open(config(outOfRange), error));
        FlexHashScreenDecision decision;
        decision.action = FlexHashScreenDecision::Deny;
        assert(writer.recordTriage(0, 0, 0, "read-0", 6, decision,
                                   false, true, 0xFF, false, false, error));
        assert(writer.recordTriage(2, 0, 2, "read-2", 6, decision,
                                   false, true, 0xFF, false, false, error));
        assert(!writer.finalize(2, error));
    }

    for (const std::string &path : {serial, threaded, singleN, anchors, missing,
                                    missing + ".tmp", outOfRange,
                                    outOfRange + ".tmp"}) {
        std::remove(path.c_str());
    }
    assert(::rmdir(directory) == 0);
    std::cout << "FlexDecisionSidecar tests passed\n";
    return 0;
}

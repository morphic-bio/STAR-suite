#include "SoloReadBarcode.h"
#include "SplitSeq.h"

// Barcode-read handling for the two assays whose cell barcode is assembled
// from separately matched segments. getCBandUMI() hands the read to one of
// these and returns; none of the ordinary whitelist paths run.

void SoloReadBarcode::getCBandUMIpipseq(char **readSeq, uint64 *readLen)
{
    bSeq.assign(readSeq[pSolo.barcodeRead], readLen[pSolo.barcodeRead]);
    pipDecoded = pSolo.pipseqDecoder->decode(bSeq);
    ++pipStats[0];
    cbQual.clear(); umiQual.clear(); cbMatchInd.clear(); cbMatchQual.clear();
    cbMatchString.clear(); urSeq.clear(); urValid = false; urPacked = UINT32_MAX;
    cbMatch = -1; umiCheck = -23; umiB = 0;
    // BAM tags: CR and UR as read at the chosen phase, CB after correction;
    // "-" where undefined (no phase, tied phases, rejected barcode).
    const bool phased = pipDecoded.status != pipseq::DecodeStatus::InvalidBarcode &&
                        pipDecoded.status != pipseq::DecodeStatus::PhaseTie;
    cbSeq = pipseq::Decoder::rawCell(bSeq, phased ? pipDecoded.phase : -1);
    umiSeq = pipseq::Decoder::rawBin(bSeq, phased ? pipDecoded.phase : -1);
    cbSeqCorrected = "-";
    if (pipDecoded.status == pipseq::DecodeStatus::Accepted) {
        cbSeqCorrected = pipseq::Decoder::cellSequence(pipDecoded.cell);
        cbMatch = pipDecoded.mismatches == 0 ? 0 : 1;
        ++pipStats[cbMatch == 0 ? 1 : 2];
        if (pipDecoded.bin >= 0) {
            umiB = pipDecoded.bin;
            umiCheck = 0; // AAA/CCC/GGG/TTT are valid BI values.
        } else ++pipStats[6];
    } else if (pipDecoded.status == pipseq::DecodeStatus::InvalidBarcode) ++pipStats[3];
    else if (pipDecoded.status == pipseq::DecodeStatus::NoMatch) ++pipStats[4];
    else if (pipDecoded.status == pipseq::DecodeStatus::PhaseTie) ++pipStats[5];
    else ++pipStats[7];
}

void SoloReadBarcode::getCBandUMIsplitseq(char **readSeq, char **readQual, uint64 *readLen,
                                          const char *readName, uint64 readNameLength, uint32 fileIndex, uint64 oneBasedReadNumber)
{
    ++splitStats[0];
    bSeq.assign(readSeq[pSolo.barcodeRead], readLen[pSolo.barcodeRead]);
    bQual.assign(readQual[pSolo.barcodeRead], readLen[pSolo.barcodeRead]);
    const auto decoded = pSolo.splitseqDecoder->decode(bSeq);
    cbSeq = decoded.rawCell;
    cbSeqCorrected = "-";
    umiSeq = decoded.umi;
    cbQual.clear(); umiQual.clear();
    cbMatchInd.clear(); cbMatchQual.clear(); cbMatchString.clear();
    urSeq.clear(); urPacked = UINT32_MAX; urValid = false;
    umiB = 0; umiCheck = -23; cbMatch = -11;
    const bool tied = pSolo.splitseqTieWriter && decoded.status == splitseq::Status::Ambiguous;
    if (pSolo.splitseqTieWriter) {
        const bool trainPrefix = oneBasedReadNumber <= splitseq::frequencyPrefixReads && pSolo.splitseqFrequencyBackup == "yes";
        splitseq::countExact(decoded, splitExactCounts, trainPrefix ? &splitPrefixCounts : nullptr);
    }
    if (decoded.status != splitseq::Status::Accepted && !tied) {
        switch (decoded.status) {
            case splitseq::Status::ShortRead: ++splitStats[3]; cbMatch = -11; break;
            case splitseq::Status::InvalidBarcode: ++splitStats[4]; cbMatch = -2; break;
            case splitseq::Status::NoMatch: ++splitStats[5]; cbMatch = -1; break;
            case splitseq::Status::Ambiguous: ++splitStats[6]; cbMatch = -3; break;
            default: break;
        }
        addStats(cbMatch);
        return;
    }
    const auto& regions = pSolo.splitseqDecoder->regions();
    for (size_t i = 0; i < regions.size(); ++i)
        cbQual.append(bQual, regions[i].offset, regions[i].length);
    const auto umi = pSolo.splitseqDecoder->umiRegion();
    umiQual = bQual.substr(umi.offset, umi.length);
    if (convertCheckUMI()) umiCheck = 0;
    uint64_t key = decoded.cell;
    if (tied) {
        string id = readNameLength == UINT64_MAX ? string(readName) : string(readName, readNameLength);
        if (!id.empty() && id[0] == '@') id.erase(0, 1);
        key = pSolo.splitseqTieWriter->defer(id, fileIndex, bSeq, cbQual, umiQual, umiCheck);
    } else {
        cbSeqCorrected = decoded.canonicalCell;
        for (unsigned char q : cbQual) ++qualHist[q];
        for (unsigned char q : umiQual) ++qualHist[q];
        countSplitseqBarcode(decoded, umiCheck);
    }
    if (umiCheck < 0) { cbMatch = umiCheck; return; }
    cbMatch = decoded.correctedRounds ? 1 : 0;
    cbMatchInd.push_back(key);
    cbMatchString = std::to_string(key);
}

void SoloReadBarcode::countSplitseqBarcode(const splitseq::Decoded& decoded, int umiStatus)
{
    if (decoded.status == splitseq::Status::Ambiguous) {
        ++splitStats[6];
        addStats(-3);
    } else if (umiStatus < 0) {
        ++splitStats[umiStatus == -24 ? 8 : 7];
        addStats(umiStatus);
    } else {
        ++splitStats[decoded.correctedRounds ? 2 : 1];
        for (size_t i = 0; i < pSolo.splitseqDecoder->regions().size(); ++i)
            if (decoded.correctedRounds & (1U << i)) ++splitStats[9 + i];
        addStats(decoded.correctedRounds ? 1 : 0);
    }
}

void SoloReadBarcode::finishSplitseqTies()
{
    pSolo.splitseqTieWriter->finish(*pSolo.splitseqDecoder, splitExactCounts, pSolo.splitseqMinPosterior,
        [this](const splitseq::Decoded& result, int umiStatus, const string& cbQ, const string& umiQ) {
            if (result.status == splitseq::Status::Accepted) {
                for (unsigned char q : cbQ) ++qualHist[q];
                for (unsigned char q : umiQ) ++qualHist[q];
            }
            countSplitseqBarcode(result, umiStatus);
        }, pSolo.splitseqFrequencyBackup == "yes" ? &splitPrefixCounts : nullptr);
    splitseq::ExactCounts().swap(splitExactCounts);
    splitseq::ExactCounts().swap(splitPrefixCounts);
}

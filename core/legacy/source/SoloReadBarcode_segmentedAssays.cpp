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

void SoloReadBarcode::getCBandUMIsplitseq(char **readSeq, char **readQual, uint64 *readLen)
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
    if (decoded.status != splitseq::Status::Accepted) {
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
    for (unsigned char q : cbQual) ++qualHist[q];
    for (unsigned char q : umiQual) ++qualHist[q];
    cbSeqCorrected = decoded.canonicalCell;
    if (!convertCheckUMI()) {
        ++splitStats[umiCheck == -24 ? 8 : 7];
        cbMatch = umiCheck;
        addStats(cbMatch);
        return;
    }
    // Counted only now, so every read falls in exactly one category.
    ++splitStats[decoded.correctedRounds ? 2 : 1];
    for (size_t i = 0; i < regions.size(); ++i)
        if (decoded.correctedRounds & (1U << i)) ++splitStats[9 + i];
    umiCheck = 0;
    cbMatch = decoded.correctedRounds ? 1 : 0;
    cbMatchInd.push_back(decoded.cell); // 64-bit sparse key, remapped after mapping.
    cbMatchString = std::to_string(decoded.cell);
    addStats(cbMatch);
}

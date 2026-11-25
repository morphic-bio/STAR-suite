#include "Genome.h"
#include "Transcript.h"
#include "ErrorWarning.h"
#include "sjAlignSplit.h"
#include <algorithm>
#include <limits>

void Genome::loadAlignExcludedRegions() {
    alignExcludedRegions.clear();
    if (P.alignExcludeRegions == "-") return;
    auto fail = [&](const string &reason) {
        exitWithError("EXITING because of invalid --alignExcludeRegions " + P.alignExcludeRegions + ": " + reason + "\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    };
    if (pGe.gType != 1 || genomeOut.convYes || pGe.transform.type != 0)
        fail("only untransformed genomic indexes are supported");
    ifstream input(P.alignExcludeRegions);
    if (!input) fail("cannot open BED file");
    string line;
    uint64 lineNumber=0, records=0;
    const auto coordinate = [](const string &text, uint64 &out) {
        if (text.empty()) return false;
        out=0;
        for (const char c : text) {
            if (c<'0' || c>'9' || out>(std::numeric_limits<uint64>::max()-(c-'0'))/10) return false;
            out=out*10+(c-'0');
        }
        return true;
    };
    while (getline(input,line)) {
        ++lineNumber;
        istringstream fields(line);
        string chrom, startText, endText;
        if (!(fields>>chrom) || chrom[0]=='#' || chrom=="track" || chrom=="browser") continue;
        uint64 start=0, end=0;
        if (!(fields>>startText>>endText) || !coordinate(startText,start) || !coordinate(endText,end) || start>=end)
            fail("line " + to_string(lineNumber) + ": expected chromosome, start < end, nonnegative integer BED coordinates");
        const auto chr=chrNameIndex.find(chrom);
        if (chr==chrNameIndex.end() || chr->second>=nChrReal)
            fail("line " + to_string(lineNumber) + ": unknown chromosome " + chrom);
        if (end>chrLength[chr->second])
            fail("line " + to_string(lineNumber) + ": interval exceeds chromosome length");
        alignExcludedRegions.push_back({{chrStart[chr->second]+start,chrStart[chr->second]+end}});
        ++records;
    }
    if (input.bad()) fail("error reading BED file");
    if (alignExcludedRegions.empty()) fail("BED contains no intervals; use - to disable exclusion");
    sort(alignExcludedRegions.begin(),alignExcludedRegions.end());
    vector<array<uint64,2>> merged;
    for (const auto &region : alignExcludedRegions) {
        if (!merged.empty() && region[0]<=merged.back()[1]) merged.back()[1]=max(merged.back()[1],region[1]);
        else merged.push_back(region);
    }
    alignExcludedRegions.swap(merged);
    uint64 bases=0;
    for (const auto &region : alignExcludedRegions) bases+=region[1]-region[0];
    P.inOut->logMain << "Runtime alignment exclusion: " << P.alignExcludeRegions << "; " << records
                    << " BED records, " << alignExcludedRegions.size() << " merged intervals, " << bases
                    << " bases. Genome and suffix array unchanged.\n" << flush;
}

bool Genome::alignmentBlockExcluded(uint64 start, uint64 length) const {
    if (length==0 || alignExcludedRegions.empty()) return false;
    // Merged disjoint intervals have monotonically increasing ends.
    const auto next=upper_bound(alignExcludedRegions.begin(),alignExcludedRegions.end(),start,
                               [](uint64 value, const array<uint64,2> &region) { return value<region[1]; });
    return next!=alignExcludedRegions.end() && (*next)[0]<start+length;
}

bool Genome::alignmentSeedExcluded(uint64 packedSA, uint64 length) const {
    uint64 start=packedSA & GstrandMask;
    if (packedSA>>GstrandBit) start=nGenome-(start+length);
    if (start<sjGstart) return alignmentBlockExcluded(start,length);
    uint donor, donorLength, acceptor, acceptorLength, junction;
    if (!sjAlignSplit(start,length,*this,donor,donorLength,acceptor,acceptorLength,junction)) {
        // Non-crossing SJ hits are ignored during stitching, but historically
        // contribute to seed multiplicity. Preserve that behavior outside the BED.
        const uint64 offset=(start-sjGstart)%sjdbLength;
        const uint64 sj=(start-sjGstart)/sjdbLength;
        if (offset<sjdbOverhang && offset+length<=sjdbOverhang)
            return alignmentBlockExcluded(sjDstart[sj]+offset,length);
        if (offset>=sjdbOverhang && offset+length<=2*sjdbOverhang)
            return alignmentBlockExcluded(sjAstart[sj]+offset-sjdbOverhang,length);
        return false;
    }
    return alignmentBlockExcluded(donor,donorLength) || alignmentBlockExcluded(acceptor,acceptorLength);
}

bool Genome::alignmentTranscriptExcluded(const Transcript &tr) const {
    for (uint exon=0; exon<tr.nExons; ++exon)
        if (alignmentBlockExcluded(tr.exons[exon][EX_G],tr.exons[exon][EX_L])) return true;
    return false;
}

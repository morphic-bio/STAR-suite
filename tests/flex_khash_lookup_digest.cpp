// Compare this same driver linked against the previous and new real loader.
#include "ParametersSolo.h"
#include "FlexHashScreen.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

static uint64_t digest = UINT64_C(1469598103934665603);
static void add(uint64_t x) { digest = (digest ^ x) * UINT64_C(1099511628211); }
static void add(const FlexHashScreenDecision& d) {
    add(d.action); add(d.geneIdx15); add(d.cacheClass); add(d.negativeCode);
    add(static_cast<uint8_t>(d.offset)); add(d.probeRegion); add(d.singleN);
    add(d.singleNCacheClass); add(d.probeHammingDistance); add(d.residualAnchorGeneIdx15);
}
int main(int argc, char** argv) {
    if (argc != 4) return 2;
    Parameters p; p.readFilesTypeN = std::stoi(argv[1]);
    ParametersSolo solo; solo.pP = &p; solo.hashScreenFile = argv[2];
    std::string error;
    auto& cache = FlexHashScreenCache::instance();
    if (!cache.ensureLoaded(solo, &error)) { std::cerr << error << '\n'; return 1; }
    add(cache.recordCount()); add(cache.h0RecordCount()); add(cache.h1DenyRecordCount());
    std::ifstream in(argv[3]); std::string s; uint64_t n = 0;
    while (in >> s) {
        ++n;
        add(cache.classifyReadH0Offset0(s.data(), s.size()));
        add(cache.classifyReadH0H1Offset0(s.data(), s.size()));
        add(cache.classifyReadH0H1Offset0SingleN(s.data(), s.size()));
        add(cache.classifyReadH1X2SeedExtend(s.data(), s.size()));
        for (unsigned sample = 0; sample < 4; ++sample) {
            add(cache.classifyRead(s.data(), s.size(), sample));
            add(cache.classifyReadH0Only(s.data(), s.size(), sample));
        }
        if (s.size() < 50) continue;
        uint64_t lo = 0, hi = 0, nmask = 0;
        for (unsigned i = 0; i < 50; ++i) {
            auto b = std::string("ACGT").find(s[i]);
            if (b == std::string::npos) { nmask |= UINT64_C(1) << i; b = 0; }
            if (i < 32) lo |= uint64_t(b) << (2*i);
            else hi |= uint64_t(b) << (2*(i-32));
        }
        if (!nmask) {
            add(cache.classifyCbqH0Offset0(lo, hi));
            add(cache.classifyCbqH0H1Offset0(lo, hi));
        }
        add(cache.classifyCbqH0H1Offset0SingleN(lo, hi, nmask));
        add(cache.classifyCbqH1X2SeedExtend(lo, hi, nmask));
    }
    if (!n) return 1;
    std::cout << "queries=" << n << " digest=" << digest << '\n';
}

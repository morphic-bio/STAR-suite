// Exercise the actual private implementation without adding a production API.
// The numerical oracle below freezes the pre-optimization recurrence. Distance
// uses an independent scalar Levenshtein implementation, including N as a byte.
#include "../source/SpatialR1Decoder.cpp"
#include <random>

namespace reference {
double errorProbability(int phred) {
    const double probability =
        std::pow(10.0, -static_cast<double>(std::max(0, phred)) / 10.0);
    return std::max(1.0e-10, std::min(0.75, probability));
}
double pairLikelihood(char observed, char candidate, int phred) {
    const char obs = static_cast<char>(std::toupper(static_cast<unsigned char>(observed)));
    const char truth = static_cast<char>(std::toupper(static_cast<unsigned char>(candidate)));
    if ((obs != 'A' && obs != 'C' && obs != 'G' && obs != 'T') ||
        (truth != 'A' && truth != 'C' && truth != 'G' && truth != 'T')) {
        return std::log(0.25);
    }
    const double error = errorProbability(phred);
    return obs == truth ? std::log1p(-error) : std::log(error / 3.0);
}
double likelihood(const std::string& observed, const std::string& qualities,
                  const std::string& candidate, int gapQ, int missingQ) {
    std::vector<double> previous(candidate.size() + 1, 0.0);
    std::vector<double> current(candidate.size() + 1, 0.0);
    const double deletion = std::log(errorProbability(gapQ));
    for (size_t j = 1; j <= candidate.size(); ++j) previous[j] = previous[j - 1] + deletion;
    for (size_t i = 1; i <= observed.size(); ++i) {
        const int q = i - 1 < qualities.size()
            ? std::max(0, static_cast<int>(static_cast<unsigned char>(qualities[i - 1])) - 33)
            : missingQ;
        const double insertion = std::log(errorProbability(q) / 4.0);
        current[0] = previous[0] + insertion;
        for (size_t j = 1; j <= candidate.size(); ++j) {
            current[j] = std::max({
                previous[j - 1] + pairLikelihood(observed[i - 1], candidate[j - 1], q),
                previous[j] + insertion,
                current[j - 1] + deletion,
            });
        }
        previous.swap(current);
    }
    return previous[candidate.size()];
}
int distance(const std::string& a, const std::string& b) {
    std::vector<int> previous(b.size() + 1), current(b.size() + 1);
    for (size_t j=0; j<=b.size(); ++j) previous[j]=static_cast<int>(j);
    for (size_t i=1; i<=a.size(); ++i) {
        current[0]=static_cast<int>(i);
        for (size_t j=1; j<=b.size(); ++j)
            current[j]=std::min({previous[j]+1,current[j-1]+1,previous[j-1]+(a[i-1]!=b[j-1])});
        previous.swap(current);
    }
    return previous.back();
}
}

void require(bool condition, const char* message, size_t trial) {
    if (!condition) {
        std::cerr << message << " at case " << trial << '\n';
        std::exit(1);
    }
}
std::string sequence(std::mt19937_64& rng, size_t length) {
    const std::string alphabet="ACGTNacgtnX";
    std::string result;
    for (size_t i=0; i<length; ++i) result+=alphabet[rng()%alphabet.size()];
    return result;
}
void checkLikelihood(const std::string& observed, const std::string& qualities,
                     const std::string& candidate, int gapQ, int missingQ, size_t trial) {
    const double expected=reference::likelihood(observed,qualities,candidate,gapQ,missingQ);
    const double actual=raw_candidate_alignment_log_likelihood(observed,qualities,candidate,gapQ,missingQ);
    require(std::memcmp(&expected,&actual,sizeof(double))==0,"bitwise likelihood mismatch",trial);
}
int main() {
    std::mt19937_64 rng(20261026);
    size_t likelihoodCases=0, distanceCases=0, scanCases=0;
    // Every possible quality byte, with canonical, lowercase and ambiguous bases.
    for (int quality=0; quality<256; ++quality) {
        for (int base=0; base<256; ++base) {
            checkLikelihood(std::string(1,static_cast<char>(base)),
                            std::string(1,static_cast<char>(quality)),"ACGTNacgtn",30,30,likelihoodCases++);
            checkLikelihood("ACGTNacgtn",std::string(10,static_cast<char>(quality)),
                            std::string(1,static_cast<char>(base)),30,30,likelihoodCases++);
        }
    }
    // Empty strings, missing/truncated qualities, insertions, deletions, ties,
    // clamp extremes and longer dynamic-program rows all retain exact bits.
    for (size_t trial=0; trial<20000; ++trial) {
        const size_t n=rng()%81, m=rng()%81;
        std::string observed=sequence(rng,n), candidate=sequence(rng,m);
        if (trial%3==0) candidate=observed;
        if (trial%5==0 && !candidate.empty()) candidate.erase(rng()%candidate.size(),1);
        if (trial%7==0) candidate.insert(rng()%(candidate.size()+1),1,'N');
        std::string qualities(rng()%(n+5),'!');
        for (char& q:qualities) q=static_cast<char>(rng()%256);
        const int qValues[]={-20,0,1,10,30,60,100,10000};
        checkLikelihood(observed,qualities,candidate,qValues[rng()%8],qValues[rng()%8],likelihoodCases++);
    }
    for (size_t m=0; m<=64; ++m) {
        for (size_t n=0; n<=64; ++n) {
            for (int repeat=0; repeat<3; ++repeat) {
                const std::string a=sequence(rng,m), b=sequence(rng,n);
                const int expected=reference::distance(a,b);
                const Myers64Pattern prepared(a.data(),static_cast<int>(a.size()));
                require(prepared.distance(b.data(),static_cast<int>(b.size()))==expected,"prepared distance mismatch",distanceCases);
                require(edit_distance_myers64(b.data(),static_cast<int>(b.size()),a.data(),static_cast<int>(a.size()))==expected,"distance symmetry mismatch",distanceCases++);
            }
        }
    }
    for (size_t trial=0; trial<400; ++trial) {
        const std::string seq=sequence(rng,rng()%90);
        const int start=static_cast<int>(rng()%100)-5;
        const int observedLength=static_cast<int>(rng()%82);
        const int maxDistance=static_cast<int>(rng()%4);
        std::vector<std::string> targets;
        for (size_t t=0;t<40;++t) targets.push_back(sequence(rng,rng()%82));
        if (start>=0 && observedLength>0 && start+observedLength<=static_cast<int>(seq.size())) {
            const std::string exact=seq.substr(start,observedLength);
            targets.push_back(exact);
            targets.push_back(exact.substr(1));
            targets.push_back(exact+"N");
            targets.push_back(std::string(1,'X')+exact.substr(1));
        }
        std::vector<std::pair<int,int>> expected;
        if (start>=0 && observedLength>0 && start+observedLength<=static_cast<int>(seq.size())) {
            for (size_t t=0;t<targets.size();++t) {
                if (std::abs(observedLength-static_cast<int>(targets[t].size()))>maxDistance) continue;
                const int d=reference::distance(seq.substr(start,std::min(observedLength,64)),targets[t].substr(0,64));
                if (d<=maxDistance) expected.emplace_back(static_cast<int>(t),d);
            }
        }
        std::map<std::pair<int,int>,std::vector<PackedHit>> cache;
        for (int repeat=0;repeat<2;++repeat) {
            const HitSpan actual=tiered_target_scan_span(seq,start,observedLength,targets,maxDistance,cache);
            require(actual.count==expected.size(),"scan candidate count mismatch",scanCases);
            for (size_t i=0;i<expected.size();++i) {
                require(actual.data[i].idx==expected[i].first && actual.data[i].distance==expected[i].second,
                        "scan candidate order/distance mismatch",scanCases);
            }
            ++scanCases;
        }
    }
    std::cout << "PASS: " << likelihoodCases << " bitwise likelihood cases, "
              << distanceCases << " independent distance cases, " << scanCases
              << " target-scan/cache cases\n";
}

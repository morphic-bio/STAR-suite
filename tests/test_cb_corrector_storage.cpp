#include "solo/CbCorrector.h"
#include "solo/LegacyCbCorrector.h"
#include <cassert>
#include <iostream>
#include <map>
#include <random>
#include <thread>
#include <atomic>

static uint32_t packNative(const std::string& s) {
    uint32_t key=0;
    for(size_t i=0;i<s.size();++i) key |= uint32_t(std::string("ACGT").find(s[i])) << (2*i);
    return key;
}
static void check(const CbCorrector& now,const LegacyCbCorrector& old,const std::string& query) {
    const auto a=now.correct(query);const auto b=old.correct(query);
    if(a.whitelistIdx!=b.whitelistIdx || a.hammingDist!=b.hammingDist ||
       a.ambiguous!=b.ambiguous || a.ambiguousIdx!=b.ambiguousIdx) {
        std::cerr<<"Mismatch for "<<query<<std::endl;std::abort();
    }
    if(query.find_first_not_of("ACGT")==std::string::npos && query.size()<=16) {
        uint32_t ai=100,bi=200;uint8_t ah=99,bh=98;
        const auto x=now.correctPackedCbq(packNative(query),ai,ah);
        const auto y=old.correctPackedCbq(packNative(query),bi,bh);
        assert(x==y && ai==bi && ah==bh);
    }
}
static void run(const std::vector<std::string>& whitelist,const std::vector<std::string>& queries) {
    for(bool native:{false,true}) for(int distance:{0,1}) {
        CbCorrector now(whitelist,distance,native);
        LegacyCbCorrector old(whitelist,distance,native);
        std::map<uint32_t,std::vector<uint32_t>> candidates;
        now.forEachAmbiguousVariant([&](uint32_t key,CbCorrector::CandidateView view) {
            assert(view.size()>1);
            candidates[key]=std::vector<uint32_t>(view.begin(),view.end());
        });
        const auto& expected=old.getAmbiguousVariants();
        assert(candidates.size()==expected.size());
        for(const auto& kv:expected) assert(candidates.at(kv.first)==kv.second);
        for(const auto& query:queries) check(now,old,query);
        for(size_t i=0;i<whitelist.size();++i) {
            check(now,old,whitelist[i]);
            if(whitelist[i].find_first_not_of("ACGT")==std::string::npos)
                assert(now.decodePackedKey(native?packNative(whitelist[i]):packNative(std::string(whitelist[i].rbegin(),whitelist[i].rend())),whitelist[i].size())==whitelist[i]);
        }
        // These tables are shared by all mapping workers after construction.
        std::vector<std::thread> threads;
        for(size_t worker=0;worker<4;++worker)
            threads.emplace_back([&,worker]{for(size_t i=worker;i<queries.size();i+=4)check(now,old,queries[i]);});
        for(auto& thread:threads)thread.join();
    }
}
int main() {
    std::vector<std::string> exhaustive;
    for(unsigned bits=0;bits<15625;++bits) {unsigned n=bits;std::string s(6,'A');for(auto& c:s){c="ACGTN"[n%5];n/=5;}exhaustive.push_back(s);}
    std::mt19937 rng(1943);
    std::vector<std::string> wl;
    for(unsigned i=0;i<90;++i){std::string s(6,'A');for(auto& c:s)c="ACGT"[rng()%4];wl.push_back(s);}
    wl.push_back(wl[4]); // Preserve duplicate candidate multiplicity and exact last-index semantics.
    exhaustive.push_back("acgtac");exhaustive.push_back("ACGT?C");
    run(wl,exhaustive);
    run({}, {"", "AAAAAA", "N", "NNNNNN"});
    std::vector<std::string> dense;
    const std::string center(16,'A');
    for(size_t pos=0;pos<16;++pos)for(char c:std::string("CGT")){auto s=center;s[pos]=c;dense.push_back(s);}
    std::vector<std::string> queries=dense;queries.push_back(center);queries.push_back("NNNNNAAAAAAAAAAA");
    for(unsigned i=0;i<20000;++i){auto s=dense[rng()%dense.size()];for(unsigned j=0;j<i%4;++j)s[rng()%16]="ACGTN"[rng()%5];queries.push_back(s);}
    run(dense,queries); // 48 candidates must survive storage even though correction rejects >5.
    for(unsigned n:{1,2,5,6}){
        std::vector<std::string> subset(dense.begin(),dense.begin()+n);run(subset,queries);
        CbCorrector c(subset);const auto m=c.correct(center);
        assert(n==1 ? (m.whitelistIdx==1 && m.hammingDist==1) : n<=5 ? (m.ambiguous && m.ambiguousIdx.size()==n) : (!m.ambiguous && m.whitelistIdx==0));
    }
    dense.insert(dense.begin(),center);run(dense,queries); // Exact match wins over ambiguous H1.
    std::cout<<"PASS: legacy parity for exact/H1/N expansion, ordered candidates, >5 limit, duplicate inputs, both packing orders and shared readers\n";
}

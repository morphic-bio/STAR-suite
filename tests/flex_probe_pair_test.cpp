#include "FlexProbePairKhash.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
static void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
static FlexProbeRecord record(const std::string& s,unsigned gene) {
    FlexProbeRecord r {};
    for(char b:s){r.key.hi=(r.key.hi<<2)|(r.key.lo>>62);r.key.lo=(r.key.lo<<2)|std::string("ACGT").find(b);}
    r.value.geneAndRegion=gene|(1u<<30);r.value.sample=1;return r;
}
int main(int argc,char** argv) {
    require(argc==2,"usage: pair_test fixture.bin");
    const std::string a(25,'A'),c(25,'C'),g(25,'G'),t(25,'T');
    std::string fourth=a+c;fourth[0]='C';fourth[25]='G';
    std::vector<FlexProbeRecord> rows={record(a+c,1),record(a+g,2),record(t+c,3),record(fourth,1)};
    auto marker=record(g+t,4);marker.value.cacheClass=4;marker.value.sample=0;rows.push_back(marker);
    std::sort(rows.begin(),rows.end(),[](const FlexProbeRecord& a,const FlexProbeRecord& b){
        return a.key.hi<b.key.hi || (a.key.hi==b.key.hi && a.key.lo<b.key.lo);
    });
    struct Header {char magic[8];uint16_t version,k;uint32_t size;uint64_t count;};
    Header h{{'F','H','0','1','S','E','Q','1'},3,50,24,rows.size()};
    {std::ofstream f(argv[1],std::ios::binary);f.write(reinterpret_cast<char*>(&h),sizeof(h));f.write(reinterpret_cast<char*>(rows.data()),rows.size()*sizeof(rows[0]));}
    FlexHashCacheStorage source;std::string error;require(source.open(argv[1],&error),error.c_str());
    FlexProbePairKhash pair;pair.build(source);
    // Both halves individually have multiple parent candidates. Intersection
    // resolves the single common full probe, including shared-gene probes.
    std::string q=a+c;q[10]='C';q[35]='A';auto r=pair.classifyRead(q.data());
    require(r.route==FlexProbePairKhash::DoubleMatch && r.decision.geneIdx15==1 && r.decision.action==FlexHashScreenDecision::Keep,"ambiguous halves must intersect");
    q=a+g;q[5]='T';q[6]='T';r=pair.classifyRead(q.data());
    require(r.route==FlexProbePairKhash::HalfMatch && r.decision.geneIdx15==2 && r.decision.probeHammingDistance==2,"half match must score full probe");
    q=a+g;for(unsigned i=0;i<11;++i)q[i]='T';r=pair.classifyRead(q.data());
    require(r.route==FlexProbePairKhash::HalfMatch && r.decision.action==FlexHashScreenDecision::Deny && r.decision.negativeCode==FlexHashNegHalfScoreFail,"Hamming floor must remain ten");
    q=t+g;r=pair.classifyRead(q.data());
    require(r.route==FlexProbePairKhash::DifferentProbes && r.decision.action==FlexHashScreenDecision::Deny,"different probes must not be joined");
    q=a+c;q[0]='G';q[25]='A';r=pair.classifyRead(q.data());
    require(r.route==FlexProbePairKhash::PairAmbiguous && r.decision.action==FlexHashScreenDecision::Deny,"two common probes sharing a gene must be denied");
    q=a+g;q[10]='N';r=pair.classifyRead(q.data());
    require(r.decision.action==FlexHashScreenDecision::Keep && r.decision.geneIdx15==2 && r.decision.singleN,"single N must preserve full-window merge policy");
    const std::string snapshot=std::string(argv[1])+".half.khash";
    require(pair.writeSnapshot(snapshot,&error),error.c_str());
    FlexProbePairKhash mapped;require(mapped.open(snapshot,&error),error.c_str());
    require(pair.sameTables(mapped),"stored half tables must be byte-identical");
    require(!pair.writeSnapshot(snapshot,&error),"must not overwrite stored half cache");
    require(!mapped.verifySource(source,&error),"partial source must not be accepted by converter");
    // Exercise every N position and both APIs after mmap, including two Ns.
    for(unsigned i=0;i<50;++i)for(unsigned count=0;count<3;++count) {
        q=a+g;q[i]=count?'N':"ACGT"[i%4];if(count==2)q[(i+7)%50]='N';
        auto x=pair.classifyRead(q.data()),y=mapped.classifyRead(q.data());
        require(x.decision.action==y.decision.action && x.decision.geneIdx15==y.decision.geneIdx15 &&
                x.decision.cacheClass==y.decision.cacheClass && x.decision.negativeCode==y.decision.negativeCode &&
                x.decision.singleN==y.decision.singleN && x.decision.probeHammingDistance==y.decision.probeHammingDistance,
                "mapped N and half-match decision parity");
    }
    FlexProbeRecord found;auto raw=record(a+g,2);
    require(mapped.find(raw.key,1,true,found) && found.value.geneAndRegion==raw.value.geneAndRegion,"sample H0 lookup");
    require(!mapped.find(raw.key,2,false,found),"sample mismatch must not become sample-agnostic H1X2");
    const std::string corrupt=snapshot+".corrupt";
    {std::ifstream in(snapshot,std::ios::binary);std::ofstream out(corrupt,std::ios::binary);out<<in.rdbuf();}
    {std::fstream f(corrupt,std::ios::binary|std::ios::in|std::ios::out);f.seekp(8);uint32_t bad=99;f.write(reinterpret_cast<char*>(&bad),4);}
    FlexProbePairKhash invalid;require(!invalid.open(corrupt,&error),"unknown half format must fail closed");
    std::cout<<"PASS: half-list intersection, direct half-match scoring, distance limit, split probes, same-gene ambiguity and single N\n";
}

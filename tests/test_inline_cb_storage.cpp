#include "InlineCBCorrection.h"
#include <algorithm>
#include <iostream>
#include <map>
#include <random>
ParametersSolo::~ParametersSolo() {}
int main() {
    ParametersSolo solo;
    solo.cbL=16;solo.cbWLyes=true;
    const std::string center(16,'A');
    for(unsigned p=0;p<16;++p)for(char b:std::string("CGT")){
        auto s=center;s[p]=b;solo.cbWLstr.push_back(s);
    }
    std::mt19937 random(5544);
    for(unsigned i=0;i<128;++i){auto s=center;for(char& c:s)c="ACGT"[random()%4];solo.cbWLstr.push_back(s);}
    solo.cbWLstr.push_back(solo.cbWLstr.back());
    InlineCBCorrection::initializeWhitelist(solo);
    std::cout<<"S "<<InlineCBCorrection::exactMapSize()<<' '<<InlineCBCorrection::variantMapSize()<<' '
             <<InlineCBCorrection::variantCollisionSize()<<' '<<InlineCBCorrection::variantCollisionMaxFanout()<<'\n';
    std::vector<std::string> queries=solo.cbWLstr;queries.push_back(center);
    for(const auto& s:solo.cbWLstr)for(unsigned p=0;p<16;++p)for(char c:std::string("ACGTN")){
        auto q=s;q[p]=c;queries.push_back(q);
    }
    for(unsigned n=0;n<20000;++n){auto s=center;for(char& c:s)c="ACGT"[random()%4];queries.push_back(s);}
    unsigned ordinal=0;
    for(const auto& query:queries){
        std::string fast,nseq;std::vector<std::string> closest;
        int f=InlineCBCorrection::fastPathCorrection(query,fast);
        int c=InlineCBCorrection::findClosestBarcodes(query,closest);
        int n=InlineCBCorrection::checkSequenceAndCorrectForN(query,1,nseq);
        uint64_t key=InlineCBCorrection::packCBForLookup(query);
        uint32_t idx=0;
        int packed=key==UINT64_MAX?-1:InlineCBCorrection::fastPathCorrectionPacked(key,idx);
        std::cout<<ordinal++<<' '<<f<<' '<<fast<<' '<<packed<<' '<<idx<<' '<<c;
        for(const auto& s:closest)std::cout<<' '<<s;
        std::cout<<" N "<<n<<' '<<nseq<<'\n';
        if(InlineCBCorrection::isAmbiguousVariant(key)){
            std::string qual(16,'I');qual[ordinal%16]='!';
            InlineCBCorrection::recordAmbiguousCB(key,query,qual);
            std::cout<<"R "<<key<<' '<<InlineCBCorrection::resolveAmbiguousVariant(key,query,qual,solo)<<'\n';
        }
    }
    std::unordered_map<uint64_t,InlineCBCorrection::MergedAmbigEntry> merged;
    InlineCBCorrection::mergeAmbiguousShards(merged);
    std::map<uint64_t,InlineCBCorrection::MergedAmbigEntry> ordered(merged.begin(),merged.end());
    for(const auto& row:ordered){
        std::cout<<"A "<<row.first<<' '<<row.second.count<<' '<<row.second.evidenceReads<<' '
                 <<row.second.cbSeq<<' '<<row.second.cbQual;
        for(auto p:row.second.parents)std::cout<<' '<<p;
        std::cout<<'\n';
    }
    InlineCBCorrection::clearEvidence();InlineCBCorrection::clearAmbiguous();InlineCBCorrection::clearWhitelist();
}

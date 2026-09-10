#include "ParametersSolo.h"
#include "FlexHashScreen.h"
#ifdef FLEX_PAIR_BENCH
#include "FlexProbePairKhash.h"
#endif
#include <chrono>
#include <fstream>
#include <iostream>
#include <cstring>
#include <sys/resource.h>
#include <stdexcept>

struct Query { uint64_t lo, hi, nmask; };
struct Verdict { uint8_t bytes[16]; };
using Clock = std::chrono::steady_clock;
static double seconds(Clock::time_point t) { return std::chrono::duration<double>(Clock::now()-t).count(); }
static uint64_t rss() { rusage r {};getrusage(RUSAGE_SELF,&r);return r.ru_maxrss; }
static Verdict verdict(const FlexHashScreenDecision& d, unsigned route = 0) {
    Verdict v {}; v.bytes[0]=d.action;v.bytes[1]=d.geneIdx15&255;v.bytes[2]=d.geneIdx15>>8;
    v.bytes[3]=d.cacheClass;v.bytes[4]=d.negativeCode;v.bytes[5]=d.probeRegion;
    v.bytes[6]=d.probeHammingDistance;v.bytes[7]=d.singleN;v.bytes[8]=d.singleNCacheClass;
    v.bytes[9]=static_cast<uint8_t>(d.offset);v.bytes[10]=route;return v;
}
static FlexHashScreenDecision classify(FlexHashScreenCache& cache, const Query& q, unsigned mode) {
    FlexHashScreenDecision d;
    if (mode==20) {
        d = q.nmask ? cache.classifyCbqH0H1Offset0SingleN(q.lo,q.hi,q.nmask)
                    : cache.classifyCbqH0H1Offset0(q.lo,q.hi);
        if (d.action==FlexHashScreenDecision::Pass && cache.h1x2ProbeIndexReady())
            d=cache.classifyCbqH1X2SeedExtend(q.lo,q.hi,q.nmask);
    } else {
        char s[50];for(unsigned i=0;i<50;++i)s[i]=(q.nmask&(UINT64_C(1)<<i))?'N':"ACGT"[((i<32?q.lo>>(2*i):q.hi>>(2*(i-32)))&3)];
        d=q.nmask?cache.classifyReadH0H1Offset0SingleN(s,50):cache.classifyReadH0H1Offset0(s,50);
        if(d.action==FlexHashScreenDecision::Pass && cache.h1x2ProbeIndexReady())d=cache.classifyReadH1X2SeedExtend(s,50);
    }
    return d;
}
int main(int argc,char** argv) {
    if(argc!=6){std::cerr<<"usage: flex_cache_bench large|pair|audit cache queries mode verdicts\n";return 2;}
    const std::string engine=argv[1];const unsigned mode=std::stoul(argv[4]);
    const auto start=Clock::now();std::string error;
    FlexHashScreenCache* cache=nullptr;
#ifdef FLEX_PAIR_BENCH
    FlexProbePairKhash pair;
    if(engine=="pair" || engine=="audit") {
        FlexHashCacheStorage source;
        if(!source.open(argv[2],&error)){std::cerr<<error<<'\n';return 1;}
        pair.build(source);
        std::cout<<"{\"half_keys_left\":"<<pair.keys(0)<<",\"half_keys_right\":"<<pair.keys(1)
                 <<",\"half_resident_bytes\":"<<pair.bytes()<<",\"parents\":"<<pair.probeCount()<<"}"<<std::endl;
    }
#endif
    if(engine!="pair") {
        Parameters p;p.readFilesTypeN=mode;ParametersSolo solo;solo.pP=&p;solo.hashScreenFile=argv[2];
        cache=&FlexHashScreenCache::instance();
        if(!cache->ensureLoaded(solo,&error)){std::cerr<<error<<'\n';return 1;}
        if(!cache->hasH1X2() || !cache->h1x2ProbeIndexReady() || !cache->h1x2ProbeCount()) {
            std::cerr<<"benchmark requires a complete H0/H1X2 cache and matching loader headers\n";return 1;
        }
        std::cout<<"{\"records\":"<<cache->recordCount()<<",\"probe_parents\":"<<cache->h1x2ProbeCount()<<"}"<<std::endl;
    }
    std::cout<<"{\"phase\":\"load\",\"seconds\":"<<seconds(start)<<",\"max_rss_kib\":"<<rss()<<"}"<<std::endl;
#ifdef FLEX_PAIR_BENCH
    if(engine=="audit") {
        FlexHashCacheStorage source;if(!source.open(argv[2],&error))return 1;
        uint64_t actionGene=0,metadata=0;uint64_t transitions[5][5] {};uint64_t routes[6] {};
        const auto t=Clock::now();
        for(uint64_t i=0;i<source.recordCount();++i) {
            const auto& r=source.record(i);auto k=FlexHashCacheStorage::cbqKey(r.key);
            const auto a=verdict(classify(*cache,Query{k.lo,k.hi,0},20));
            const auto p=pair.classify(k);const auto b=verdict(p.decision);
            ++routes[p.route];++transitions[a.bytes[0]][b.bytes[0]];
            if(memcmp(a.bytes,b.bytes,3)) {
                ++actionGene;
                if(actionGene<=20)std::cout<<"{\"difference_record\":"<<i<<",\"old_action\":"<<unsigned(a.bytes[0])<<",\"new_action\":"<<unsigned(b.bytes[0])
                 <<",\"old_gene\":"<<(unsigned(a.bytes[1])+(unsigned(a.bytes[2])<<8))<<",\"new_gene\":"<<(unsigned(b.bytes[1])+(unsigned(b.bytes[2])<<8))<<"}"<<std::endl;
            }
            if(memcmp(a.bytes,b.bytes,10))++metadata;
        }
        std::cout<<"{\"audit_records\":"<<source.recordCount()<<",\"action_gene_differences\":"<<actionGene<<",\"any_decision_differences\":"<<metadata
          <<",\"seconds\":"<<seconds(t)<<",\"transitions\":[";
        for(unsigned a=0;a<5;++a)for(unsigned b=0;b<5;++b)std::cout<<(a||b?",":"")<<transitions[a][b];
        std::cout<<"]}"<<std::endl;return 0;
    }
#endif
    std::ifstream in(argv[3],std::ios::binary);std::ofstream out(argv[5],std::ios::binary);
    if(!in||!out)return 1;
    std::vector<Query> queries(65536);std::vector<Verdict> result(queries.size());
    uint64_t count=0,actions[5]={},routes[6]={};double classifySeconds=0;
    const auto begin=Clock::now();
    while(in.read(reinterpret_cast<char*>(queries.data()),queries.size()*sizeof(Query)) || in.gcount()) {
        if(in.gcount()%sizeof(Query))return 1;
        const size_t n=in.gcount()/sizeof(Query);const auto t=Clock::now();
        for(size_t i=0;i<n;++i) {
#ifdef FLEX_PAIR_BENCH
            if(engine=="pair") {
                const auto& q=queries[i];FlexProbePairKhash::Outcome p;
                if(mode==20)p=pair.classify(FlexProbeKey{q.lo,q.hi},q.nmask);
                else {
                    char s[50];for(unsigned j=0;j<50;++j)s[j]=(q.nmask&(UINT64_C(1)<<j))?'N':"ACGT"[((j<32?q.lo>>(2*j):q.hi>>(2*(j-32)))&3)];
                    p=pair.classifyRead(s);
                }
                result[i]=verdict(p.decision,p.route);++routes[p.route];
            } else
#endif
            result[i]=verdict(classify(*cache,queries[i],mode));
            ++actions[result[i].bytes[0]];
        }
        classifySeconds+=seconds(t);count+=n;
        out.write(reinterpret_cast<char*>(result.data()),n*sizeof(Verdict));if(!out)return 1;
    }
    out.close();if(!out)return 1;
    std::cout<<"{\"queries\":"<<count<<",\"classify_seconds\":"<<classifySeconds<<",\"query_wall_seconds\":"<<seconds(begin)
      <<",\"max_rss_kib\":"<<rss()<<",\"actions\":[";
    for(unsigned i=0;i<5;++i)std::cout<<(i?",":"")<<actions[i];
    std::cout<<"],\"routes\":[";
    for(unsigned i=0;i<6;++i)std::cout<<(i?",":"")<<routes[i];
    std::cout<<"]}"<<std::endl;
}

#include "pf_api.h"
#include "common.h"
#include <atomic>
#include <cassert>
#include <cstring>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <zlib.h>
struct PermitCounts { std::atomic<int> used{0}, acquired{0}, released{0}; };
static uint64_t acquire(void* p) { auto& c=*static_cast<PermitCounts*>(p); ++c.used; ++c.acquired; return 0; }
static void release(void* p,uint64_t,uint64_t,uint64_t,uint64_t) { auto& c=*static_cast<PermitCounts*>(p); --c.used; ++c.released; }
static std::string sequence(const std::string& path) {
    gzFile f=gzopen(path.c_str(),"rb"); assert(f); char line[65536]; assert(gzgets(f,line,sizeof(line))); assert(gzgets(f,line,sizeof(line))); gzclose(f);
    std::string s(line); while(!s.empty()&&(s.back()=='\n'||s.back()=='\r'))s.pop_back(); return s;
}
int main(int argc,char** argv) {
    assert(argc==3);const std::string input=argv[1],out=argv[2];
    const auto barcode=sequence(input+"/R1.gz"), feature=sequence(input+"/R2.gz");
    for(int kind=0;kind<4;++kind) {
        PermitCounts permits;auto* config=pf_config_create();pf_config_set_skip_emptydrops(config,1);pf_config_set_skip_qc_outputs(config,1);
        pf_config_set_search_threads(config,1);pf_config_set_feature_offset(config,0);pf_config_set_permit_hooks(config,acquire,release,&permits);
        auto* context=pf_init(config);pf_config_destroy(config);assert(context);
        assert(pf_load_whitelist(context,(input+"/whitelist.txt").c_str())==PF_OK);
        assert(pf_load_feature_ref(context,(input+"/features.csv").c_str())==PF_OK);
        pf_direct_range_job* job=nullptr;assert(pf_direct_range_begin(context,(out+"/case"+std::to_string(kind)).c_str(),"sample",2,2,&job)==PF_OK);
        pf_read_record_view valid{};valid.barcode_sequence={barcode.data(),barcode.size()};valid.feature_sequence={feature.data(),feature.size()};
        auto invalid=valid;
        if(kind==0) invalid.barcode_quality={"I",1};
        if(kind==1) invalid.feature_sequence2={"ACGT",4};
        if(kind==2) invalid.feature_sequence.length=SIZE_MAX;
        // A maximum-sized LF-terminated sequence with default qualities is
        // valid. Its synthesized quality must exclude the sequence newline.
        std::string boundary=feature+std::string(LINE_LENGTH-2-feature.size(),'A')+'\n';
        if(kind==3) invalid.feature_sequence={boundary.data(),boundary.size()};
        const pf_read_record_view records[]={valid,invalid};std::vector<std::thread> workers;
        for(int i=0;i<2;++i)workers.emplace_back([&,i]{assert(pf_direct_range_process_record_views(job,i,records,2)==(kind==3?PF_OK:PF_ERR_INVALID_ARG));});
        for(auto& worker:workers)worker.join();
        assert(permits.used==0); // failed batches must release before job cleanup
        assert(pf_direct_range_end(job,nullptr)==(kind==3?PF_OK:PF_ERR_INVALID_ARG));
        assert(permits.used==0&&permits.acquired==permits.released);if(kind!=3)assert(std::strlen(pf_get_error(context))>0);
        pf_destroy(context);std::cout<<"PASS direct batch boundary case "<<kind<<" permits_returned="<<permits.released<<"\n";
    }
}

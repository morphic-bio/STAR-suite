#include "vb_engine.h"
#include <cassert>
#include <cmath>
#include <iostream>
static bool close(double a,double b,double eps=1e-6){return std::abs(a-b)<=eps||(std::isnan(a)&&std::isnan(b));}
int main(){
    for(unsigned kind=0;kind<4;++kind){
        ECTable table;TranscriptState state;state.resize(kind==0?0:64);table.n_transcripts=state.n;
        for(size_t i=0;i<state.n;++i){state.lengths[i]=1000+i;state.eff_lengths[i]=900+i;state.names[i]="t"+std::to_string(i);}
        if(kind>1)for(unsigned e=0;e<500;++e){EC ec;unsigned x=kind==2?e%12:(e%8)*8;ec.transcript_ids={x,x+1};ec.weights={0.3+(e%2)*0.1,0.7-(e%2)*0.1};ec.count=1+e%17;table.ecs.push_back(ec);}
        table.n_ecs=table.ecs.size();EMParams p;p.use_vb=true;p.threads=8;unsigned updatesA=0,updatesB=0;
        p.effective_length_update=[&](uint32_t,TranscriptState& s,const std::vector<double>&){++updatesA;for(double& x:s.eff_lengths)x-=5;return true;};
        auto aState=state;auto a=run_vb(table,aState,p);p.component_partition=true;
        p.effective_length_update=[&](uint32_t,TranscriptState& s,const std::vector<double>&){++updatesB;for(double& x:s.eff_lengths)x-=5;return true;};
        auto bState=state;auto b=run_vb(table,bState,p);assert(a.iterations==b.iterations&&a.converged==b.converged&&updatesA==updatesB);
        for(size_t i=0;i<state.n;++i){assert(close(a.counts[i],b.counts[i]));assert(close(a.tpm[i],b.tpm[i]));assert(aState.eff_lengths[i]==bState.eff_lengths[i]);}
        assert(close(a.final_ll,b.final_ll));std::cout<<"PASS component case="<<kind<<" iterations="<<a.iterations<<"\n";
    }
}

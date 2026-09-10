#include "scrna_api.h"
#include "OrdMagStage.h"
#include "OrdMagRank.h"
#ifndef CALLER_BASELINE
#include "ScrnaTrace.h"
#include "BoundedSampleTasks.h"
#endif
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <cassert>
#include <sys/resource.h>
int main(int argc, char** argv) {
    if (argc != 6) return 2; // matrix directory, output, MC workers, layout, simulations
    const std::string dir=argv[1]; std::ifstream m(dir+"/matrix.mtx");
    std::string line; uint32_t rows=0,cols=0; size_t nnz=0;
    while(std::getline(m,line)) if(!line.empty()&&line[0]!='%') {std::istringstream s(line);s>>rows>>cols>>nnz;break;}
    if (!cols) return 3;
    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> cells(cols);
    uint32_t g,c,v; while(m>>g>>c>>v) cells.at(c-1).push_back({g-1,v});
    std::vector<uint32_t> genes,counts,offsets(cols),entries(cols),umis(cols),strided,wordOffsets(cols);
    std::vector<std::string> names;std::ifstream b(dir+"/barcodes.tsv");while(std::getline(b,line))names.push_back(line);
    assert(names.size()==cols);std::vector<char*> ptrs;for(auto& n:names)ptrs.push_back(&n[0]);
    for(uint32_t i=0;i<cols;++i) {offsets[i]=genes.size(); wordOffsets[i]=strided.size();entries[i]=cells[i].size();for(auto x:cells[i]){genes.push_back(x.first);counts.push_back(x.second);umis[i]+=x.second;strided.insert(strided.end(),{x.first,999,x.second,888});}}
    cells.clear();cells.shrink_to_fit();
    scrna_matrix_input input{};input.n_cells=cols;input.n_features=rows;input.umi_counts=umis.data();input.barcodes=ptrs.data();input.sparse_gene_ids=genes.data();input.sparse_counts=counts.data();input.sparse_cell_index=offsets.data();input.n_genes_per_cell=entries.data();input.sparse_nnz=genes.size();
    auto* cfg=scrna_ed_config_create();cfg->use_bootstrap=1;cfg->n_expected_cells=0;cfg->umi_min=100;cfg->sim_n=std::stoul(argv[5]);cfg->mc_threads=std::stoul(argv[3]);cfg->use_fdr_gate=cfg->apply_bh_correction=1;
    scrna_ed_result result{};auto start=std::chrono::steady_clock::now();int rc;
#ifdef CALLER_BASELINE
    rc=scrna_emptydrops_run_with_rank_options(&input,cfg,nullptr,4,&result);
#else
    ScrnaTrace trace;SparseCountView view; view.genes=strided.data();view.counts=strided.data()+2;view.geneWords=strided.size();view.countWords=strided.size()-2;view.stride=4;view.offsets=wordOffsets.data();view.entries=entries.data();view.cells=cols;
    rc=scrnaEmptyDropsTrace(&input,cfg,nullptr,4,&result,&trace,std::string(argv[4])=="strided"?&view:nullptr);
#endif
    const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();if(rc){std::cerr<<result.error_message;return 4;}
    std::ofstream out(argv[2]);out<<std::setprecision(17)<<result.n_simple_cells<<'\t'<<result.n_tail_cells<<'\t'<<result.n_ed_passers<<'\t'<<result.retain_threshold<<'\t'<<result.min_umi<<'\n';
    for(size_t i=0;i<result.n_barcodes;++i)out<<"cell\t"<<result.barcodes[i]<<'\n';
    for(size_t i=0;i<result.n_candidates;++i){auto& x=result.candidates[i];out<<x.barcode<<'\t'<<x.cell_index<<'\t'<<x.umi_count<<'\t'<<x.p_value<<'\t'<<x.p_adjusted<<'\t'<<x.obs_log_prob<<'\t'<<x.passes_raw_p<<'\t'<<x.passes_fdr<<'\t'<<x.is_simple_cell<<'\n';}
    struct rusage usage;getrusage(RUSAGE_SELF,&usage);std::cout<<"caller_seconds="<<wall<<" peak_rss_kib="<<usage.ru_maxrss<<" input_nnz="<<nnz<<"\n";
#ifndef CALLER_BASELINE
    // Validate wide-offset rejection before dereferencing, and shared scheduling cleanup.
    uint64_t big=UINT64_MAX;uint32_t one=2;SparseCountView bad=view;bad.cells=1;bad.offsets64=&big;bad.entries=&one;bool rejected=false;try{bad.validate(1);}catch(const std::exception&){rejected=true;}assert(rejected);
    std::vector<uint64_t> weights{10,20,0,100,5,15};std::vector<int> visits(6);auto stats=scrna::boundedSampleTasks(weights,3,35,[&](size_t i){++visits[i];});assert(stats.peakWorkers<=3&&stats.permitsReturned==3&&stats.peakEstimatedBytes==100);for(int n:visits)assert(n==1);
    bool threw=false;try{scrna::boundedSampleTasks(weights,2,35,[&](size_t i){if(i==1)throw std::runtime_error("injected");});}catch(const std::exception&){threw=true;}assert(threw);
#endif
    scrna_ed_result_free(&result);scrna_ed_config_destroy(cfg);
}

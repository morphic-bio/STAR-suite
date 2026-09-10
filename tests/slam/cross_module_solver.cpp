#include "SlamSolver.h"
#include "libem/slam_vb_overdisp.h"
#include <iostream>
#include <iomanip>
int main() {
    std::cout << std::hexfloat;
    for (unsigned gene=0;gene<1500;++gene) {
        MismatchHistogram h;
        for (unsigned k=0;k<gene%61;++k) h[slamPackMismatchKey(25+k,k%9)] = (gene*17+k*31)%71+0.5;
        const double pe=gene%3?0.001:0.002,pc=gene%2?0.05:0.03;
        auto a=SlamSolver(pe,pc).solve(h);auto b=SlamVbOverdispSolver(pe,pc,50,1.1,1.2).solve(h);
        std::cout<<gene<<'\t'<<a.ntr<<'\t'<<a.sigma<<'\t'<<a.log_likelihood<<'\t'<<a.converged<<'\t'<<b.ntr_map<<'\t'<<b.ntr_mean<<'\t'<<b.log_likelihood<<'\t'<<b.iters<<'\t'<<b.converged<<'\n';
    }
    for(double p:{0.,1.,-0.1,0.5})for(unsigned k:{0u,5u,6u}){
        MismatchHistogram h{{slamPackMismatchKey(5,k),4}};
        auto a=SlamSolver(p,p).solve(h);auto b=SlamVbOverdispSolver(p,p,0,0,0).solve(h);
        std::cout<<p<<'\t'<<a.ntr<<'\t'<<a.log_likelihood<<'\t'<<a.converged<<'\t'<<b.ntr_mean<<'\t'<<b.ntr_map<<'\t'<<b.log_likelihood<<'\t'<<b.iters<<'\t'<<b.converged<<'\n';
    }
}

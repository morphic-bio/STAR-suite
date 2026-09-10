#include "SlamQuant.h"
#include <cassert>
#include <iostream>
#include <cmath>
static bool equal(double a,double b){return a==b||(std::isnan(a)&&std::isnan(b));}
int main() {
    SlamQuant serial(129), parallel(129);parallel.setFitWorkers(5);
    for(unsigned g=0;g<128;++g)for(unsigned k=0;k<(g%17)+1;++k){serial.addRead(g,30+k,k%5,1+g%9);parallel.addRead(g,30+k,k%5,1+g%9);}
    for(bool vb:{false,true}) {
        SlamFitParameters p{0.001,0.05,50,1,1,vb};
        auto a=serial.fits(p);auto b=parallel.fits(p);
        for(size_t i=0;i<a.size();++i){assert(equal(a[i].mean,b[i].mean)&&equal(a[i].map,b[i].map)&&equal(a[i].sigma,b[i].sigma)&&equal(a[i].likelihood,b[i].likelihood)&&a[i].converged==b[i].converged);}
        auto passes=parallel.fitPasses();parallel.fits(p);assert(parallel.fitPasses()==passes);
        p.error=0.002;parallel.fits(p);assert(parallel.fitPasses()==++passes);
        parallel.addRead(128,35,2,4);parallel.fits(p);assert(parallel.fitPasses()==++passes);
        SlamQuant other(129);other.addRead(128,35,1,8);parallel.merge(other);parallel.fits(p);assert(parallel.fitPasses()==++passes);
        // Keep the next model's biological input identical after invalidation checks.
        serial.addRead(128,35,2,4);serial.merge(other);
    }
    std::cout<<"PASS: exact serial/parallel fits, reuse, parameter/addRead/merge invalidation\n";
}

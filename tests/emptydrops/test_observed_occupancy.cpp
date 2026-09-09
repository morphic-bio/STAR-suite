#include "ObservedTagOccupancy.h"
#include <cassert>
#include <iostream>
#include <algorithm>
std::string dna(unsigned n, unsigned len) { std::string s(len,'A'); for (auto& c:s) { c="ACGT"[n%4]; n/=4; } return s; }
int main() {
  assert(fitObservedTagOccupancy({}).rejectedGems.empty());
  assert(fitObservedTagOccupancy({dna(0,24),dna(0,24)+"-1"}).rejectedGems.empty());
  std::vector<std::string> calls;
  for(unsigned i=0;i<1000;++i) for(unsigned j=0;j<2;++j) calls.push_back(dna(i,16)+dna(j,8));
  for(unsigned j=0;j<8;++j) calls.push_back(dna(1000,16)+dna(j,8));
  auto r=fitObservedTagOccupancy(calls);
  assert(r.cutoff==7 && r.rejectedGems.size()==1 && r.rejectedGems.count(dna(1000,16)));
  std::reverse(calls.begin(),calls.end()); calls.push_back(calls.front()+"-1");
  auto reversed=fitObservedTagOccupancy(calls);
  assert(r.lambda==reversed.lambda && r.rejectedGems==reversed.rejectedGems);
  bool failed=false; try { fitObservedTagOccupancy({"ACGT"}); } catch(const std::runtime_error&) {failed=true;}
  assert(failed);
  std::cout << "Observed occupancy tests PASS: lambda="<<r.lambda<<" cutoff="<<r.cutoff<<"\n";
}

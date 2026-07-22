#ifndef CODE_signalFromBAM
#define CODE_signalFromBAM
#if defined(WITH_CHROMAP) && WITH_CHROMAP
#include <htslib/sam.h>
#else
#include "htslib/htslib/sam.h"
#endif
#include  <fstream>
#include <string>
#include "Stats.h"
#include "Parameters.h"

using namespace std;

void signalFromBAM(const string bamFileName, const string sigFileName, Parameters P);

#endif

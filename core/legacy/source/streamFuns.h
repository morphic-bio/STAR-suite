#ifndef CODE_streamFuns
#define CODE_streamFuns

#include "Parameters.h"
#include <fstream>

unsigned long long fstreamReadBig(std::ifstream &S, char* A, unsigned long long N);

// Parallel equivalent of fstreamReadBig for one flat, contiguous destination
// buffer. The genome index files are a byte-for-byte image of the array they
// load into, so file byte k always lands at A[k] and any byte range can be
// filled independently. S is used only for its current offset, which becomes
// the base of the parallel reads, and is advanced to match on return. Falls
// back to the serial path when the offset or the file cannot be obtained.
unsigned long long fstreamReadBigParallel(const std::string &path, std::ifstream &S,
                                          char* A, unsigned long long N, int nThreads);

void fstreamWriteBig(std::ofstream &S, char* A, unsigned long long N, std::string fileName, std::string errorID, Parameters &P) ;

fstream  &fstrOpen  (std::string fileName, std::string errorID, Parameters &P, bool flagDelete);
fstream  &fstrOpenBinary  (std::string fileName, std::string errorID, Parameters &P, bool flagDelete);
ofstream &ofstrOpen (std::string fileName, std::string errorID, Parameters &P);
ifstream &ifstrOpen (std::string fileName, std::string errorID, std::string solutionString, Parameters &P);
ifstream &ifstrOpenGenomeFile (std::string fileName, std::string errorID, Parameters &P);
void createDirectory(const string dirPathIn, const mode_t dirPerm, const string dirParameter, Parameters &P);

void copyFile(string fileIn, string fileOut);
#endif

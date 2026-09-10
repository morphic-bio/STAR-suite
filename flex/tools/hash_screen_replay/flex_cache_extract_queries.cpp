#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <zlib.h>
struct Query { uint64_t lo,hi,nmask; };
int main(int argc,char** argv) {
    if(argc<4){std::cerr<<"usage: extract output.bin reads_per_lane R2.fastq.gz ...\n";return 2;}
    std::ofstream out(argv[1],std::ios::binary);if(!out)return 1;
    const uint64_t limit=std::stoull(argv[2]);std::vector<Query> buffer;buffer.reserve(65536);
    uint8_t lut[256];for(auto& b:lut)b=4;lut['A']=0;lut['C']=1;lut['G']=2;lut['T']=3;
    for(int lane=3;lane<argc;++lane) {
        gzFile in=gzopen(argv[lane],"rb");if(!in)return 1;gzbuffer(in,1<<20);
        char name[1024],s[1024],plus[1024],qual[1024];uint64_t count=0;
        for(;count<limit;++count) {
            if(!gzgets(in,name,sizeof(name)) || !gzgets(in,s,sizeof(s)) || !gzgets(in,plus,sizeof(plus)) || !gzgets(in,qual,sizeof(qual)))return 1;
            Query q{};
            for(unsigned i=0;i<50;++i){if(s[i]=='\n'||s[i]=='\0')return 1;uint64_t b=lut[static_cast<unsigned char>(s[i])];if(b>3){q.nmask|=UINT64_C(1)<<i;b=0;}if(i<32)q.lo|=b<<(2*i);else q.hi|=b<<(2*(i-32));}
            buffer.push_back(q);
            if(buffer.size()==buffer.capacity()){out.write(reinterpret_cast<const char*>(buffer.data()),buffer.size()*sizeof(Query));buffer.clear();}
        }
        gzclose(in);std::cerr<<argv[lane]<<": "<<count<<" windows\n";
    }
    out.write(reinterpret_cast<const char*>(buffer.data()),buffer.size()*sizeof(Query));return out?0:1;
}

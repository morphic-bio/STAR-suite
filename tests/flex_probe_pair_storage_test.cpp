#include "FlexProbePairKhash.h"
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>

static void require(bool ok,const char* message) { if(!ok)throw std::runtime_error(message); }
static FlexProbeRecord record(const std::string& s,uint8_t cls,uint16_t sample) {
    FlexProbeRecord r {};
    for(char c:s) { r.key.hi=(r.key.hi<<2)|(r.key.lo>>62);r.key.lo=(r.key.lo<<2)|std::string("ACGT").find(c); }
    r.value={1u|(2u<<30),cls,0,sample};return r;
}
static void writeRecords(const std::string& path,std::vector<FlexProbeRecord> rows) {
    std::sort(rows.begin(),rows.end(),[](const FlexProbeRecord& a,const FlexProbeRecord& b) {
        if(a.key.hi!=b.key.hi)return a.key.hi<b.key.hi;
        if(a.key.lo!=b.key.lo)return a.key.lo<b.key.lo;
        return a.value.sample<b.value.sample;
    });
    struct Header {char magic[8];uint16_t version,k;uint32_t size;uint64_t count;};
    Header h{{'F','H','0','1','S','E','Q','1'},3,50,24,rows.size()};
    std::ofstream f(path,std::ios::binary);f.write(reinterpret_cast<const char*>(&h),sizeof(h));
    f.write(reinterpret_cast<const char*>(rows.data()),rows.size()*sizeof(rows[0]));require(bool(f),"write fixture");
}
template<class T> static T readAt(const std::string& path,uint64_t offset) {
    T v {};std::ifstream f(path,std::ios::binary);f.seekg(offset);f.read(reinterpret_cast<char*>(&v),sizeof(v));require(bool(f),"read stored fixture");return v;
}
template<class T> static void writeAt(const std::string& path,uint64_t offset,T v) {
    std::fstream f(path,std::ios::binary|std::ios::in|std::ios::out);f.seekp(offset);f.write(reinterpret_cast<const char*>(&v),sizeof(v));require(bool(f),"mutate stored fixture");
}
static void copyFile(const std::string& from,const std::string& to) {
    std::ifstream in(from,std::ios::binary);std::ofstream out(to,std::ios::binary);out<<in.rdbuf();require(bool(out),"copy fixture");
}
int main(int argc,char** argv) {
    require(argc==2,"usage: pair_storage_test fresh-output-directory");
    const std::string dir=argv[1];require(mkdir(dir.c_str(),0755)==0,"output must be fresh");
    std::string probe;for(unsigned i=0;i<50;++i)probe+="ACGT"[i%4];
    std::vector<std::string> halves[2];
    for(unsigned s=0;s<2;++s) {
        const auto exact=probe.substr(s*25,25);halves[s].push_back(exact);
        for(unsigned p=0;p<25;++p)for(char b:std::string("ACGT"))if(b!=exact[p]) {
            auto variant=exact;variant[p]=b;halves[s].push_back(variant);
        }
    }
    std::vector<FlexProbeRecord> rows {record(probe,0,1),record(probe,0,2)};
    for(unsigned l=0;l<76;++l)for(unsigned r=0;r<76;++r)if(l||r)rows.push_back(record(halves[0][l]+halves[1][r],4,0));
    writeRecords(dir+"/complete.bin",rows);
    std::string error;FlexHashCacheStorage source;require(source.open(dir+"/complete.bin",&error),error.c_str());
    FlexProbePairKhash built;built.build(source);require(built.verifySource(source,&error),error.c_str());
    const std::string packed=dir+"/complete.half.khash";require(built.writeSnapshot(packed,&error),error.c_str());
    FlexProbePairKhash mapped;require(mapped.open(packed,&error),error.c_str());require(built.sameTables(mapped),"exact stored arrays");
    require(mapped.verifySource(source,&error),error.c_str());
    FlexProbeRecord got;
    require(mapped.find(rows[1].key,2,true,got) && got.value.sample==2,"preserve second H0 sample");
    require(mapped.find(rows.back().key,2,false,got) && got.value.sample==0 && got.value.cacheClass==4,"preserve global H1X2 fallback");
    rows.pop_back();writeRecords(dir+"/partial.bin",rows);
    FlexHashCacheStorage partial;require(partial.open(dir+"/partial.bin",&error),error.c_str());
    FlexProbePairKhash incomplete;incomplete.build(partial);
    require(!incomplete.verifySource(partial,&error) && error.find("omits a half-H1 variant")!=std::string::npos,
            "reverse verification must reject missing variants even when every present record agrees");
    // Public on-disk v1 offsets: parent-list offset at 80, first half table at 144.
    const uint64_t parentOffset=readAt<uint64_t>(packed,80),flagsOffset=readAt<uint64_t>(packed,144),valuesOffset=readAt<uint64_t>(packed,160);
    const uint32_t buckets=readAt<uint32_t>(packed,168);
    uint32_t occupied=0;
    for(;occupied<buckets;++occupied) {
        const uint32_t flags=readAt<uint32_t>(packed,flagsOffset+(occupied/16)*4);
        if(((flags>>(2*(occupied%16)))&3)==0)break;
    }
    require(occupied<buckets,"fixture needs an occupied half bucket");
    for(unsigned defect=0;defect<4;++defect) {
        const std::string bad=dir+"/bad"+std::to_string(defect)+".khash";copyFile(packed,bad);
        if(defect==0) {struct stat st {};require(stat(bad.c_str(),&st)==0 && truncate(bad.c_str(),st.st_size-1)==0,"truncate fixture");}
        if(defect==1)writeAt<uint64_t>(bad,80,UINT64_MAX);
        if(defect==2)writeAt<uint32_t>(bad,parentOffset,1); // only parent ID 0 exists
        if(defect==3)writeAt<uint64_t>(bad,valuesOffset+uint64_t(occupied)*8,UINT64_MAX);
        FlexProbePairKhash invalid;require(!invalid.open(bad,&error) && !invalid.ready(),"invalid stored half cache must fail closed");
    }
    std::cout<<"PASS: complete source equivalence, missing-variant rejection, sample/global lookup, truncated file, bad sections, parent IDs and spans\n";
}

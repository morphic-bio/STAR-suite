#include "FlexProbePairKhash.h"
#include <cerrno>
#include <fstream>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace {
const char magic[8] = {'F','H','2','5','K','H','0','1'};
struct Table { uint64_t flags, keys, values; uint32_t buckets, entries; };
struct Header {
    char magic[8];
    uint32_t version, headerBytes, endian, encoding, hashAlgorithm, sourceVersion;
    uint64_t fileBytes, records, h0Count, probeCount, probes, rawH0;
    uint64_t parents[2], parentCounts[2];
    Table tables[3];
    uint64_t reserved[6];
};
static_assert(sizeof(Header)==256,"half cache header layout");
bool fail(std::string* error, const std::string& s) { if(error)*error=s;return false; }
uint64_t align(uint64_t x) { return (x+4095)&~UINT64_C(4095); }
bool writeAt(int fd,uint64_t off,const void* p,uint64_t n) {
    const char* c=static_cast<const char*>(p);
    while(n) {
        const ssize_t k=pwrite(fd,c,std::min<uint64_t>(n,64*1024*1024),off);
        if(k<0 && errno==EINTR)continue;
        if(k<=0)return false;
        n-=k;off+=k;c+=k;
    }
    return true;
}
bool keyLess(FlexProbeKey a,FlexProbeKey b) { return a.hi<b.hi || (a.hi==b.hi && a.lo<b.lo); }
}

FlexProbePairKhash::~FlexProbePairKhash() { close(); }
void FlexProbePairKhash::close() {
    if(mapping_) munmap(mapping_,mappingBytes_);
    else {
        for(auto& h:halves_) { free(h.flags);free(h.keys);free(h.vals); }
        free(h0_.flags);free(h0_.keys);free(h0_.vals);
    }
    h0_={};for(auto& h:halves_)h={};mapping_=nullptr;mappingBytes_=0;
    probesData_=rawH0Data_=nullptr;probeCount_=recordCount_=h0Count_=0;
    for(unsigned s=0;s<2;++s) { parentsData_[s]=nullptr;parentCounts_[s]=0;parents_[s].clear(); }
    probes_.clear();rawH0_.clear();
}
bool FlexProbePairKhash::isSnapshot(const std::string& path) {
    std::ifstream f(path,std::ios::binary);char b[8] {};
    return bool(f.read(b,8)) && memcmp(b,magic,8)==0;
}

bool FlexProbePairKhash::open(const std::string& path,std::string* error) {
    close();const uint32_t endian=1;
    if(sizeof(size_t)!=8 || *reinterpret_cast<const uint8_t*>(&endian)!=1)
        return fail(error,"half cache requires a 64-bit little-endian host");
    const int fd=::open(path.c_str(),O_RDONLY);struct stat st {};
    if(fd<0)return fail(error,"cannot open half cache");
    if(fstat(fd,&st) || st.st_size<256) { ::close(fd);return fail(error,"half cache header truncated"); }
    mappingBytes_=st.st_size;mapping_=mmap(nullptr,mappingBytes_,PROT_READ,MAP_PRIVATE,fd,0);::close(fd);
    if(mapping_==MAP_FAILED) { mapping_=nullptr;return fail(error,"cannot mmap half cache"); }
    auto reject=[&](const char* s) {close();return fail(error,s);};
    const auto& h=*static_cast<const Header*>(mapping_);
    if(memcmp(h.magic,magic,8) || h.version!=1 || h.headerBytes!=256 || h.endian!=0x01020304 ||
       h.encoding!=1 || h.hashAlgorithm!=1 || h.sourceVersion!=3 || h.fileBytes!=mappingBytes_ ||
       !h.probeCount || h.probeCount>UINT32_MAX || h.h0Count<h.probeCount || h.records<h.h0Count)
        return reject("invalid or incompatible half cache header");
    for(auto v:h.reserved)if(v)return reject("unsupported half cache fields");
    uint64_t cursor=4096;
    auto section=[&](uint64_t off,uint64_t n,uint64_t width) {
        if(off!=cursor || off>mappingBytes_ || n>(mappingBytes_-off)/width)return false;
        cursor=align(off+n*width);return true;
    };
    if(!section(h.probes,h.probeCount,24) || !section(h.rawH0,h.h0Count,24))return reject("invalid half probe section");
    for(unsigned s=0;s<2;++s) {
        if(h.parentCounts[s]>UINT32_MAX || !section(h.parents[s],h.parentCounts[s],4))return reject("invalid half parent section");
    }
    for(unsigned t=0;t<3;++t) {
        const auto& d=h.tables[t];
        if(d.buckets<4 || d.buckets>(UINT32_C(1)<<31) || (d.buckets&(d.buckets-1)) ||
           d.entries>=uint64_t(d.buckets*__ac_HASH_UPPER+0.5) ||
           !section(d.flags,__ac_fsize(d.buckets),4) || !section(d.keys,d.buckets,t?8:16) ||
           !section(d.values,d.buckets,8))return reject("invalid half hash section");
    }
    if(cursor!=mappingBytes_ || h.tables[0].entries!=h.probeCount)return reject("half cache length/count mismatch");
    char* base=static_cast<char*>(mapping_);
    probesData_=reinterpret_cast<const FlexProbeRecord*>(base+h.probes);
    rawH0Data_=reinterpret_cast<const FlexProbeRecord*>(base+h.rawH0);
    for(unsigned s=0;s<2;++s) {
        parentsData_[s]=reinterpret_cast<const uint32_t*>(base+h.parents[s]);parentCounts_[s]=h.parentCounts[s];
        for(uint64_t i=0;i<parentCounts_[s];++i)if(parentsData_[s][i]>=h.probeCount)return reject("half parent ID out of bounds");
    }
    const auto& d=h.tables[0];
    h0_.flags=reinterpret_cast<khint32_t*>(base+d.flags);h0_.keys=reinterpret_cast<FlexProbeKey*>(base+d.keys);
    h0_.vals=reinterpret_cast<FlexProbeValue*>(base+d.values);h0_.n_buckets=d.buckets;h0_.size=h0_.n_occupied=d.entries;
    h0_.upper_bound=uint64_t(d.buckets*__ac_HASH_UPPER+0.5);
    uint64_t h0Used=0;
    for(khint_t k=0;k<h0_.n_buckets;++k) {
        if(__ac_isdel(h0_.flags,k))return reject("deleted bucket in immutable half H0 hash");
        if(kh_exist(&h0_,k))++h0Used;
    }
    if(h0Used!=d.entries)return reject("half H0 occupancy mismatch");
    for(unsigned s=0;s<2;++s) {
        const auto& d=h.tables[s+1];auto& a=halves_[s];
        a.flags=reinterpret_cast<khint32_t*>(base+d.flags);a.keys=reinterpret_cast<uint64_t*>(base+d.keys);
        a.vals=reinterpret_cast<uint64_t*>(base+d.values);a.n_buckets=d.buckets;a.size=a.n_occupied=d.entries;
        a.upper_bound=uint64_t(d.buckets*__ac_HASH_UPPER+0.5);
        uint64_t used=0;
        for(khint_t k=0;k<a.n_buckets;++k) {
            if(__ac_isdel(a.flags,k))return reject("deleted bucket in immutable half hash");
            if(!kh_exist(&a,k))continue;
            ++used;const uint64_t off=a.vals[k]>>32,n=uint32_t(a.vals[k]);
            if(!n || off>parentCounts_[s] || n>parentCounts_[s]-off || a.keys[k]>>50)return reject("invalid half hash parent span");
            for(uint64_t i=1;i<n;++i)if(parentsData_[s][off+i-1]>=parentsData_[s][off+i])return reject("unordered half parent span");
        }
        if(used!=d.entries)return reject("half hash occupancy mismatch");
    }
    probeCount_=h.probeCount;recordCount_=h.records;h0Count_=h.h0Count;
    return true;
}

bool FlexProbePairKhash::writeSnapshot(const std::string& path,std::string* error) const {
    if(!ready())return fail(error,"half tables are not built");
    Header h {};memcpy(h.magic,magic,8);h.version=1;h.headerBytes=256;h.endian=0x01020304;
    h.encoding=h.hashAlgorithm=1;h.sourceVersion=3;h.records=recordCount_;h.h0Count=h0Count_;h.probeCount=probeCount_;
    uint64_t cursor=4096;
    auto section=[&](uint64_t n,uint64_t width) {uint64_t off=cursor;cursor=align(cursor+n*width);return off;};
    h.probes=section(probeCount_,24);h.rawH0=section(h0Count_,24);
    for(unsigned s=0;s<2;++s) {h.parentCounts[s]=parentCounts_[s];h.parents[s]=section(parentCounts_[s],4);}
    for(unsigned t=0;t<3;++t) {
        auto& d=h.tables[t];d.buckets=t?halves_[t-1].n_buckets:h0_.n_buckets;d.entries=t?halves_[t-1].size:h0_.size;
        d.flags=section(__ac_fsize(d.buckets),4);d.keys=section(d.buckets,t?8:16);d.values=section(d.buckets,8);
    }
    h.fileBytes=cursor;
    const std::string tmp=path+".tmp."+std::to_string(getpid());
    const int fd=::open(tmp.c_str(),O_WRONLY|O_CREAT|O_EXCL,0644);
    if(fd<0)return fail(error,"cannot create half cache temporary file");
    bool ok=ftruncate(fd,h.fileBytes)==0 && writeAt(fd,0,&h,sizeof(h)) &&
        writeAt(fd,h.probes,probesData_,probeCount_*24) && writeAt(fd,h.rawH0,rawH0Data_,h0Count_*24);
    for(unsigned s=0;s<2 && ok;++s)ok=writeAt(fd,h.parents[s],parentsData_[s],parentCounts_[s]*4);
    for(unsigned t=0;t<3 && ok;++t) {
        const auto& d=h.tables[t];
        ok=writeAt(fd,d.flags,t?halves_[t-1].flags:h0_.flags,uint64_t(__ac_fsize(d.buckets))*4) &&
           writeAt(fd,d.keys,t?static_cast<const void*>(halves_[t-1].keys):h0_.keys,uint64_t(d.buckets)*(t?8:16)) &&
           writeAt(fd,d.values,t?static_cast<const void*>(halves_[t-1].vals):h0_.vals,uint64_t(d.buckets)*8);
    }
    if(ok)ok=fsync(fd)==0;
    if(::close(fd)!=0)ok=false;
    // Atomic publish, never overwrite an existing paired artifact.
    if(ok)ok=link(tmp.c_str(),path.c_str())==0;
    const int saved=errno;unlink(tmp.c_str());
    return ok || fail(error,"cannot publish half cache: "+std::string(strerror(saved)));
}

bool FlexProbePairKhash::find(FlexProbeKey rawKey,uint16_t sample,bool h0Only,FlexProbeRecord& out) const {
    uint64_t l=0,r=h0Count_;
    while(l<r) {const uint64_t m=l+(r-l)/2;if(keyLess(rawH0Data_[m].key,rawKey))l=m+1;else r=m;}
    const FlexProbeRecord* fallback=nullptr;bool exact=false;
    for(uint64_t i=l;i<h0Count_ && flexProbeEqual(rawH0Data_[i].key,rawKey);++i) {
        const auto& a=rawH0Data_[i];exact=true;
        if(a.value.sample==sample) {out=a;return true;}
        if(a.value.sample==0)fallback=&a;
    }
    if(fallback) {out=*fallback;return true;}
    if(exact || h0Only)return false;
    const auto d=primary(FlexHashCacheStorage::cbqKey(rawKey)).decision;
    if(d.action==FlexHashScreenDecision::Pass)return false;
    out={rawKey,{uint32_t(d.geneIdx15)|(uint32_t(d.probeRegion)<<30),d.cacheClass,d.negativeCode,0}};return true;
}

bool FlexProbePairKhash::verifySource(const FlexHashCacheStorage& source,std::string* error) const {
    if(source.sourceVersion()!=3 || source.recordCount()!=recordCount_ || source.h0Count()!=h0Count_)
        return fail(error,"half cache/source metadata mismatch");
    // Require the source's complete H0 + H1X2 universe. H1/H2-only caches and
    // sample-specific non-exact records remain on the general cache format.
    for(uint64_t i=0;i<source.recordCount();++i) {
        const auto& r=source.record(i);const auto& v=r.value;
        if((v.cacheClass!=0 && v.cacheClass!=4 && v.cacheClass!=2) ||
           (v.cacheClass!=0 && v.sample) || (v.cacheClass==0 && (v.negativeCode || !(v.geneAndRegion&0x7fff))))
            return fail(error,"source is not a complete H0/H1X2 cache supported by half tables");
        const auto d=primary(FlexHashCacheStorage::cbqKey(r.key)).decision;
        const bool deny=v.cacheClass==2;
        if(d.action!=(deny?FlexHashScreenDecision::Deny:FlexHashScreenDecision::Keep) ||
           d.geneIdx15!=(v.geneAndRegion&0x7fff) || d.cacheClass!=v.cacheClass ||
           d.negativeCode!=v.negativeCode || uint32_t(d.probeRegion)!=(v.geneAndRegion>>30))
            return fail(error,"half primary differs from source record "+std::to_string(i));
    }
    // Check the reverse direction too: missing variants in a partial source
    // must not silently become newly accepted sequences after conversion.
    for(uint64_t id=0;id<probeCount_;++id) {
        uint64_t variants[2][76];
        for(unsigned side=0;side<2;++side) {
            unsigned n=0;const auto k=half(probesData_[id].key,side);variants[side][n++]=k;
            for(unsigned p=0;p<25;++p)for(uint64_t b=0;b<4;++b)if(b!=((k>>(2*p))&3))
                variants[side][n++]=(k&~(UINT64_C(3)<<(2*p)))|(b<<(2*p));
        }
        for(unsigned l=0;l<76;++l)for(unsigned r=0;r<76;++r) {
            const FlexProbeKey key{variants[0][l]|(variants[1][r]<<50),variants[1][r]>>14};
            if(!source.lookup(0,key) && !source.lookup(1,key))
                return fail(error,"source omits a half-H1 variant; keep the general cache format");
        }
    }
    return true;
}

bool FlexProbePairKhash::sameTables(const FlexProbePairKhash& b) const {
    if(!ready() || probeCount_!=b.probeCount_ || recordCount_!=b.recordCount_ || h0Count_!=b.h0Count_ ||
       memcmp(probesData_,b.probesData_,probeCount_*24) || memcmp(rawH0Data_,b.rawH0Data_,h0Count_*24))return false;
    if(h0_.n_buckets!=b.h0_.n_buckets || h0_.size!=b.h0_.size ||
       memcmp(h0_.flags,b.h0_.flags,__ac_fsize(h0_.n_buckets)*4) ||
       memcmp(h0_.keys,b.h0_.keys,uint64_t(h0_.n_buckets)*16) || memcmp(h0_.vals,b.h0_.vals,uint64_t(h0_.n_buckets)*8))return false;
    for(unsigned s=0;s<2;++s) {
        const auto& x=halves_[s];const auto& y=b.halves_[s];
        if(parentCounts_[s]!=b.parentCounts_[s] || x.n_buckets!=y.n_buckets || x.size!=y.size ||
           memcmp(parentsData_[s],b.parentsData_[s],parentCounts_[s]*4) || memcmp(x.flags,y.flags,__ac_fsize(x.n_buckets)*4) ||
           memcmp(x.keys,y.keys,uint64_t(x.n_buckets)*8) || memcmp(x.vals,y.vals,uint64_t(x.n_buckets)*8))return false;
    }
    return true;
}

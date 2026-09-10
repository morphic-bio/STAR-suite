#include "FlexHashCacheStorage.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <limits>
#include <new>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace {
const char recordMagic[8] = {'F','H','0','1','S','E','Q','1'};
const char snapshotMagic[8] = {'F','H','0','1','K','H','0','1'};
struct RecordHeader {
    char magic[8]; uint16_t version, kmerLength;
    uint32_t recordBytes; uint64_t count;
};
struct DiskTable {
    uint64_t flagsOffset, keysOffset, valuesOffset;
    uint32_t buckets, entries;
};
struct SnapshotHeader {
    char magic[8];
    uint32_t version, headerBytes, endian, keyEncoding, hashAlgorithm;
    uint32_t sourceVersion, flags, recordBytes;
    uint64_t fileBytes, recordCount, h0Count, h1Count, recordsOffset, h0Offset;
    DiskTable tables[2];
    uint64_t reserved[13];
};
static_assert(sizeof(RecordHeader) == 24, "legacy header layout");
static_assert(sizeof(SnapshotHeader) == 256, "snapshot header layout");
static_assert(sizeof(khint_t) == 4, "snapshot requires 32-bit khash indices");

bool fail(std::string* error, const std::string& message) {
    if (error) *error = message;
    return false;
}
bool littleEndian() {
    const uint32_t n = 1;
    return *reinterpret_cast<const uint8_t*>(&n) == 1;
}
int tierOf(const FlexProbeValue& v) {
    const bool gene = (v.geneAndRegion & 0x7fff) != 0;
    if (v.cacheClass == 0 && gene) return 0;
    if (!gene || v.cacheClass == 1 || v.cacheClass == 4) return 1;
    return -1; // H2 remains available to the sample-aware full classifier.
}
uint64_t alignPage(uint64_t n) { return (n + 4095) & ~UINT64_C(4095); }
uint64_t flagBytes(uint32_t buckets) { return uint64_t(__ac_fsize(buckets)) * 4; }
bool writeAt(int fd, uint64_t offset, const void* src, uint64_t bytes) {
    const char* p = static_cast<const char*>(src);
    while (bytes) {
        const size_t chunk = static_cast<size_t>(std::min<uint64_t>(bytes, 64 * 1024 * 1024));
        const ssize_t n = pwrite(fd, p, chunk, static_cast<off_t>(offset));
        if (n < 0 && errno == EINTR) continue;
        if (n <= 0) return false;
        p += n; bytes -= n; offset += n;
    }
    return true;
}
uint64_t reverseGroups(uint64_t v) {
    v = ((v & UINT64_C(0x3333333333333333)) << 2) | ((v >> 2) & UINT64_C(0x3333333333333333));
    v = ((v & UINT64_C(0x0f0f0f0f0f0f0f0f)) << 4) | ((v >> 4) & UINT64_C(0x0f0f0f0f0f0f0f0f));
    return __builtin_bswap64(v);
}
}

FlexProbeKey FlexHashCacheStorage::cbqKey(FlexProbeKey key) {
    const uint64_t first18 = reverseGroups(key.hi) >> 28;
    const uint64_t last32 = reverseGroups(key.lo);
    return FlexProbeKey{first18 | ((last32 & UINT64_C(0x0fffffff)) << 36),
                        (last32 >> 28) & UINT64_C(0xfffffffff)};
}

FlexHashCacheStorage::~FlexHashCacheStorage() { close(); }
void FlexHashCacheStorage::close() {
    if (!persisted_) {
        for (auto& h : maps_) { free(h.flags); free(h.keys); free(h.vals); }
    }
    for (auto& h : maps_) h = khash_t(flex_probe){};
    if (mapping_) munmap(mapping_, static_cast<size_t>(mappingBytes_));
    mapping_ = nullptr; mappingBytes_ = 0; records_ = nullptr; h0Indices_ = nullptr;
    ownedH0Indices_.clear(); recordCount_ = h0Count_ = h1Count_ = 0;
    sourceVersion_ = 0; persisted_ = hasH1X2_ = false;
}

FlexProbeValue FlexHashCacheStorage::normalized(FlexProbeValue v) const {
    v.geneAndRegion &= sourceVersion_ >= 3 ? 0xc0007fffu : 0x7fffu;
    if (sourceVersion_ < 2) v.sample = 0;
    return v;
}

bool FlexHashCacheStorage::open(const std::string& path, std::string* error) {
    close();
    if (!littleEndian() || sizeof(size_t) < 8)
        return fail(error, "Flex khash caches require a 64-bit little-endian host");
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) return fail(error, "cannot open cache: " + std::string(strerror(errno)));
    struct stat st {};
    if (fstat(fd, &st) != 0 || st.st_size < 24) {
        ::close(fd); return fail(error, "cannot read cache header");
    }
    mappingBytes_ = static_cast<uint64_t>(st.st_size);
    mapping_ = mmap(nullptr, static_cast<size_t>(mappingBytes_), PROT_READ, MAP_PRIVATE, fd, 0);
    ::close(fd);
    if (mapping_ == MAP_FAILED) {
        mapping_ = nullptr; return fail(error, "cannot mmap cache: " + std::string(strerror(errno)));
    }
    bool ok = false;
    try {
        if (memcmp(mapping_, snapshotMagic, 8) == 0) {
            // Set ownership before any table pointers are rebound, including failures.
            persisted_ = true;
            ok = openSnapshot(error);
        } else if (memcmp(mapping_, recordMagic, 8) == 0) {
            const auto& h = *static_cast<const RecordHeader*>(mapping_);
            if (h.version < 1 || h.version > 3 || h.kmerLength != 50 || h.recordBytes != 24 ||
                h.count > (mappingBytes_ - 24) / 24 || mappingBytes_ != 24 + h.count * 24) {
                fail(error, "cache format or length mismatch");
            } else {
                sourceVersion_ = h.version; recordCount_ = h.count;
                records_ = reinterpret_cast<const FlexProbeRecord*>(static_cast<const char*>(mapping_) + 24);
                ok = build(error);
            }
        } else fail(error, "cache magic mismatch");
    } catch (const std::bad_alloc&) { fail(error, "out of memory building legacy khash cache; use a packed cache"); }
    if (!ok) close();
    return ok;
}

bool FlexHashCacheStorage::build(std::string* error) {
    // Input is already sorted by (hi, lo, sample). Filtering preserves order.
    uint64_t unique[2] = {0, 0}; FlexProbeKey last[2] {}; bool seen[2] = {false, false};
    for (uint64_t i = 0; i < recordCount_; ++i) {
        const auto& r = records_[i];
        if (r.value.cacheClass == 4) hasH1X2_ = true;
        const int t = tierOf(r.value); if (t < 0) continue;
        if (t == 0) ++h0Count_; else ++h1Count_;
        if (!seen[t] || !flexProbeEqual(last[t], r.key)) ++unique[t];
        seen[t] = true; last[t] = r.key;
    }
    ownedH0Indices_.reserve(static_cast<size_t>(h0Count_));
    for (unsigned t = 0; t < 2; ++t) {
        uint64_t buckets = 4;
        while (unique[t] >= uint64_t(buckets * __ac_HASH_UPPER + 0.5)) buckets *= 2;
        if (buckets > (UINT64_C(1) << 31)) return fail(error, "cache exceeds khash capacity");
        if (kh_resize(flex_probe, &maps_[t], static_cast<khint_t>(buckets)) != 0)
            return fail(error, "cannot allocate khash cache");
        // Empty slots are persisted too: initialize them for deterministic files
        // and to avoid writing allocator contents to disk.
        memset(maps_[t].keys, 0, buckets * sizeof(FlexProbeKey));
        memset(maps_[t].vals, 0, buckets * sizeof(FlexProbeValue));
    }
    for (uint64_t i = 0; i < recordCount_; ++i) {
        const auto& r = records_[i]; const int t = tierOf(r.value); if (t < 0) continue;
        if (t == 0) ownedH0Indices_.push_back(i);
        int inserted = 0;
        const auto k = kh_put(flex_probe, &maps_[t], cbqKey(r.key), &inserted);
        if (inserted < 0) return fail(error, "cannot insert into khash cache");
        // Match unordered_map::emplace: first record in sorted sample order wins.
        if (inserted) kh_val(&maps_[t], k) = normalized(r.value);
    }
    h0Indices_ = ownedH0Indices_.data();
    return true;
}

bool FlexHashCacheStorage::openSnapshot(std::string* error) {
    if (mappingBytes_ < sizeof(SnapshotHeader)) return fail(error, "khash header truncated");
    const auto& h = *static_cast<const SnapshotHeader*>(mapping_);
    if (h.version != 1 || h.headerBytes != sizeof(h) || h.endian != 0x01020304 ||
        h.keyEncoding != 1 || h.hashAlgorithm != 1 || h.sourceVersion < 1 || h.sourceVersion > 3 ||
        h.recordBytes != 24 || h.flags > 1 || h.fileBytes != mappingBytes_ ||
        h.recordCount > mappingBytes_ / 24 || h.h0Count > h.recordCount ||
        h.h1Count > h.recordCount - h.h0Count)
        return fail(error, "incompatible or invalid khash header");
    for (uint64_t reserved : h.reserved)
        if (reserved) return fail(error, "unsupported khash header fields");
    const char* base = static_cast<const char*>(mapping_);
    uint64_t next = alignPage(sizeof(h));
    auto section = [&](uint64_t offset, uint64_t count, uint64_t width) {
        if (offset != next || offset > mappingBytes_ || count > (mappingBytes_ - offset) / width) return false;
        next = alignPage(offset + count * width); return true;
    };
    if (!section(h.recordsOffset, h.recordCount, 24) || !section(h.h0Offset, h.h0Count, 8))
        return fail(error, "invalid khash record sections");
    records_ = reinterpret_cast<const FlexProbeRecord*>(base + h.recordsOffset);
    h0Indices_ = reinterpret_cast<const uint64_t*>(base + h.h0Offset);
    for (unsigned t = 0; t < 2; ++t) {
        const auto& d = h.tables[t];
        if (d.buckets < 4 || (d.buckets & (d.buckets - 1)) || d.buckets > (UINT32_C(1) << 31) ||
            d.entries > (t == 0 ? h.h0Count : h.h1Count) ||
            d.entries >= uint64_t(d.buckets * __ac_HASH_UPPER + 0.5) ||
            !section(d.flagsOffset, flagBytes(d.buckets), 1) ||
            !section(d.keysOffset, d.buckets, sizeof(FlexProbeKey)) ||
            !section(d.valuesOffset, d.buckets, sizeof(FlexProbeValue)))
            return fail(error, "invalid khash table sections");
        auto& m = maps_[t];
        m.n_buckets = d.buckets; m.size = m.n_occupied = d.entries;
        m.upper_bound = static_cast<khint_t>(d.buckets * __ac_HASH_UPPER + 0.5);
        m.flags = reinterpret_cast<khint32_t*>(const_cast<char*>(base + d.flagsOffset));
        m.keys = reinterpret_cast<FlexProbeKey*>(const_cast<char*>(base + d.keysOffset));
        m.vals = reinterpret_cast<FlexProbeValue*>(const_cast<char*>(base + d.valuesOffset));
    }
    if (next != mappingBytes_) return fail(error, "invalid khash file length");
    // Only the small H0 index is scanned on startup; never scan/reinsert the
    // hundreds of millions of H1 records or re-sort any records.
    for (uint64_t i = 0; i < h.h0Count; ++i)
        if (h0Indices_[i] >= h.recordCount || (i && h0Indices_[i] <= h0Indices_[i-1]))
            return fail(error, "invalid khash H0 record index");
    sourceVersion_ = static_cast<uint16_t>(h.sourceVersion); recordCount_ = h.recordCount;
    h0Count_ = h.h0Count; h1Count_ = h.h1Count; hasH1X2_ = h.flags & 1;
    return true;
}

bool FlexHashCacheStorage::find(FlexProbeKey key, uint16_t sample, bool h0Only, FlexProbeRecord& out) const {
    const uint64_t count = h0Only ? h0Count_ : recordCount_;
    auto search = [&](uint16_t target) {
        uint64_t lo = 0, hi = count;
        while (lo < hi) {
            const uint64_t mid = lo + (hi - lo) / 2;
            const auto& r = h0Only ? h0Record(mid) : record(mid);
            const uint16_t s = sourceVersion_ >= 2 ? r.value.sample : 0;
            if (r.key.hi < key.hi || (r.key.hi == key.hi &&
                (r.key.lo < key.lo || (r.key.lo == key.lo && s < target)))) lo = mid + 1;
            else hi = mid;
        }
        if (lo == count) return false;
        const auto& r = h0Only ? h0Record(lo) : record(lo);
        const auto v = normalized(r.value);
        if (!flexProbeEqual(r.key, key) || v.sample != target) return false;
        out = FlexProbeRecord{r.key, v}; return true;
    };
    return search(sample) || (sample != 0 && search(0));
}

bool FlexHashCacheStorage::verify(std::string* error) const {
    FlexProbeKey last[2] {}; bool seen[2] = {false, false};
    uint64_t unique[2] = {0, 0}, h0 = 0, h1 = 0; bool h1x2 = false;
    for (uint64_t i = 0; i < recordCount_; ++i) {
        const auto& r = records_[i]; const int t = tierOf(r.value);
        h1x2 |= r.value.cacheClass == 4;
        if (t < 0) continue;
        if (t == 0) {
            if (h0 >= h0Count_ || h0Indices_[h0] != i) return fail(error, "H0 index differs from records");
            ++h0;
        } else ++h1;
        if (seen[t] && flexProbeEqual(last[t], r.key)) continue;
        seen[t] = true; last[t] = r.key; ++unique[t];
        const auto* value = lookup(t, cbqKey(r.key)); const auto expected = normalized(r.value);
        if (!value || memcmp(value, &expected, sizeof(expected)))
            return fail(error, "khash payload differs from record " + std::to_string(i));
    }
    if (h0 != h0Count_ || h1 != h1Count_ || h1x2 != hasH1X2_ ||
        unique[0] != maps_[0].size || unique[1] != maps_[1].size)
        return fail(error, "khash table counts differ from records");
    return true;
}

bool FlexHashCacheStorage::writeSnapshot(const std::string& path, std::string* error) const {
    if (!mapping_) return fail(error, "no cache loaded");
    SnapshotHeader h {}; memcpy(h.magic, snapshotMagic, 8);
    h.version = 1; h.headerBytes = sizeof(h); h.endian = 0x01020304;
    h.keyEncoding = 1; h.hashAlgorithm = 1; h.sourceVersion = sourceVersion_;
    h.flags = hasH1X2_ ? 1 : 0; h.recordBytes = 24;
    h.recordCount = recordCount_; h.h0Count = h0Count_; h.h1Count = h1Count_;
    uint64_t next = alignPage(sizeof(h));
    auto section = [&](uint64_t bytes) { const uint64_t offset = next; next = alignPage(next + bytes); return offset; };
    h.recordsOffset = section(recordCount_ * 24); h.h0Offset = section(h0Count_ * 8);
    for (unsigned t = 0; t < 2; ++t) {
        auto& d = h.tables[t]; const auto& m = maps_[t]; d.buckets = m.n_buckets; d.entries = m.size;
        d.flagsOffset = section(flagBytes(d.buckets));
        d.keysOffset = section(uint64_t(d.buckets) * sizeof(FlexProbeKey));
        d.valuesOffset = section(uint64_t(d.buckets) * sizeof(FlexProbeValue));
    }
    h.fileBytes = next;
    std::string temporary = path + ".tmp.XXXXXX";
    std::vector<char> name(temporary.begin(), temporary.end()); name.push_back(0);
    const int fd = mkstemp(name.data());
    if (fd < 0) return fail(error, "cannot create khash output: " + std::string(strerror(errno)));
    bool ok = ftruncate(fd, static_cast<off_t>(h.fileBytes)) == 0 &&
        writeAt(fd, h.recordsOffset, records_, recordCount_ * 24) &&
        writeAt(fd, h.h0Offset, h0Indices_, h0Count_ * 8);
    for (unsigned t = 0; ok && t < 2; ++t) {
        const auto& d = h.tables[t]; const auto& m = maps_[t];
        ok = writeAt(fd, d.flagsOffset, m.flags, flagBytes(d.buckets)) &&
             writeAt(fd, d.keysOffset, m.keys, uint64_t(d.buckets) * sizeof(FlexProbeKey)) &&
             writeAt(fd, d.valuesOffset, m.vals, uint64_t(d.buckets) * sizeof(FlexProbeValue));
    }
    if (ok) ok = writeAt(fd, 0, &h, sizeof(h)) && fsync(fd) == 0;
    int saved = errno;
    if (::close(fd) != 0 && ok) { ok = false; saved = errno; }
    // Atomic publication with no replacement of an existing paired artifact.
    if (ok && link(name.data(), path.c_str()) != 0) { ok = false; saved = errno; }
    unlink(name.data());
    return ok || fail(error, "cannot write khash cache: " + std::string(strerror(saved)));
}

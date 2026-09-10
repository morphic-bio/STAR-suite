#ifndef FLEX_HASH_CACHE_STORAGE_H
#define FLEX_HASH_CACHE_STORAGE_H

#include <cstdint>
#include <string>
#include <vector>
#include "klib/khash.h"

// Full 100-bit probe identity. Hash collisions never stand in for key equality.
struct FlexProbeKey { uint64_t lo, hi; };
struct FlexProbeValue {
    uint32_t geneAndRegion;
    uint8_t cacheClass, negativeCode;
    uint16_t sample;
};
struct FlexProbeRecord { FlexProbeKey key; FlexProbeValue value; };
static_assert(sizeof(FlexProbeKey) == 16, "probe key layout");
static_assert(sizeof(FlexProbeValue) == 8, "probe value layout");
static_assert(sizeof(FlexProbeRecord) == 24, "probe record layout");

// Persisted hash algorithm 1. Avalanche both words before khash's power-of-two
// bucket mask; the old unordered_map hash's low bits alone are insufficient.
inline khint_t flexProbeHash(FlexProbeKey key) {
    uint64_t h = key.lo ^ (key.hi * UINT64_C(0x9e3779b97f4a7c15));
    h ^= h >> 30; h *= UINT64_C(0xbf58476d1ce4e5b9);
    h ^= h >> 27; h *= UINT64_C(0x94d049bb133111eb);
    h ^= h >> 31;
    return static_cast<khint_t>(h);
}
inline bool flexProbeEqual(FlexProbeKey a, FlexProbeKey b) {
    return a.lo == b.lo && a.hi == b.hi;
}
KHASH_INIT(flex_probe, FlexProbeKey, FlexProbeValue, 1, flexProbeHash, flexProbeEqual)

// Owns either legacy record mmap + built khash arrays, or one read-only mmap
// containing records, H0 record indices and both already-built khash tables.
// The mapped tables are immutable: never pass them to kh_put/resize/destroy.
class FlexHashCacheStorage {
public:
    FlexHashCacheStorage() = default;
    ~FlexHashCacheStorage();
    FlexHashCacheStorage(const FlexHashCacheStorage&) = delete;
    FlexHashCacheStorage& operator=(const FlexHashCacheStorage&) = delete;

    bool open(const std::string& path, std::string* error);
    bool writeSnapshot(const std::string& path, std::string* error) const;
    // Offline validation: compare every distinct tier key/payload with records.
    bool verify(std::string* error) const;
    bool persisted() const { return persisted_; }
    uint16_t sourceVersion() const { return sourceVersion_; }
    uint64_t recordCount() const { return recordCount_; }
    uint64_t h0Count() const { return h0Count_; }
    uint64_t h1Count() const { return h1Count_; }
    bool hasH1X2() const { return hasH1X2_; }
    uint64_t tableSize(unsigned tier) const { return maps_[tier].size; }
    uint64_t tableBuckets(unsigned tier) const { return maps_[tier].n_buckets; }
    const FlexProbeRecord& record(uint64_t i) const { return records_[i]; }
    const FlexProbeRecord& h0Record(uint64_t i) const { return records_[h0Indices_[i]]; }
    const FlexProbeValue* lookup(unsigned tier, FlexProbeKey key) const {
        const auto* h = &maps_[tier];
        const khiter_t i = kh_get(flex_probe, h, key);
        return i == kh_end(h) ? nullptr : &kh_val(h, i);
    }
    bool find(FlexProbeKey key, uint16_t sample, bool h0Only, FlexProbeRecord& out) const;
    static FlexProbeKey cbqKey(FlexProbeKey key);

private:
    void close();
    bool build(std::string* error);
    bool openSnapshot(std::string* error);
    FlexProbeValue normalized(FlexProbeValue value) const;
    void* mapping_ = nullptr;
    uint64_t mappingBytes_ = 0;
    const FlexProbeRecord* records_ = nullptr;
    const uint64_t* h0Indices_ = nullptr;
    std::vector<uint64_t> ownedH0Indices_;
    khash_t(flex_probe) maps_[2] {};
    uint64_t recordCount_ = 0, h0Count_ = 0, h1Count_ = 0;
    uint16_t sourceVersion_ = 0;
    bool persisted_ = false, hasH1X2_ = false;
};
#endif

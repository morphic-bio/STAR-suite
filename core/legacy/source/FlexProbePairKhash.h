#ifndef FLEX_PROBE_PAIR_KHASH_H
#define FLEX_PROBE_PAIR_KHASH_H

// Replace the enumerated full H1X2 table with two
// Hamming-1 half tables, retaining probe lists so their intersection is exact.
#include "FlexHashCacheStorage.h"
#include "FlexHashScreenDecision.h"
#include <cstring>
#include <algorithm>
#include <stdexcept>

inline khint_t flexHalfHash(uint64_t key) { return flexProbeHash(FlexProbeKey{key, 0}); }
KHASH_INIT(flex_half_span, uint64_t, uint64_t, 1, flexHalfHash, kh_int64_hash_equal)

class FlexProbePairKhash {
public:
    enum Route { NoMatch, ExactMatch, DoubleMatch, HalfMatch, PairAmbiguous, DifferentProbes };
    struct Outcome { FlexHashScreenDecision decision; Route route = NoMatch; };
    FlexProbePairKhash() = default;
    FlexProbePairKhash(const FlexProbePairKhash&) = delete;
    FlexProbePairKhash& operator=(const FlexProbePairKhash&) = delete;
    ~FlexProbePairKhash();
    bool open(const std::string& path, std::string* error);
    bool writeSnapshot(const std::string& path, std::string* error) const;
    bool verifySource(const FlexHashCacheStorage& source, std::string* error) const;
    static bool isSnapshot(const std::string& path);
    bool sameTables(const FlexProbePairKhash& other) const;
    bool ready() const { return probeCount_ != 0; }
    uint64_t recordCount() const { return recordCount_; }
    uint64_t h0Count() const { return h0Count_; }
    const FlexProbeValue* lookupH0(FlexProbeKey key) const {
        const auto k = kh_get(flex_probe, &h0_, key);
        return k == kh_end(&h0_) ? nullptr : &kh_val(&h0_, k);
    }
    Outcome primaryOnly(FlexProbeKey key) const { return primary(key); }
    Outcome seedOnly(FlexProbeKey key, uint64_t nmask) const { return extend(key,nmask); }
    bool find(FlexProbeKey rawKey, uint16_t sample, bool h0Only, FlexProbeRecord& out) const;

    void build(const FlexHashCacheStorage& cache) {
        if (cache.sourceVersion() != 3 || !cache.hasH1X2())
            throw std::runtime_error("half tables require an H0/H1X2 v3 cache");
        if (!probes_.empty()) throw std::runtime_error("pair hash already built");
        recordCount_ = cache.recordCount(); h0Count_ = cache.h0Count();
        rawH0_.reserve(h0Count_);
        for (uint64_t i = 0; i < cache.h0Count(); ++i) {
            const auto& r = cache.h0Record(i);
            rawH0_.push_back(r);
            const auto key = FlexHashCacheStorage::cbqKey(r.key);
            if (!probes_.empty() && flexProbeEqual(probes_.back().key, key)) continue;
            probes_.push_back(FlexProbeRecord{key, r.value});
        }
        probesData_ = probes_.data(); probeCount_ = probes_.size(); rawH0Data_ = rawH0_.data();
        if (probes_.empty()) throw std::runtime_error("pair hash needs H0 probes");
        if (kh_resize(flex_probe, &h0_, capacity(probes_.size()))) throw std::bad_alloc();
        memset(h0_.keys, 0, size_t(h0_.n_buckets)*sizeof(*h0_.keys));
        memset(h0_.vals, 0, size_t(h0_.n_buckets)*sizeof(*h0_.vals));
        for (const auto& p : probes_) {
            int ret; const auto k = kh_put(flex_probe, &h0_, p.key, &ret);
            if (ret < 0) throw std::bad_alloc();
            if (ret) kh_val(&h0_, k) = p.value;
        }
        struct Entry { uint64_t key; uint32_t parent; };
        for (unsigned side = 0; side < 2; ++side) {
            std::vector<Entry> entries; entries.reserve(probes_.size()*76);
            for (uint32_t id = 0; id < probes_.size(); ++id) {
                const uint64_t key = half(probesData_[id].key, side);
                entries.push_back(Entry{key, id});
                for (unsigned pos = 0; pos < 25; ++pos) {
                    const unsigned shift = 2*pos; const uint64_t ref = (key >> shift) & 3;
                    for (uint64_t alt = 0; alt < 4; ++alt) if (alt != ref)
                        entries.push_back(Entry{(key & ~(UINT64_C(3)<<shift)) | (alt<<shift), id});
                }
            }
            std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
                return a.key < b.key || (a.key == b.key && a.parent < b.parent);
            });
            uint64_t unique = 0;
            for (size_t i = 0; i < entries.size(); ++i) if (!i || entries[i].key != entries[i-1].key) ++unique;
            auto& h = halves_[side];
            if (kh_resize(flex_half_span, &h, capacity(unique))) throw std::bad_alloc();
            memset(h.keys, 0, size_t(h.n_buckets)*sizeof(*h.keys));
            memset(h.vals, 0, size_t(h.n_buckets)*sizeof(*h.vals));
            auto& parents = parents_[side]; parents.reserve(entries.size());
            for (size_t begin = 0; begin < entries.size();) {
                size_t end = begin+1;
                while (end < entries.size() && entries[end].key == entries[begin].key) ++end;
                const uint64_t offset = parents.size();
                for (size_t i = begin; i < end; ++i)
                    if (i == begin || entries[i].parent != entries[i-1].parent) parents.push_back(entries[i].parent);
                const uint64_t count = parents.size()-offset;
                if (parents.size() > UINT32_MAX || count > UINT32_MAX) throw std::runtime_error("half span overflow");
                int ret; const auto k = kh_put(flex_half_span, &h, entries[begin].key, &ret);
                if (ret < 0) throw std::bad_alloc();
                kh_val(&h, k) = (offset << 32) | count;
                begin = end;
            }
            parentsData_[side] = parents.data(); parentCounts_[side] = parents.size();
        }
    }

    Outcome classifyRead(const char* s, bool singleN = true) const {
        struct BaseCodes {
            uint8_t value[256];
            BaseCodes() {
                for(auto& v:value)v=4;
                value['A']=value['a']=0;value['C']=value['c']=1;
                value['G']=value['g']=2;value['T']=value['t']=3;
            }
        };
        static const BaseCodes lut;
        FlexProbeKey k {0,0}; uint64_t nmask = 0;
        for(unsigned i=0;i<32;++i) {
            uint64_t b=lut.value[static_cast<unsigned char>(s[i])];
            if(b==4){b=0;nmask|=UINT64_C(1)<<i;}
            k.lo|=b<<(2*i);
        }
        for(unsigned i=0;i<18;++i) {
            uint64_t b=lut.value[static_cast<unsigned char>(s[i+32])];
            if(b==4){b=0;nmask|=UINT64_C(1)<<(i+32);}
            k.hi|=b<<(2*i);
        }
        return classify(k,nmask,singleN);
    }
    Outcome classify(FlexProbeKey key, uint64_t nmask = 0, bool singleN = true) const {
        // Match the production single-N policy: four full-window substitutions
        // before the seed stage, not four independent relaxed half matches.
        if (singleN && nmask && !(nmask & (nmask-1))) {
            const unsigned pos = __builtin_ctzll(nmask); bool kept = false, ambiguous = false, denied = false;
            Outcome result; uint8_t cls = 0xff;
            for (uint64_t b = 0; b < 4; ++b) {
                auto k = key;
                if (pos < 32) k.lo = (k.lo & ~(UINT64_C(3)<<(2*pos))) | (b<<(2*pos));
                else k.hi = (k.hi & ~(UINT64_C(3)<<(2*(pos-32)))) | (b<<(2*(pos-32)));
                const Outcome next = primary(k);
                if (next.decision.action == FlexHashScreenDecision::Deny) denied = true;
                if (next.decision.action != FlexHashScreenDecision::Keep) continue;
                if (!kept) { result = next; kept = true; cls = next.decision.cacheClass; }
                else if (next.decision.geneIdx15 != result.decision.geneIdx15) ambiguous = true;
                else if (next.decision.cacheClass != cls) cls = 0xfe;
            }
            if (kept && !ambiguous && !denied) {
                result.decision.singleN = true; result.decision.singleNCacheClass = cls;
                result.decision.cacheClass = FlexHashCacheH1; return result;
            }
        } else if (!nmask) {
            return primary(key, true);
        }
        return extend(key, nmask);
    }

    size_t probeCount() const { return probeCount_; }
    uint64_t keys(unsigned side) const { return halves_[side].size; }
    uint64_t bytes() const {
        uint64_t n = (probeCount_+h0Count_)*sizeof(FlexProbeRecord) + uint64_t(h0_.n_buckets)*24 + __ac_fsize(h0_.n_buckets)*4;
        for (unsigned s = 0; s < 2; ++s) n += uint64_t(halves_[s].n_buckets)*16 + uint64_t(__ac_fsize(halves_[s].n_buckets))*4 + parentCounts_[s]*4;
        return n;
    }

private:
    static khint_t capacity(uint64_t n) {
        uint64_t b = 4; while (n >= uint64_t(b*__ac_HASH_UPPER+0.5)) b *= 2;
        if (b > (UINT64_C(1)<<31)) throw std::runtime_error("half hash capacity overflow");
        return static_cast<khint_t>(b);
    }
    static uint64_t half(FlexProbeKey key, unsigned side) {
        const uint64_t mask = (UINT64_C(1)<<50)-1;
        return side ? ((key.lo>>50) | (key.hi<<14)) & mask : key.lo & mask;
    }
    static unsigned distance(uint64_t a, uint64_t b, uint32_t nmask = 0) {
        uint64_t d = a ^ b;
        while (nmask) { const unsigned p = __builtin_ctz(nmask); d |= UINT64_C(1)<<(2*p); nmask &= nmask-1; }
        return __builtin_popcountll((d | (d>>1)) & UINT64_C(0x5555555555555));
    }
    struct Span { const uint32_t* ids = nullptr; uint32_t count = 0; };
    Span lookup(unsigned side, uint64_t key) const {
        const auto& h = halves_[side]; const auto k = kh_get(flex_half_span, &h, key);
        Span s; if (k == kh_end(&h)) return s;
        const uint64_t v = kh_val(&h, k); s.ids = parentsData_[side]+(v>>32); s.count = uint32_t(v); return s;
    }
    Outcome primary(FlexProbeKey key, bool finish = false) const {
        Outcome out; out.decision.action = FlexHashScreenDecision::Pass;
        const auto k = kh_get(flex_probe, &h0_, key);
        if (k != kh_end(&h0_)) {
            const auto& v = kh_val(&h0_, k); out.route = ExactMatch;
            out.decision.action = v.negativeCode == FlexHashNegProbeAmbig ? FlexHashScreenDecision::Deny : FlexHashScreenDecision::Keep;
            out.decision.geneIdx15 = out.decision.action == FlexHashScreenDecision::Keep ? v.geneAndRegion & 0x7fff : 0;
            out.decision.cacheClass = 0; out.decision.negativeCode = v.negativeCode == FlexHashNegProbeAmbig ? v.negativeCode : 0;
            out.decision.probeRegion = static_cast<FlexGdnaRegion>(v.geneAndRegion >> 30); return out;
        }
        const Span l = lookup(0, half(key,0)), r = lookup(1,half(key,1));
        uint32_t i = 0, j = 0, count = 0, parent = 0;
        while (i < l.count && j < r.count) {
            if (l.ids[i] < r.ids[j]) ++i;
            else if (r.ids[j] < l.ids[i]) ++j;
            else { parent = l.ids[i]; ++i; ++j; if (++count > 1) break; }
        }
        if (count > 1) {
            out.route = PairAmbiguous; out.decision.action = FlexHashScreenDecision::Deny;
            out.decision.cacheClass = 2; out.decision.negativeCode = FlexHashNegProbeAmbig;
        } else if (count == 1) {
            const auto& p = probesData_[parent]; out.route = DoubleMatch;
            out.decision.action = FlexHashScreenDecision::Keep; out.decision.geneIdx15 = p.value.geneAndRegion & 0x7fff;
            out.decision.cacheClass = 4; out.decision.probeRegion = static_cast<FlexGdnaRegion>(p.value.geneAndRegion >> 30);
        }
        if (!count && finish) {
            Seed left, right;
            left.count = std::min(l.count, 2u); right.count = std::min(r.count, 2u);
            if (l.count) left.id = l.ids[0]; if (r.count) right.id = r.ids[0];
            return score(key, 0, left, right);
        }
        return out;
    }
    struct Seed { uint32_t id = 0, count = 0; };
    Seed seed(unsigned side, uint64_t key, uint32_t nmask) const {
        Seed result;
        if (nmask && (nmask & (nmask-1))) return result;
        auto merge = [&](Span s, uint64_t query, bool exactOnly) {
            for (uint32_t i = 0; i < s.count; ++i) {
                const auto id = s.ids[i];
                if (exactOnly && half(probesData_[id].key, side) != query) continue;
                if (result.count && result.id != id) result.count = 2;
                else if (!result.count) { result.count = 1; result.id = id; }
            }
        };
        if (!nmask) { const Span s = lookup(side,key); result.count = std::min(s.count,2u); if(s.count)result.id=s.ids[0]; }
        else {
            const unsigned shift = 2*__builtin_ctz(nmask);
            for (uint64_t b = 0; b < 4; ++b) {
                const auto q = (key & ~(UINT64_C(3)<<shift)) | (b<<shift);
                merge(lookup(side,q),q,true);
            }
        }
        return result;
    }
    Outcome extend(FlexProbeKey key, uint64_t nmask) const {
        const uint64_t lkey = half(key,0), rkey = half(key,1);
        const uint32_t ln = nmask & ((1u<<25)-1), rn = (nmask>>25) & ((1u<<25)-1);
        return score(key, nmask, seed(0,lkey,ln), seed(1,rkey,rn));
    }
    Outcome score(FlexProbeKey key, uint64_t nmask, Seed l, Seed r) const {
        const uint64_t lkey = half(key,0), rkey = half(key,1);
        const uint32_t ln = nmask & ((1u<<25)-1), rn = (nmask>>25) & ((1u<<25)-1);
        Outcome out; out.decision.action = FlexHashScreenDecision::Deny;
        if (l.count > 1 || r.count > 1) { out.route=PairAmbiguous;out.decision.negativeCode=FlexHashNegHalfGeneAmbig;return out; }
        if (!l.count && !r.count) { out.decision.negativeCode=FlexHashNegHalfNoAnchor;return out; }
        if (l.count && r.count && l.id != r.id) { out.route=DifferentProbes;out.decision.negativeCode=FlexHashNegHalfSplitProbe;return out; }
        // HALF_MATCH is explicit here: it carries the parent directly into
        // full-probe Hamming scoring, with no second half-hash search.
        out.route = l.count && r.count ? DoubleMatch : HalfMatch;
        const auto& p = probesData_[l.count ? l.id : r.id];
        const unsigned d = distance(lkey,half(p.key,0),ln) + distance(rkey,half(p.key,1),rn);
        out.decision.probeHammingDistance = d;
        if (d > 10) { out.decision.negativeCode=FlexHashNegHalfScoreFail;return out; }
        out.decision.action=FlexHashScreenDecision::Keep;out.decision.cacheClass=4;
        out.decision.geneIdx15=p.value.geneAndRegion & 0x7fff;return out;
    }
    void close();
    void* mapping_ = nullptr;
    uint64_t mappingBytes_ = 0;
    const FlexProbeRecord* probesData_ = nullptr;
    const FlexProbeRecord* rawH0Data_ = nullptr;
    const uint32_t* parentsData_[2] {};
    uint64_t probeCount_ = 0, recordCount_ = 0, h0Count_ = 0, parentCounts_[2] {};
    std::vector<FlexProbeRecord> rawH0_;
    std::vector<FlexProbeRecord> probes_;
    std::vector<uint32_t> parents_[2];
    khash_t(flex_probe) h0_ {};
    khash_t(flex_half_span) halves_[2] {};
};
#endif

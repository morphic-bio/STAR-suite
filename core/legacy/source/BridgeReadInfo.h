#ifndef STAR_BRIDGE_READ_INFO_H
#define STAR_BRIDGE_READ_INFO_H

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

// Optional read identities for consumers of a direct bridge's final CB/UMI
// decisions (Velocyto). No per-read allocations or state on GEX-only runs.
class BridgeReadInfo {
public:
    struct Read {
        uint32_t cb, umi, gene, readId;
    };
    struct Pending {
        uint64_t cbKey;
        Read read;
    };
    struct Decision {
        uint32_t gene, original, corrected;
    };
    std::vector<Read> reads;
    std::vector<Pending> pending;
    std::vector<std::pair<uint64_t, uint32_t>> resolvedCBs;

    void append(BridgeReadInfo &other) {
        reads.insert(reads.end(), other.reads.begin(), other.reads.end());
        pending.insert(pending.end(), other.pending.begin(), other.pending.end());
        std::vector<Read>().swap(other.reads);
        std::vector<Pending>().swap(other.pending);
    }

    // Resolve with the same key-level CB decisions that populated the count
    // hash. Orphan/unresolved ambiguous keys remain absent. Match CountingSink:
    // reads without an assignment to the selected gene feature are excluded.
    template<class Allowed>
    void prepare(Allowed allowed) {
        std::sort(resolvedCBs.begin(), resolvedCBs.end());
        for (auto &p : pending) {
            auto it = std::lower_bound(resolvedCBs.begin(), resolvedCBs.end(),
                                      std::make_pair(p.cbKey, uint32_t{0}));
            if (it != resolvedCBs.end() && it->first == p.cbKey) {
                p.read.cb = it->second;
                reads.push_back(p.read);
            }
        }
        std::vector<Pending>().swap(pending);
        decltype(resolvedCBs)().swap(resolvedCBs);
        size_t kept = 0;
        for (const Read &r : reads) {
            if (!allowed(r.cb) || r.gene == UINT32_MAX)
                continue;
            reads[kept++] = r;
        }
        reads.resize(kept);
        std::sort(reads.begin(), reads.end(), [](const Read &a, const Read &b) {
            return a.cb < b.cb;
        });
    }

    // Called once per CB by its collapse worker. Distinct CBs own disjoint
    // read IDs, and PackedReadInfo stores one independent word per read.
    template<class Sink>
    void replay(uint32_t cb, std::vector<Decision> &decisions, Sink sink) const {
        std::sort(decisions.begin(), decisions.end(), [](const Decision &a, const Decision &b) {
            return a.gene < b.gene || (a.gene == b.gene && a.original < b.original);
        });
        auto r = std::lower_bound(reads.begin(), reads.end(), cb,
                                 [](const Read &a, uint32_t b) { return a.cb < b; });
        for (; r != reads.end() && r->cb == cb; ++r) {
            const auto key = std::make_pair(r->gene, r->umi);
            auto d = std::lower_bound(decisions.begin(), decisions.end(), key,
                [](const Decision &a, const std::pair<uint32_t, uint32_t> &b) {
                    return std::make_pair(a.gene, a.original) < b;
                });
            const bool accepted = d != decisions.end() && d->gene == r->gene
                && d->original == r->umi && d->corrected != UINT32_MAX;
            sink(r->readId, cb, accepted ? d->corrected : UINT32_MAX,
                 accepted ? uint8_t{1} : uint8_t{2});
        }
    }
};
#endif

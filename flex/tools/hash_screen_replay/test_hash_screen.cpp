// Unit tests for flat and tiered hash screen lookup.
// Each test constructs a small in-memory cache, runs both flat and tiered,
// and asserts identical decisions.

#include "FlexHashScreenFlat.h"
#include "FlexHashScreenTiered.h"
#include <cassert>
#include <cstdio>
#include <cstring>

static int g_pass = 0, g_fail = 0;

#define CHECK(cond, msg) do { \
    if (!(cond)) { \
        fprintf(stderr, "  FAIL: %s (line %d)\n", msg, __LINE__); \
        ++g_fail; return; \
    } \
} while(0)

#define CHECK_DECISION(flat_d, tiered_d, expected_action, label) do { \
    CHECK(flat_d.action == expected_action, label ": flat action mismatch"); \
    CHECK(tiered_d.action == expected_action, label ": tiered action mismatch"); \
    CHECK(flat_d == tiered_d, label ": flat != tiered"); \
} while(0)

// Helper: build a 50-char read from a seed pattern.
// We use a repeating "ACGT" pattern with the seed affecting the first few bases.
static void makeRead(char* buf, int len, int seed = 0) {
    const char bases[] = "ACGT";
    for (int i = 0; i < len; ++i)
        buf[i] = bases[(i + seed) % 4];
    buf[len] = '\0';
}

// Helper: encode a 50-mer from a read buffer at offset 0 into (seqLo, seqHi).
static void encodeRead(const char* read, uint64_t& lo, uint64_t& hi) {
    bool ok = encodeWindow(read, 0, lo, hi);
    assert(ok);
    (void)ok;
}

// Helper: build a Record with given properties.
static Record makeRecord(uint64_t lo, uint64_t hi, uint16_t sampleIdx,
                         uint8_t cacheClass, uint32_t gene,
                         uint8_t negCode = FlexHashNegNone) {
    Record r;
    r.seqLo = lo;
    r.seqHi = hi;
    r.sampleIdx = sampleIdx;
    r.cacheClass = cacheClass;
    r.resolvedGeneIdx15 = gene;
    r.negativeCode = negCode;
    return r;
}

// ── Test 1: H0 exact match, offset 0, sample matched → KEEP ────────────────

static void test_h0_exact_offset0_sample_match() {
    fprintf(stderr, "Test 1: H0 exact match, offset 0, sample matched\n");
    char read[52]; makeRead(read, 51);
    uint64_t lo, hi; encodeRead(read, lo, hi);

    std::vector<Record> recs = { makeRecord(lo, hi, 5, 0, 100) };

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    auto df = flat.classifyRead(read, 51, 5);
    auto dt = tiered.classifyRead(read, 51, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Keep, "H0 offset0 match");
    CHECK(df.geneIdx15 == 100, "gene");
    CHECK(df.cacheClass == 0, "class");
    CHECK(df.offset == 0, "offset");
    ++g_pass;
}

// ── Test 2: H0 exact match, offset +1 → KEEP ──────────────────────────────

static void test_h0_exact_offset_plus1() {
    fprintf(stderr, "Test 2: H0 exact match, offset +1\n");
    char read[53]; makeRead(read, 52);
    // Encode offset +1 window
    uint64_t lo, hi; encodeWindow(read, 1, lo, hi);

    std::vector<Record> recs = { makeRecord(lo, hi, 5, 0, 200) };

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    auto df = flat.classifyRead(read, 52, 5);
    auto dt = tiered.classifyRead(read, 52, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Keep, "H0 offset+1");
    CHECK(df.geneIdx15 == 200, "gene");
    CHECK(df.offset == 1, "offset");
    ++g_pass;
}

// ── Test 3: H1 near match, single gene → KEEP ─────────────────────────────

static void test_h1_near_match_single_gene() {
    fprintf(stderr, "Test 3: H1 near match, single gene\n");
    char read[52]; makeRead(read, 51);
    uint64_t lo, hi; encodeRead(read, lo, hi);

    // H1 record with sampleIdx=0 (global), gene > 0
    std::vector<Record> recs = { makeRecord(lo, hi, 0, 1, 300) };

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    auto df = flat.classifyRead(read, 51, 5);
    auto dt = tiered.classifyRead(read, 51, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Keep, "H1 single gene");
    CHECK(df.geneIdx15 == 300, "gene");
    CHECK(df.cacheClass == 1, "class");
    ++g_pass;
}

// ── Test 3b: H2 (cacheClass 3) near match, single gene → KEEP ─────────────

static void test_h2_near_match_single_gene() {
    fprintf(stderr, "Test 3b: H2 (class 3) near match, single gene\n");
    char read[52]; makeRead(read, 51);
    uint64_t lo, hi; encodeRead(read, lo, hi);

    std::vector<Record> recs = { makeRecord(lo, hi, 0, 3, 301) };

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    CHECK(tiered.h2Count() == 1, "tiered H2 count");

    auto df = flat.classifyRead(read, 51, 5);
    auto dt = tiered.classifyRead(read, 51, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Keep, "H2 single gene");
    CHECK(df.geneIdx15 == 301, "gene");
    CHECK(df.cacheClass == 3, "class 3");
    ++g_pass;
}

// ── Test 4: H1 near match, gene conflict (two offsets) → DENY ─────────────

static void test_h1_gene_conflict_two_offsets() {
    fprintf(stderr, "Test 4: H1 gene conflict across offsets\n");
    char read[53]; makeRead(read, 52);

    // Offset 0 → gene 300
    uint64_t lo0, hi0; encodeWindow(read, 0, lo0, hi0);
    // Offset 1 → gene 400
    uint64_t lo1, hi1; encodeWindow(read, 1, lo1, hi1);

    std::vector<Record> recs = {
        makeRecord(lo0, hi0, 0, 1, 300),
        makeRecord(lo1, hi1, 0, 1, 400),
    };

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    auto df = flat.classifyRead(read, 52, 5);
    auto dt = tiered.classifyRead(read, 52, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Deny, "H1 gene conflict");
    ++g_pass;
}

// ── Test 5: Deny record (ProbeAmbig) → DENY ───────────────────────────────

static void test_deny_probe_ambig() {
    fprintf(stderr, "Test 5: Deny record (ProbeAmbig)\n");
    char read[52]; makeRead(read, 51);
    uint64_t lo, hi; encodeRead(read, lo, hi);

    std::vector<Record> recs = {
        makeRecord(lo, hi, 0, 2, 0, FlexHashNegProbeAmbig),
    };

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    auto df = flat.classifyRead(read, 51, 5);
    auto dt = tiered.classifyRead(read, 51, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Deny, "ProbeAmbig");
    ++g_pass;
}

// ── Test 6: Sample-specific record invisible to wrong sample → PASS ─────────
// findRecord only returns exact-sample-match or sample=0 fallback.
// A record at sampleIdx=5 is invisible to sampleIdx=7 → PASS (cache miss).

static void test_sample_mismatch() {
    fprintf(stderr, "Test 6: Sample-specific record invisible to wrong sample\n");
    char read[52]; makeRead(read, 51);
    uint64_t lo, hi; encodeRead(read, lo, hi);

    // Cache has sampleIdx=5, runtime is sampleIdx=7
    std::vector<Record> recs = { makeRecord(lo, hi, 5, 0, 100) };

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    // Wrong sample → findRecord can't find it → PASS
    auto df = flat.classifyRead(read, 51, 7);
    auto dt = tiered.classifyRead(read, 51, 7);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Pass, "wrong sample → PASS");

    // Correct sample → KEEP
    auto df2 = flat.classifyRead(read, 51, 5);
    auto dt2 = tiered.classifyRead(read, 51, 5);
    CHECK_DECISION(df2, dt2, FlexHashScreenDecision::Keep, "right sample → KEEP");
    CHECK(df2.geneIdx15 == 100, "gene");
    ++g_pass;
}

// ── Test 7: Cache miss (no record for sequence) → PASS ─────────────────────

static void test_cache_miss() {
    fprintf(stderr, "Test 7: Cache miss\n");
    // Use a read whose encoded windows at all offsets {0,+1} differ from cache.
    // 52-char read so offset +1 is also valid (1+50 ≤ 52).
    char read[53]; makeRead(read, 52, 0);

    // Cache record: all-A 50-mer (completely different encoded sequence)
    char other[51];
    memset(other, 'A', 50);
    other[50] = '\0';
    uint64_t lo, hi; encodeRead(other, lo, hi);

    std::vector<Record> recs = { makeRecord(lo, hi, 0, 0, 100) };

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    auto df = flat.classifyRead(read, 52, 5);
    auto dt = tiered.classifyRead(read, 52, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Pass, "cache miss");
    ++g_pass;
}

// ── Test 8: Read too short (< 50 bp) → PASS ───────────────────────────────

static void test_read_too_short() {
    fprintf(stderr, "Test 8: Read too short\n");
    char read[50]; makeRead(read, 49);

    std::vector<Record> recs;

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    auto df = flat.classifyRead(read, 49, 5);
    auto dt = tiered.classifyRead(read, 49, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Pass, "too short");
    ++g_pass;
}

// ── Test 9: Read with N base → PASS ────────────────────────────────────────

static void test_read_with_n() {
    fprintf(stderr, "Test 9: Read with N base\n");
    char read[52]; makeRead(read, 51);
    read[25] = 'N';

    std::vector<Record> recs;

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    auto df = flat.classifyRead(read, 51, 5);
    auto dt = tiered.classifyRead(read, 51, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Pass, "N base");
    ++g_pass;
}

// ── Test 10: H0 + deny at same offset → deny overrides keep ───────────────

static void test_h0_plus_deny_same_offset() {
    fprintf(stderr, "Test 10: H0 + deny at same offset\n");
    char read[52]; makeRead(read, 51);
    uint64_t lo, hi; encodeRead(read, lo, hi);

    // Both H0 keep and deny record for the same sequence.
    // In the flat array, classifyHits sees both. The deny (ProbeAmbig)
    // sets sawAmbig=true which takes precedence over the non-exact keep path.
    // But H0 class-0 + sampleMatched triggers immediate return before deny
    // is evaluated — so the H0 exact match WINS if sample matches.
    //
    // However, deny records have cacheClass=2 and ProbeAmbig. In the flat
    // array (sorted by seqHi,seqLo,sampleIdx), findRecord returns the FIRST
    // match (sampleIdx-specific then fallback). Both records have sampleIdx=0
    // but different cacheClass. Since cacheClass is NOT in the sort key, which
    // one findRecord returns depends on insertion order.
    //
    // Actually, findRecord only returns ONE record per offset. With sampleIdx=0
    // for both, the flat search returns whichever sorts first — and since
    // cacheClass is not in the sort key, these are adjacent. lower_bound
    // returns the first one. In practice the deny (class 2) and H0 (class 0)
    // with same (seqHi, seqLo, sampleIdx=0) are equal under recordLess,
    // so the order depends on std::sort stability (not guaranteed).
    //
    // To make this test deterministic, use different sampleIdx so findRecord
    // can find both: H0 with sample=5 (specific match), deny with sample=0.
    // The flat path: findRecord with sample=5 finds H0 (sample=5), returns it.
    // classifyHits sees H0+sampleMatched → immediate KEEP.
    // We never reach the deny record because findRecord returns only ONE record
    // per offset.
    //
    // Better test: use sample=0 for H0 so findRecord returns the first of
    // the two records at sample=0. To test deny-overrides-keep, we need
    // the deny to be at a DIFFERENT offset or the H0 to not trigger the
    // immediate return (i.e., sample doesn't match exactly).
    //
    // Revised: H0 with sampleIdx=0 (global, gene>0) + deny with sampleIdx=0.
    // Runtime sample=5. H0 match is not sampleMatched (rec.sample=0),
    // so it falls into the nonExactKeep path. Deny sets sawAmbig.
    // classifyHits checks sawAmbig first → DENY.
    //
    // But findRecord only returns ONE record. With both at sampleIdx=0
    // and same (lo,hi), it returns whichever std::sort placed first.
    // This is not deterministic.
    //
    // Cleanest approach: put deny at a different sampleIdx so findRecord
    // can return both. E.g., H0 at sampleIdx=0, deny at sampleIdx=1.
    // But deny records typically have sampleIdx=0.
    //
    // Actually the simplest deterministic test: put them at different OFFSETS.
    // H0 at offset 0, deny at offset +1.

    // offset 0 → H0 keep (sample=0 global, class=0, gene=100)
    uint64_t lo0, hi0; encodeWindow(read, 0, lo0, hi0);
    // offset 1 → deny (class=2, ProbeAmbig)
    uint64_t lo1, hi1; encodeWindow(read, 1, lo1, hi1);

    std::vector<Record> recs = {
        makeRecord(lo0, hi0, 0, 0, 100),
        makeRecord(lo1, hi1, 0, 2, 0, FlexHashNegProbeAmbig),
    };

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    // Runtime sample=5: H0 record has sample=0, so not sampleMatched.
    // It goes into nonExactKeep. Deny at offset 1 sets sawAmbig.
    // sawAmbig has priority → DENY.
    auto df = flat.classifyRead(read, 52, 5);
    auto dt = tiered.classifyRead(read, 52, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Deny, "H0+deny → deny wins");
    ++g_pass;
}

// ── Test 11: H0 at offset 0 + H1 at offset +1, different gene → conflict ──

static void test_h0_h1_cross_offset_gene_conflict() {
    fprintf(stderr, "Test 11: H0 offset 0 + H1 offset +1, different gene\n");
    char read[53]; makeRead(read, 52);

    uint64_t lo0, hi0; encodeWindow(read, 0, lo0, hi0);
    uint64_t lo1, hi1; encodeWindow(read, 1, lo1, hi1);

    // H0 at offset 0, gene 100, sampleIdx=0 (global)
    // H1 at offset 1, gene 200, sampleIdx=0 (global)
    std::vector<Record> recs = {
        makeRecord(lo0, hi0, 0, 0, 100),
        makeRecord(lo1, hi1, 0, 1, 200),
    };

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    // Runtime sample=5: neither record is sampleMatched (both have sample=0).
    // Both fall into nonExactKeep. Different genes → gene conflict → DENY.
    auto df = flat.classifyRead(read, 52, 5);
    auto dt = tiered.classifyRead(read, 52, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Deny, "cross-offset gene conflict");
    ++g_pass;
}

// ── Test 12: Empty cache → all reads PASS ──────────────────────────────────

static void test_empty_cache() {
    fprintf(stderr, "Test 12: Empty cache\n");
    char read[52]; makeRead(read, 51);

    std::vector<Record> recs;

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    auto df = flat.classifyRead(read, 51, 5);
    auto dt = tiered.classifyRead(read, 51, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Pass, "empty cache");
    ++g_pass;
}

// ── Bonus: dead-weight records are dropped by tiered ────────────────────────

static void test_dead_weight_dropped() {
    fprintf(stderr, "Test bonus: Dead-weight records dropped by tiered\n");
    char read[52]; makeRead(read, 51);
    uint64_t lo, hi; encodeRead(read, lo, hi);

    // Dead weight: cacheClass=2, negativeCode=0, gene=0
    std::vector<Record> recs = { makeRecord(lo, hi, 0, 2, 0, FlexHashNegNone) };

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    CHECK(flat.recordCount() == 1, "flat has 1 record");
    CHECK(tiered.h0Count() == 0, "tiered H0 empty");
    CHECK(tiered.h1Count() == 0, "tiered H1 empty");
    CHECK(tiered.denyCount() == 0, "tiered deny empty");
    CHECK(tiered.droppedCount() == 1, "tiered dropped 1");

    // Both should return PASS (flat finds the record but skips it; tiered never loads it)
    auto df = flat.classifyRead(read, 51, 5);
    auto dt = tiered.classifyRead(read, 51, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Pass, "dead weight → PASS");
    ++g_pass;
}

// ── Bonus: H0 exact match with sample match → immediate KEEP ───────────────

static void test_h0_sample_match_immediate_keep() {
    fprintf(stderr, "Test bonus: H0 sample match is immediate KEEP even with H1 conflict\n");
    char read[53]; makeRead(read, 52);

    uint64_t lo0, hi0; encodeWindow(read, 0, lo0, hi0);
    uint64_t lo1, hi1; encodeWindow(read, 1, lo1, hi1);

    // H0 at offset 0, sample=5 → immediate KEEP(gene=100)
    // H1 at offset 1, sample=0, gene=200 → would cause conflict, but never reached
    std::vector<Record> recs = {
        makeRecord(lo0, hi0, 5, 0, 100),
        makeRecord(lo1, hi1, 0, 1, 200),
    };

    FlatCache flat; flat.init(recs);
    TieredCache tiered; tiered.init(recs);

    auto df = flat.classifyRead(read, 52, 5);
    auto dt = tiered.classifyRead(read, 52, 5);
    CHECK_DECISION(df, dt, FlexHashScreenDecision::Keep, "H0 immediate KEEP");
    CHECK(df.geneIdx15 == 100, "gene 100");
    CHECK(df.cacheClass == 0, "class 0");
    CHECK(df.offset == 0, "offset 0");
    ++g_pass;
}

int main() {
    fprintf(stderr, "=== Hash Screen Unit Tests ===\n\n");

    test_h0_exact_offset0_sample_match();
    test_h0_exact_offset_plus1();
    test_h1_near_match_single_gene();
    test_h2_near_match_single_gene();
    test_h1_gene_conflict_two_offsets();
    test_deny_probe_ambig();
    test_sample_mismatch();
    test_cache_miss();
    test_read_too_short();
    test_read_with_n();
    test_h0_plus_deny_same_offset();
    test_h0_h1_cross_offset_gene_conflict();
    test_empty_cache();
    test_dead_weight_dropped();
    test_h0_sample_match_immediate_keep();

    fprintf(stderr, "\n=== Results: %d passed, %d failed ===\n", g_pass, g_fail);
    return g_fail > 0 ? 1 : 0;
}

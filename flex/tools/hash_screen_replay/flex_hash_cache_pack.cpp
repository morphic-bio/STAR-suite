#include "FlexHashCacheStorage.h"
#include "FlexProbePairKhash.h"
#include <chrono>
#include <iostream>

int main(int argc, char** argv) {
    const bool half = argc==4 && std::string(argv[1])=="--half";
    if (argc != 3 && !half) {
        std::cerr << "Usage: flex_hash_cache_pack [--half] sorted_cache.bin output.khash\n"
                     "Build and persist H0/H1 tables once; output must not exist.\n";
        return 2;
    }
    using Clock = std::chrono::steady_clock;
    if (half) { --argc; ++argv; }
    const auto start = Clock::now();
    std::string error;
    if (half) {
        FlexHashCacheStorage source; FlexProbePairKhash pair;
        if (!source.open(argv[1],&error)) { std::cerr<<error<<'\n';return 1; }
        try { pair.build(source); } catch(const std::exception& e) { std::cerr<<e.what()<<'\n';return 1; }
        std::cerr<<"Built half tables: "<<pair.probeCount()<<" probes; verifying both directions against "<<source.recordCount()<<" source records\n";
        if (!pair.verifySource(source,&error) || !pair.writeSnapshot(argv[2],&error)) { std::cerr<<error<<'\n';return 1; }
        FlexProbePairKhash stored;
        if(!stored.open(argv[2],&error) || !pair.sameTables(stored)) { std::cerr<<"Stored half table verification failed: "<<error<<'\n';return 1; }
        std::cout<<"PASS: complete source/half equivalence and exact stored table roundtrip; "<<stored.bytes()<<" bytes in "<<std::chrono::duration<double>(Clock::now()-start).count()<<" s\n";
        return 0;
    }
    {
        FlexHashCacheStorage source;
        if (!source.open(argv[1], &error) || !source.writeSnapshot(argv[2], &error)) {
            std::cerr << error << '\n'; return 1;
        }
        std::cerr << "Built/stored " << source.recordCount() << " records, H0 keys "
                  << source.tableSize(0) << ", H1/DENY keys " << source.tableSize(1)
                  << " in " << std::chrono::duration<double>(Clock::now()-start).count() << " s\n";
    }
    // Reopen the exact on-disk tables and exhaustively compare with their
    // embedded records. This cost belongs to generation, never startup.
    FlexHashCacheStorage stored;
    if (!stored.open(argv[2], &error) || !stored.verify(&error)) {
        std::cerr << "Stored-table verification failed: " << error << '\n'; return 1;
    }
    std::cout << "PASS: persisted khash tables match every distinct tier key and payload\n";
    return 0;
}

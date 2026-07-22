#include "SoloBarcodeWhitelistLookup.h"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

int main()
{
    const std::vector<uint64> componentWhitelist {0x001122, 0x001132, 0x002233};

    assert(soloBarcodeFindSortedWhitelistIndex(0x001122, componentWhitelist) == 0);
    assert(soloBarcodeFindSortedWhitelistIndex(0x001132, componentWhitelist) == 1);
    assert(soloBarcodeFindSortedWhitelistIndex(0x002233, componentWhitelist) == 2);
    assert(soloBarcodeFindSortedWhitelistIndex(0x001123, componentWhitelist) == -1);

    // Reproduce the complex-barcode control flow that exposed the regression:
    // an exact miss is followed by H1 probes while no global hash exists.
    const uint64 observed = 0x001123;
    int64 recovered = -1;
    for (uint32 position = 0; position < 6 && recovered < 0; ++position) {
        for (uint32 alternative = 1; alternative < 4; ++alternative) {
            const uint64 candidate = observed ^ (static_cast<uint64>(alternative) << (2 * position));
            recovered = soloBarcodeFindSortedWhitelistIndex(candidate, componentWhitelist);
            if (recovered >= 0) {
                break;
            }
        }
    }
    assert(recovered >= 0);

    std::cout << "Complex barcode component whitelist lookup tests passed\n";
    return 0;
}

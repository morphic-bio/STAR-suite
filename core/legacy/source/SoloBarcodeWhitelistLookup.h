#ifndef H_SoloBarcodeWhitelistLookup
#define H_SoloBarcodeWhitelistLookup

#include "IncludeDefine.h"

#include <algorithm>
#include <vector>

// Component whitelists used by CB_UMI_Complex are sorted by
// SoloBarcode::sortWhiteList and do not have the global simple-CB hash.
inline int64 soloBarcodeFindSortedWhitelistIndex(
    uint64 packed,
    const std::vector<uint64>& whitelist)
{
    const auto it = std::lower_bound(whitelist.begin(), whitelist.end(), packed);
    if (it == whitelist.end() || *it != packed) {
        return -1;
    }
    return static_cast<int64>(it - whitelist.begin());
}

#endif

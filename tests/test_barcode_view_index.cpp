#include "BarcodeViewIndex.h"
#include <cassert>
#include <iostream>
#include <map>
#include <vector>
#include <utility>

int main() {
    std::vector<std::string> keys = {"", "A", "AC", std::string("A\0C", 3),
                                   std::string(16, 'A'), std::string(80, 'T')};
    for (unsigned i = 0; i < 100000; ++i) keys.push_back("barcode-" + std::to_string(i));
    BarcodeViewIndex index;
    index.reserve(keys.size());
    for (size_t i = 0; i < keys.size(); ++i) index.insert(keys[i], i);
    for (size_t i = 0; i < keys.size(); ++i) {
        const auto* value = index.find(std::string(keys[i]));
        assert(value && *value == i);
    }
    assert(index.find("absent") == nullptr);
    index.insert(keys[3], UINT32_MAX);
    assert(*index.find(keys[3]) == UINT32_MAX);
    assert(*index.find(BarcodeView(keys[3].data(), 1)) == 1);
    index.clear();
    assert(index.size() == 0 && index.find(keys[0]) == nullptr);
    bool rejected = false;
    try { index.reserve(SIZE_MAX); } catch (const std::length_error&) { rejected = true; }
    assert(rejected);
    // Same stable-vector compaction used by table import: includes duplicate
    // long strings and SSO strings, with keys referring only to retained slots.
    keys = {"AC", "AC", std::string(16, 'C'), "G", "AC", "G", std::string(16, 'C'), "N"};
    size_t kept = 0;
    for (size_t i = 0; i < keys.size(); ++i) {
        if (index.find(keys[i])) continue;
        if (i != kept) keys[kept] = std::move(keys[i]);
        index.insert(keys[kept], kept);
        ++kept;
    }
    keys.resize(kept);
    assert(kept == 4);
    for (size_t i = 0; i < kept; ++i) assert(*index.find(keys[i]) == i);
    std::cout << "Barcode view index PASS\n";
}

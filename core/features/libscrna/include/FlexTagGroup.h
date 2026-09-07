#ifndef STAR_SUITE_FLEX_TAG_GROUP_H
#define STAR_SUITE_FLEX_TAG_GROUP_H

#include <cstddef>
#include <string>
#include <vector>

inline bool barcodeHasFlexTag(const std::string& barcode, const std::string& tag) {
    size_t end = barcode.size();
    if (end >= 2 && barcode.compare(end - 2, 2, "-1") == 0) {
        end -= 2;
    }
    return tag.size() <= end && barcode.compare(end - tag.size(), tag.size(), tag) == 0;
}

inline bool barcodeHasAnyFlexTag(
    const std::string& barcode,
    const std::vector<std::string>& tags)
{
    if (tags.empty()) return true;
    for (const std::string& tag : tags) {
        if (barcodeHasFlexTag(barcode, tag)) return true;
    }
    return false;
}

#endif

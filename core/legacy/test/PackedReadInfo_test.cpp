#include "PackedReadInfo.h"

#include <cstdint>
#include <iostream>
#include <stdexcept>

namespace {

bool check(bool condition, const char *label) {
    if (!condition) {
        std::cerr << "FAIL: " << label << '\n';
    }
    return condition;
}

}  // namespace

int main() {
    PackedReadInfo reads;
    reads.init(4, 8, 12);
    reads.set(0, 1, 7, 3);

    bool ok = true;
    ok &= check(reads.getCB(0) == 1, "CB round trip");
    ok &= check(reads.getUMI(0) == 7, "UMI round trip");
    ok &= check(reads.getStatus(0) == 3, "status round trip");

    reads.setStatus(0, 5);
    ok &= check(reads.getStatus(0) == 5, "status update");

    uint32_t cb = 0;
    uint32_t umi = 0;
    uint8_t status = 0;
    reads.unpack(reads.pack(3, 12, 1), cb, umi, status);
    ok &= check(cb == 3 && umi == 12 && status == 1, "explicit pack/unpack");

    bool rejected = false;
    try {
        PackedReadInfo invalid;
        invalid.init(1, 1, 17);
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    ok &= check(rejected, "over-width UMI rejection");

    if (!ok) {
        return 1;
    }
    std::cout << "PackedReadInfo tests passed\n";
    return 0;
}

#include "BorrowedBarcodeIndex.h"
#include <cassert>
#include <thread>
#include <vector>

int main() {
    std::vector<std::string> barcodes;
    for (uint32_t i = 0; i < 10000; ++i) {
        std::string barcode(24, 'A');
        uint32_t value = i;
        for (size_t j = 0; j < 16; ++j) {
            barcode[j] = "ACGT"[value & 3];
            value >>= 2;
        }
        barcodes.push_back(barcode);
    }
    {
        BorrowedBarcodeIndex index(barcodes.size());
        for (uint32_t i = 0; i < barcodes.size(); ++i) index.insert(barcodes[i], i);
        index.insert(barcodes[3], 9001); // last assignment wins
        assert(index.find("missing") == UINT32_MAX);
        std::vector<std::thread> readers;
        for (int worker = 0; worker < 4; ++worker) {
            readers.emplace_back([&]() {
                for (uint32_t i = 0; i < barcodes.size(); ++i) {
                    // Different backing storage must match by content.
                    const std::string query = barcodes[i];
                    assert(index.find(query) == (i == 3 ? 9001 : i));
                }
            });
        }
        for (auto& reader : readers) reader.join();
    }
    // Index destruction must not release or alter the owner's strings.
    for (auto& barcode : barcodes) {
        assert(barcode.size() == 24);
        barcode[23] = 'T';
    }
    BorrowedBarcodeIndex empty(0);
    assert(empty.find("absent") == UINT32_MAX);
}

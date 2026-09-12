#include "MexWriter.h"
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

static std::string readFile(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    assert(input.good());
    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

int main() {
    char pattern[] = "/tmp/star_mex_csr_XXXXXX";
    const char* root = mkdtemp(pattern);
    assert(root);
    const std::string base = std::string(root) + "/";
    const std::vector<std::string> barcodes = {"empty0", "cell1", "empty2", "cell3", "empty4"};
    const std::vector<std::string> genes = {"gene0", "gene1"};
    std::vector<MexWriter::Triplet> entries;
    // A cell spans two writer blocks; empty cells occur before, between and
    // after nonempty cells. Include decimal boundaries and maximum counts.
    for (uint32_t i = 0; i < (1u << 18) + 3; ++i)
        entries.push_back({1, i % 2, i % 3 == 0 ? UINT32_MAX : i});
    entries.push_back({3, 0, 10});
    entries.push_back({3, 1, 100});
    const std::string reference = base + "reference/";
    assert(mkdir(reference.c_str(), 0700) == 0);
    assert(MexWriter::writeMex(reference, barcodes, genes, entries, -1, 1) == 0);
    for (uint32_t stride : {2u, 3u}) {
        std::vector<uint32_t> words, offsets;
        size_t entry = 0;
        for (size_t cell = 0; cell < barcodes.size(); ++cell) {
            offsets.push_back(static_cast<uint32_t>(words.size()));
            while (entry < entries.size() && entries[entry].cell_idx == cell) {
                words.push_back(entries[entry].gene_idx);
                words.push_back(entries[entry++].count);
                if (stride == 3) words.push_back(987); // ignored payload
            }
        }
        offsets.push_back(static_cast<uint32_t>(words.size()));
        for (unsigned int threads : {1u, 8u}) {
            const std::string output = base + std::to_string(stride) + "_" + std::to_string(threads) + "/";
            assert(mkdir(output.c_str(), 0700) == 0);
            assert(MexWriter::writeMexCsr(output, barcodes, genes, words, offsets, stride, threads) == 0);
            for (const char* file : {"matrix.mtx", "barcodes.tsv", "features.tsv"}) {
                assert(readFile(output + file) == readFile(reference + file));
                assert(unlink((output + file).c_str()) == 0);
            }
            assert(rmdir(output.c_str()) == 0);
        }
        auto malformed = offsets;
        malformed[2] = 1;
        assert(MexWriter::writeMexCsr(base, barcodes, genes, words, malformed, stride, 8) != 0);
        malformed = offsets;
        malformed.back() -= stride;
        assert(MexWriter::writeMexCsr(base, barcodes, genes, words, malformed, stride, 8) != 0);
    }
    for (const char* file : {"matrix.mtx", "barcodes.tsv", "features.tsv"})
        assert(unlink((reference + file).c_str()) == 0);
    assert(rmdir(reference.c_str()) == 0);
    const std::vector<uint32_t> empty, offsets(barcodes.size() + 1, 0);
    assert(MexWriter::writeMexCsr(base, barcodes, genes, empty, offsets, 2, 8) == 0);
    assert(readFile(base + "matrix.mtx") == "%%MatrixMarket matrix coordinate integer general\n%\n2 5 0\n");
    for (const char* file : {"matrix.mtx", "barcodes.tsv", "features.tsv"})
        assert(unlink((base + file).c_str()) == 0);
    assert(rmdir(root) == 0);
}

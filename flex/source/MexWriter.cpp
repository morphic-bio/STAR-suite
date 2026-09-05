#include "MexWriter.h"
#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <thread>
#include <unordered_set>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>

namespace {

struct MatrixBlock {
    size_t begin;
    size_t end;
    uint64_t fileOffset;
    uint64_t textBytes;
};

inline unsigned int decimalDigits(uint32_t value)
{
    if (value < 10u) return 1;
    if (value < 100u) return 2;
    if (value < 1000u) return 3;
    if (value < 10000u) return 4;
    if (value < 100000u) return 5;
    if (value < 1000000u) return 6;
    if (value < 10000000u) return 7;
    if (value < 100000000u) return 8;
    if (value < 1000000000u) return 9;
    return 10;
}

inline char* appendUnsigned(char* output, uint32_t value)
{
    char reverse[10];
    char* cursor = reverse + sizeof(reverse);
    do {
        *--cursor = static_cast<char>('0' + value % 10u);
        value /= 10u;
    } while (value != 0u);
    const size_t length = static_cast<size_t>(reverse + sizeof(reverse) - cursor);
    std::memcpy(output, cursor, length);
    return output + length;
}

inline bool validTriplet(const MexWriter::Triplet& triplet,
                         size_t nBarcodes,
                         size_t nFeatures)
{
    return triplet.cell_idx < nBarcodes && triplet.gene_idx < nFeatures;
}

inline uint64_t tripletTextBytes(const MexWriter::Triplet& triplet)
{
    return static_cast<uint64_t>(decimalDigits(triplet.gene_idx + 1u))
        + decimalDigits(triplet.cell_idx + 1u)
        + decimalDigits(triplet.count)
        + 3u;
}

template <typename Function>
void parallelBlocks(size_t nBlocks, unsigned int requestedThreads, Function function)
{
    if (nBlocks == 0)
        return;
    const unsigned int nThreads = std::max(
        1u, std::min(requestedThreads, static_cast<unsigned int>(nBlocks)));
    if (nThreads == 1u) {
        for (size_t block = 0; block < nBlocks; ++block)
            function(block);
        return;
    }

    std::atomic<size_t> nextBlock(0);
    auto worker = [&]() {
        for (;;) {
            const size_t block = nextBlock.fetch_add(1, std::memory_order_relaxed);
            if (block >= nBlocks)
                break;
            function(block);
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(nThreads - 1u);
    for (unsigned int thread = 1; thread < nThreads; ++thread)
        workers.emplace_back(worker);
    worker();
    for (std::thread& thread : workers)
        thread.join();
}

bool writeAllAt(int fd, const char* data, size_t length, uint64_t offset)
{
    while (length != 0u) {
        const ssize_t written = ::pwrite(
            fd, data, length, static_cast<off_t>(offset));
        if (written < 0 && errno == EINTR)
            continue;
        if (written == 0) {
            errno = EIO;
            return false;
        }
        if (written < 0)
            return false;
        data += written;
        length -= static_cast<size_t>(written);
        offset += static_cast<uint64_t>(written);
    }
    return true;
}

int writeMatrix(const std::string& path,
                const std::vector<MexWriter::Triplet>& triplets,
                size_t nBarcodes,
                size_t nFeatures,
                unsigned int matrixThreads)
{
    char dimensions[96];
    const int dimensionsLength = std::snprintf(
        dimensions, sizeof(dimensions), "%u %u %llu\n",
        static_cast<uint32_t>(nFeatures),
        static_cast<uint32_t>(nBarcodes),
        static_cast<unsigned long long>(triplets.size()));
    if (dimensionsLength < 0
        || static_cast<size_t>(dimensionsLength) >= sizeof(dimensions)) {
        std::fprintf(stderr, "[MexWriter] ERROR: failed to format matrix header\n");
        return -1;
    }
    const std::string header =
        "%%MatrixMarket matrix coordinate integer general\n%\n"
        + std::string(dimensions, static_cast<size_t>(dimensionsLength));

    static const size_t kEntriesPerBlock = 1u << 18;
    const size_t nBlocks = triplets.empty()
        ? 0u : 1u + (triplets.size() - 1u) / kEntriesPerBlock;
    std::vector<MatrixBlock> blocks(nBlocks);
    for (size_t block = 0; block < nBlocks; ++block) {
        blocks[block].begin = block * kEntriesPerBlock;
        blocks[block].end = std::min(
            triplets.size(), blocks[block].begin + kEntriesPerBlock);
        blocks[block].fileOffset = 0;
        blocks[block].textBytes = 0;
    }

    parallelBlocks(nBlocks, matrixThreads, [&](size_t block) {
        uint64_t bytes = 0;
        for (size_t index = blocks[block].begin;
             index < blocks[block].end; ++index) {
            if (validTriplet(triplets[index], nBarcodes, nFeatures))
                bytes += tripletTextBytes(triplets[index]);
        }
        blocks[block].textBytes = bytes;
    });

    uint64_t fileBytes = header.size();
    for (MatrixBlock& block : blocks) {
        block.fileOffset = fileBytes;
        if (block.textBytes > std::numeric_limits<uint64_t>::max() - fileBytes) {
            std::fprintf(stderr, "[MexWriter] ERROR: matrix output size overflow\n");
            return -1;
        }
        fileBytes += block.textBytes;
    }
    if (fileBytes > static_cast<uint64_t>(std::numeric_limits<off_t>::max())) {
        std::fprintf(stderr, "[MexWriter] ERROR: matrix output exceeds off_t range\n");
        return -1;
    }

    const int fd = ::open(path.c_str(), O_CREAT | O_TRUNC | O_RDWR, 0666);
    if (fd < 0) {
        std::fprintf(stderr, "[MexWriter] ERROR: failed to open %s: %s\n",
                     path.c_str(), std::strerror(errno));
        return -1;
    }
    if (::ftruncate(fd, static_cast<off_t>(fileBytes)) != 0
        || !writeAllAt(fd, header.data(), header.size(), 0u)) {
        const int savedErrno = errno;
        std::fprintf(stderr, "[MexWriter] ERROR: failed to size/write %s: %s\n",
                     path.c_str(), std::strerror(savedErrno));
        ::close(fd);
        return -1;
    }

    long pageSizeValue = ::sysconf(_SC_PAGESIZE);
    if (pageSizeValue <= 0)
        pageSizeValue = 4096;
    const uint64_t pageSize = static_cast<uint64_t>(pageSizeValue);
    std::atomic<int> firstError(0);
    parallelBlocks(nBlocks, matrixThreads, [&](size_t blockIndex) {
        if (firstError.load(std::memory_order_relaxed) != 0)
            return;
        const MatrixBlock& block = blocks[blockIndex];
        if (block.textBytes == 0u)
            return;
        const uint64_t mapOffset = (block.fileOffset / pageSize) * pageSize;
        const size_t prefix = static_cast<size_t>(block.fileOffset - mapOffset);
        if (block.textBytes > std::numeric_limits<size_t>::max() - prefix) {
            int expected = 0;
            firstError.compare_exchange_strong(expected, EOVERFLOW);
            return;
        }
        const size_t mapLength = prefix + static_cast<size_t>(block.textBytes);
        void* mapping = ::mmap(nullptr, mapLength, PROT_READ | PROT_WRITE,
                               MAP_SHARED, fd, static_cast<off_t>(mapOffset));
        if (mapping == MAP_FAILED) {
            const int mapError = errno;
            int expected = 0;
            firstError.compare_exchange_strong(expected, mapError);
            return;
        }

        char* output = static_cast<char*>(mapping) + prefix;
        char* const expectedEnd = output + block.textBytes;
        for (size_t index = block.begin; index < block.end; ++index) {
            const MexWriter::Triplet& triplet = triplets[index];
            if (!validTriplet(triplet, nBarcodes, nFeatures))
                continue;
            output = appendUnsigned(output, triplet.gene_idx + 1u);
            *output++ = ' ';
            output = appendUnsigned(output, triplet.cell_idx + 1u);
            *output++ = ' ';
            output = appendUnsigned(output, triplet.count);
            *output++ = '\n';
        }
        if (output != expectedEnd) {
            int expected = 0;
            firstError.compare_exchange_strong(expected, EIO);
        }
        if (::munmap(mapping, mapLength) != 0) {
            const int unmapError = errno;
            int expected = 0;
            firstError.compare_exchange_strong(expected, unmapError);
        }
    });

    const int writeError = firstError.load(std::memory_order_relaxed);
    const int closeResult = ::close(fd);
    if (writeError != 0 || closeResult != 0) {
        const int savedErrno = writeError != 0 ? writeError : errno;
        std::fprintf(stderr, "[MexWriter] ERROR: failed to write %s: %s\n",
                     path.c_str(), std::strerror(savedErrno));
        return -1;
    }
    return 0;
}

} // namespace

namespace MexWriter {

int writeMex(const std::string& outputPrefix,
             const std::vector<std::string>& barcodes,
             const std::vector<Feature>& features,
             const std::vector<Triplet>& triplets,
             int cb_len)
{
    return writeMex(outputPrefix, barcodes, features, triplets, cb_len, 1u);
}

int writeMex(const std::string& outputPrefix,
             const std::vector<std::string>& barcodes,
             const std::vector<Feature>& features,
             const std::vector<Triplet>& triplets,
             int cb_len,
             unsigned int matrix_threads)
{
    if (barcodes.empty() || features.empty()) {
        std::fprintf(stderr, "[MexWriter] ERROR: empty input (barcodes=%zu, features=%zu)\n",
                     barcodes.size(), features.size());
        return -1; // No data to write
    }
    
    // Normalize barcodes if cb_len > 0 (truncate to cb_len characters)
    std::vector<std::string> truncatedBarcodes;
    const std::vector<std::string> *outputBarcodesPtr = &barcodes;
    if (cb_len > 0) {
        truncatedBarcodes.reserve(barcodes.size());
        std::unordered_set<std::string> seen;
        seen.reserve(barcodes.size());
        
        for (size_t i = 0; i < barcodes.size(); i++) {
            const std::string& bc = barcodes[i];
            std::string truncated;
            
            if (static_cast<int>(bc.size()) <= cb_len) {
                truncated = bc;
            } else {
                truncated = bc.substr(0, static_cast<size_t>(cb_len));
            }
            
            if (seen.find(truncated) != seen.end()) {
                std::fprintf(stderr, "[MexWriter] ERROR: duplicate barcode after truncation to %d bp\n", cb_len);
                std::fprintf(stderr, "  barcode[%zu] = '%s' -> '%s' (already exists)\n", 
                            i, bc.c_str(), truncated.c_str());
                std::fprintf(stderr, "  This should not happen in per-sample MEX output.\n");
                std::fprintf(stderr, "  SOLUTION: Check for barcode collisions or disable truncation:\n");
                std::fprintf(stderr, "    - External tool: --keep-cb-tag\n");
                std::fprintf(stderr, "    - STAR inline:   --soloFlexKeepCBTag yes\n");
                return -1;
            }
            
            seen.insert(truncated);
            truncatedBarcodes.push_back(truncated);
        }
        outputBarcodesPtr = &truncatedBarcodes;
    }
    const std::vector<std::string> &outputBarcodes = *outputBarcodesPtr;
    
    // Create output file paths
    bool dirMode = (!outputPrefix.empty() && outputPrefix.back() == '/');
    std::string base = outputPrefix;
    if (dirMode && base.size() > 1 && base.back() == '/') {
        // leave trailing slash so we emit files inside the directory
    }
    std::string mtx_path = dirMode ? (base + "matrix.mtx") : (outputPrefix + "_matrix.mtx");
    std::string barcodes_path = dirMode ? (base + "barcodes.tsv") : (outputPrefix + "_barcodes.tsv");
    std::string features_path = dirMode ? (base + "features.tsv") : (outputPrefix + "_features.tsv");
    
    FILE *barcodes_fp = fopen(barcodes_path.c_str(), "w");
    FILE *features_fp = fopen(features_path.c_str(), "w");
    
    if (!barcodes_fp || !features_fp) {
        std::fprintf(stderr, "[MexWriter] ERROR: failed to open output files\n");
        if (!barcodes_fp) std::fprintf(stderr, "  failed: %s\n", barcodes_path.c_str());
        if (!features_fp) std::fprintf(stderr, "  failed: %s\n", features_path.c_str());
        if (barcodes_fp) fclose(barcodes_fp);
        if (features_fp) fclose(features_fp);
        return -1;
    }

    if (writeMatrix(mtx_path, triplets, outputBarcodes.size(), features.size(),
                    std::max(1u, matrix_threads)) != 0) {
        fclose(barcodes_fp);
        fclose(features_fp);
        return -1;
    }
    
    // Write barcodes (ordered by column) - using normalized barcodes
    for (const auto& bc : outputBarcodes) {
        fprintf(barcodes_fp, "%s\n", bc.c_str());
    }
    
    // Write features (ordered by row)
    for (const auto& feat : features) {
        fprintf(features_fp, "%s\t%s\t%s\n", 
               feat.id.c_str(), 
               feat.name.c_str(),
               feat.featureType.c_str());
    }
    
    const int barcodeClose = fclose(barcodes_fp);
    const int featureClose = fclose(features_fp);
    return barcodeClose == 0 && featureClose == 0 ? 0 : -1;
}

int writeMex(const std::string& outputPrefix,
             const std::vector<std::string>& barcodes,
             const std::vector<std::string>& featureIds,
             const std::vector<Triplet>& triplets,
             int cb_len)
{
    return writeMex(outputPrefix, barcodes, featureIds, triplets, cb_len, 1u);
}

int writeMex(const std::string& outputPrefix,
             const std::vector<std::string>& barcodes,
             const std::vector<std::string>& featureIds,
             const std::vector<Triplet>& triplets,
             int cb_len,
             unsigned int matrix_threads)
{
    // Convert simple feature IDs to Feature structs
    std::vector<Feature> features;
    features.reserve(featureIds.size());
    for (const auto& id : featureIds) {
        features.emplace_back(id, id, "Gene Expression");
    }
    
    return writeMex(
        outputPrefix, barcodes, features, triplets, cb_len, matrix_threads);
}

} // namespace MexWriter

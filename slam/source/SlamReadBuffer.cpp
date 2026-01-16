#include "SlamReadBuffer.h"

SlamReadBuffer::SlamReadBuffer(uint64_t maxReads)
    : maxReads_(maxReads)
{
    reads_.reserve(std::min(maxReads_, static_cast<uint64_t>(100000)));
}

bool SlamReadBuffer::addRead(SlamBufferedRead&& read) {
    if (reads_.size() >= maxReads_) {
        return false;
    }
    reads_.push_back(std::move(read));
    return true;
}

void SlamReadBuffer::clear() {
    reads_.clear();
    reads_.shrink_to_fit();
}

uint64_t SlamReadBuffer::memoryBytes() const {
    uint64_t total = sizeof(*this);
    total += reads_.capacity() * sizeof(SlamBufferedRead);
    
    for (const auto& read : reads_) {
        total += read.readName.capacity();
        total += read.geneIds.capacity() * sizeof(uint32_t);
        total += read.positions.capacity() * sizeof(SlamBufferedPosition);
        total += read.mismatchGenomicPositions.capacity() * sizeof(uint32_t);
    }
    
    return total;
}

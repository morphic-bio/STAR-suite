#ifndef SAMTOOLS_SORTER_H
#define SAMTOOLS_SORTER_H

#include "IncludeDefine.h"
#include "Parameters.h"
#include <vector>
#include <string>
#include <pthread.h>
#include <cstdint>
#include <climits>
#include <cassert>
#include <cstring>
#include <utility>

// Simplified coordinate sorting: legacy-style coordinate order (tid<<32|pos, unmapped→INT32_MAX) then readId
// No QNAME tie-break, no SortKey storage - recompute coord from bamData in comparator

// Structure to hold BAM record data with minimal metadata
struct BAMRecord {
    char* data;
    uint32_t size;
    bool hasY;
    uint32_t readId;  // readId extracted from iRead >> 32
    
    BAMRecord() : data(nullptr), size(0), hasY(false), readId(0) {}
    
    // Move constructor
    BAMRecord(BAMRecord&& other) noexcept 
        : data(other.data), size(other.size), hasY(other.hasY), readId(other.readId) {
        other.data = nullptr;
        other.size = 0;
        other.readId = 0;
    }
    
    // Move assignment
    BAMRecord& operator=(BAMRecord&& other) noexcept {
        if (this != &other) {
            if (data) delete[] data;
            data = other.data;
            size = other.size;
            hasY = other.hasY;
            readId = other.readId;
            other.data = nullptr;
            other.size = 0;
            other.readId = 0;
        }
        return *this;
    }
    
    // Delete copy constructor and assignment to prevent accidental copies
    BAMRecord(const BAMRecord& other) = delete;
    BAMRecord& operator=(const BAMRecord& other) = delete;
    
    ~BAMRecord() {
        if (data) delete[] data;
    }
};

// Helper to compute coordinate key from raw BAM buffer (tid<<32|pos, unmapped→INT32_MAX)
namespace SamtoolsSorterHelpers {
    inline uint64_t computeCoordKey(const char* bamData) {
        const uint32_t* bam32 = reinterpret_cast<const uint32_t*>(bamData);
        int32_t tid = static_cast<int32_t>(bam32[1]);
        int32_t pos = static_cast<int32_t>(bam32[2]);
        
        // Unmapped reads: set to INT32_MAX so they sort last
        uint64_t tid64 = (tid == -1) ? static_cast<uint64_t>(INT32_MAX) : static_cast<uint64_t>(tid);
        uint64_t pos64 = (pos == -1) ? static_cast<uint64_t>(INT32_MAX) : static_cast<uint64_t>(pos);
        
        return (tid64 << 32) | pos64;
    }
}

// Comparator for sorting BAM records: coord (tid<<32|pos) then readId
struct BAMRecordComparator {
    bool operator()(const BAMRecord& a, const BAMRecord& b) const {
        // Compute coord keys from raw BAM buffers
        uint64_t coordA = SamtoolsSorterHelpers::computeCoordKey(a.data);
        uint64_t coordB = SamtoolsSorterHelpers::computeCoordKey(b.data);
        
        if (coordA != coordB) {
            return coordA < coordB;
        }
        
        // Tie-break by readId
        return a.readId < b.readId;
    }
};

class SamtoolsSorter {
public:
    // Initialize with memory limit and thread count
    SamtoolsSorter(uint64_t maxRAM, int nThreads, const string& tmpDir, Parameters& P);
    
    // Thread-safe: called from coordOneAlign replacement path
    // Takes raw BAM record bytes (same format as coordOneAlign input)
    // readId: uint32 readId extracted from iRead >> 32
    // hasY: separate Y-chromosome flag
    void addRecord(const char* bamData, uint32_t bamSize, uint32_t readId, bool hasY);
    
    // Called after all threads complete; performs sort + merge
    void finalize();
    
    // Iterate sorted records for output (post-sort streaming via k-way merge)
    // Returns false when no more records
    // Caller does NOT own the returned data pointer - it's valid until next call
    // readId: returns uint32 readId for tag injection
    bool nextRecord(const char** bamData, uint32_t* bamSize, uint32_t* readId, bool* hasY);
    
    // Cleanup temp files
    ~SamtoolsSorter();

private:
    uint64_t maxRAM_;
    int nThreads_;
    string tmpDir_;
    Parameters& P_;
    
    // Thread-safe buffer
    pthread_mutex_t bufferMutex_;
    std::vector<BAMRecord> records_;
    uint64_t currentRAM_;
    
    // Spill files for when memory limit is exceeded
    std::vector<string> spillFiles_;
    pthread_mutex_t spillFileCounterMutex_;  // Protects spillFileCounter_ from concurrent access
    int spillFileCounter_;
    
    // Forward declaration - SpillFileReader defined in .cpp file
    struct SpillFileReader;
    std::vector<SpillFileReader*> spillReaders_;
    bool finalized_;
    BAMRecord returnedRecord_;  // Spill record buffer returned to caller (valid until next call)
    
    // Min-heap for k-way merge: coord (tid<<32|pos) then readId
    // source_id: -1 = in-memory records, >=0 = spill file index
    struct HeapEntry {
        uint64_t coordKey;  // tid<<32|pos (computed from bamPtr)
        uint32_t readId;    // readId for tie-breaking
        int sourceId;       // -1 = in-memory, >=0 spill file index
        size_t recordIdx;   // index into records_ when sourceId == -1
        const char* bamPtr; // pointer to current BAM record (for recomputing coord)
        BAMRecord record;   // Owned copy of BAM record (for spill files, ensures data persists after readNext)
        
        HeapEntry() : coordKey(0), readId(0), sourceId(0), recordIdx(0), bamPtr(nullptr) {}
        
        // Move constructor
        HeapEntry(HeapEntry&& other) noexcept
            : coordKey(other.coordKey), readId(other.readId), sourceId(other.sourceId),
              recordIdx(other.recordIdx), bamPtr(other.bamPtr), record(std::move(other.record)) {
            other.bamPtr = nullptr;
        }
        
        // Move assignment
        HeapEntry& operator=(HeapEntry&& other) noexcept {
            if (this != &other) {
                coordKey = other.coordKey;
                readId = other.readId;
                sourceId = other.sourceId;
                recordIdx = other.recordIdx;
                bamPtr = other.bamPtr;
                record = std::move(other.record);
                other.bamPtr = nullptr;
            }
            return *this;
        }
        
        // Copy constructor (needed for heap operations)
        HeapEntry(const HeapEntry& other) 
            : coordKey(other.coordKey), readId(other.readId), sourceId(other.sourceId),
              recordIdx(other.recordIdx), bamPtr(other.bamPtr) {
            // Deep copy record data (only if it exists - in-memory records have nullptr)
            record.size = other.record.size;
            record.hasY = other.record.hasY;
            record.readId = other.record.readId;
            if (other.record.data != nullptr && other.record.size > 0) {
                record.data = new char[other.record.size];
                memcpy(record.data, other.record.data, other.record.size);
            } else {
                record.data = nullptr;
            }
        }
        
        // Copy assignment
        HeapEntry& operator=(const HeapEntry& other) {
            if (this != &other) {
                coordKey = other.coordKey;
                readId = other.readId;
                sourceId = other.sourceId;
                recordIdx = other.recordIdx;
                bamPtr = other.bamPtr;
                // Deep copy record data (only if it exists - in-memory records have nullptr)
                if (record.data) delete[] record.data;
                record.size = other.record.size;
                record.hasY = other.record.hasY;
                record.readId = other.record.readId;
                if (other.record.data != nullptr && other.record.size > 0) {
                    record.data = new char[other.record.size];
                    memcpy(record.data, other.record.data, other.record.size);
                } else {
                    record.data = nullptr;
                }
            }
            return *this;
        }
    };
    
    // Comparator for min-heap: coord then readId
    struct HeapLess {
        bool operator()(const HeapEntry& a, const HeapEntry& b) const {
            if (a.coordKey != b.coordKey) {
                return a.coordKey > b.coordKey;  // min-heap: smaller coord has higher priority
            }
            // Tie-break by readId
            return a.readId > b.readId;  // min-heap: smaller readId has higher priority
        }
    };
    
    std::vector<HeapEntry> mergeHeap_;
    HeapLess heapLess_;  // Comparator instance
    
    // Helper functions
    void sortAndSpill();
    void initializeKWayMerge();
    void cleanupSpillFiles();
    void clearReturnedRecord();
};

#endif // SAMTOOLS_SORTER_H

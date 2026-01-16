#include "SamtoolsSorter.h"
#include "BAMfunctions.h"
#include "ErrorWarning.h"
#include "Parameters.h"
#include SAMTOOLS_BGZF_H
#include SAMTOOLS_SAM_H
#include <algorithm>
#include <fstream>
#include <cstring>
#include <pthread.h>
#include <queue>
#include <climits>

// SpillFileReader implementation for k-way merge
struct SamtoolsSorter::SpillFileReader {
    ifstream stream;
    BAMRecord currentRecord;
    bool hasRecord;
    string filename;
    int sourceId;
    
    SpillFileReader(const string& fname, int id) : hasRecord(false), filename(fname), sourceId(id) {
        stream.open(fname.c_str(), std::ios::binary);
    }
    
    bool readNext() {
        if (!stream.good()) {
            hasRecord = false;
            return false;
        }
        uint32_t size;
        uint8_t hasYFlag;
        uint32_t readId;
        
        stream.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
        if (stream.eof()) {
            hasRecord = false;
            return false;
        }
        
        stream.read(reinterpret_cast<char*>(&hasYFlag), sizeof(uint8_t));
        if (!stream.good()) {
            hasRecord = false;
            return false;
        }
        
        stream.read(reinterpret_cast<char*>(&readId), sizeof(uint32_t));
        if (!stream.good()) {
            hasRecord = false;
            return false;
        }
        
        currentRecord.size = size;
        currentRecord.hasY = (hasYFlag != 0);
        currentRecord.readId = readId;
        if (currentRecord.data) delete[] currentRecord.data;
        currentRecord.data = new char[size];
        stream.read(currentRecord.data, size);
        
        if (stream.good()) {
            hasRecord = true;
            return true;
        } else {
            hasRecord = false;
            return false;
        }
    }
    
    ~SpillFileReader() {
        if (stream.is_open()) stream.close();
    }
};

SamtoolsSorter::SamtoolsSorter(uint64_t maxRAM, int nThreads, const string& tmpDir, Parameters& P)
    : maxRAM_(maxRAM), nThreads_(nThreads), tmpDir_(tmpDir), P_(P),
      currentRAM_(0), spillFileCounter_(0), finalized_(false)
{
    records_.reserve(100000); // Pre-allocate space
    pthread_mutex_init(&bufferMutex_, nullptr);
    pthread_mutex_init(&spillFileCounterMutex_, nullptr);
}

SamtoolsSorter::~SamtoolsSorter() {
    cleanupSpillFiles();
    pthread_mutex_destroy(&bufferMutex_);
    pthread_mutex_destroy(&spillFileCounterMutex_);
}


void SamtoolsSorter::addRecord(const char* bamData, uint32_t bamSize, uint32_t readId, bool hasY) {
    if (bamSize == 0) return;
    
    BAMRecord record;
    record.size = bamSize;
    record.hasY = hasY;
    record.readId = readId;
    record.data = new char[bamSize];
    memcpy(record.data, bamData, bamSize);
    
    // Quick check without lock for common case
    bool needSpill = false;
    {
        pthread_mutex_lock(&bufferMutex_);
        records_.push_back(std::move(record));
        // Estimate RAM usage: record data + overhead
        currentRAM_ += bamSize + sizeof(BAMRecord) + 64; // 64 bytes overhead per record
        
        // Check if we need to spill to disk
        if (currentRAM_ > maxRAM_ && maxRAM_ > 0 && !records_.empty()) {
            needSpill = true;
        }
        pthread_mutex_unlock(&bufferMutex_);
    }
    
    // Spill outside the lock to avoid blocking writers
    if (needSpill) {
        sortAndSpill();
    }
}

void SamtoolsSorter::sortAndSpill() {
    // Lock only to extract records, then release for sorting/writing
    std::vector<BAMRecord> recordsToSpill;
    {
        pthread_mutex_lock(&bufferMutex_);
        if (records_.empty()) {
            pthread_mutex_unlock(&bufferMutex_);
            return;
        }
        // Move all records out (no copy due to move semantics)
        recordsToSpill = std::move(records_);
        records_.clear();
        currentRAM_ = 0;
        pthread_mutex_unlock(&bufferMutex_);
    }
    
    // Sort records by SortKey + QNAME (outside lock)
    std::sort(recordsToSpill.begin(), recordsToSpill.end(), BAMRecordComparator());
    
    // Write sorted chunk to temp file with SortKey fields serialized explicitly
    // Protect spillFileCounter_ from concurrent access
    int spillFileNum;
    {
        pthread_mutex_lock(&spillFileCounterMutex_);
        spillFileNum = spillFileCounter_++;
        pthread_mutex_unlock(&spillFileCounterMutex_);
    }
    string spillFile = tmpDir_ + "/samtools_sort_spill_" + to_string(spillFileNum) + ".dat";
    
    ofstream spillStream(spillFile.c_str(), std::ios::binary);
    if (!spillStream) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not open spill file: " << spillFile << "\n";
        errOut << "SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
        exitWithError(errOut.str(), std::cerr, P_.inOut->logMain, EXIT_CODE_PARAMETER, P_);
    }
    
    // Write records to spill file with minimal metadata
    // Format: [bamSize:uint32][hasY:uint8][readId:uint32][bamData:bytes]
    for (auto& record : recordsToSpill) {
        uint32_t size = record.size;
        uint8_t hasYFlag = record.hasY ? 1 : 0;
        spillStream.write(reinterpret_cast<const char*>(&size), sizeof(uint32_t));
        spillStream.write(reinterpret_cast<const char*>(&hasYFlag), sizeof(uint8_t));
        spillStream.write(reinterpret_cast<const char*>(&record.readId), sizeof(uint32_t));
        spillStream.write(record.data, record.size);
    }
    
    spillStream.close();
    
    // Store spill file name (will be used for k-way merge)
    // Use spillFileCounterMutex_ since we're already protecting counter access
    {
        pthread_mutex_lock(&spillFileCounterMutex_);
        spillFiles_.push_back(spillFile);
        pthread_mutex_unlock(&spillFileCounterMutex_);
    }
}

void SamtoolsSorter::finalize() {
    if (finalized_) return;
    
    pthread_mutex_lock(&bufferMutex_);
    
    // Sort remaining in-memory records
    if (!records_.empty()) {
        std::sort(records_.begin(), records_.end(), BAMRecordComparator());
    }
    
    // Initialize k-way merge with spill files
    initializeKWayMerge();
    
    finalized_ = true;
    pthread_mutex_unlock(&bufferMutex_);
}

void SamtoolsSorter::clearReturnedRecord() {
    if (returnedRecord_.data) {
        delete[] returnedRecord_.data;
        returnedRecord_.data = nullptr;
        returnedRecord_.size = 0;
        returnedRecord_.hasY = false;
        returnedRecord_.readId = 0;
    }
}

void SamtoolsSorter::initializeKWayMerge() {
    mergeHeap_.clear();
    
    // Add in-memory records to heap (if any)
    if (!records_.empty()) {
        HeapEntry entry;
        entry.coordKey = SamtoolsSorterHelpers::computeCoordKey(records_[0].data);
        entry.readId = records_[0].readId;
        entry.sourceId = -1;  // -1 indicates in-memory source
        entry.recordIdx = 0;   // First record in sorted vector
        entry.bamPtr = records_[0].data;  // Pointer to BAM buffer for recomputing coord
        // Initialize record to empty for in-memory (we use recordIdx to index into records_)
        entry.record.data = nullptr;
        entry.record.size = 0;
        entry.record.hasY = false;
        entry.record.readId = 0;
        mergeHeap_.push_back(std::move(entry));
    }
    
    // Open all spill files for reading and add first record from each to heap
    for (size_t i = 0; i < spillFiles_.size(); i++) {
        SpillFileReader* reader = new SpillFileReader(spillFiles_[i], static_cast<int>(i));
        if (reader->readNext()) {
            HeapEntry entry;
            entry.coordKey = SamtoolsSorterHelpers::computeCoordKey(reader->currentRecord.data);
            entry.readId = reader->currentRecord.readId;
            entry.sourceId = reader->sourceId;
            entry.recordIdx = 0;  // Unused for spill files
            entry.bamPtr = reader->currentRecord.data;  // Pointer to BAM buffer for recomputing coord
            // CRITICAL: Copy record data into heap entry so it persists after readNext() advances
            entry.record.size = reader->currentRecord.size;
            entry.record.hasY = reader->currentRecord.hasY;
            entry.record.readId = reader->currentRecord.readId;
            entry.record.data = new char[reader->currentRecord.size];
            memcpy(entry.record.data, reader->currentRecord.data, reader->currentRecord.size);
            mergeHeap_.push_back(std::move(entry));
            spillReaders_.push_back(reader);
        } else {
            delete reader;
        }
    }
    
    // Build min-heap using explicit comparator: coord then readId
    std::make_heap(mergeHeap_.begin(), mergeHeap_.end(), heapLess_);
}

bool SamtoolsSorter::nextRecord(const char** bamData, uint32_t* bamSize, uint32_t* readId, bool* hasY) {
    if (!finalized_) {
        finalize();
    }
    
    // Free previous spill-return buffer (contract: valid until next call)
    clearReturnedRecord();
    
    // K-way merge using min-heap: coord (tid<<32|pos) then readId
    while (!mergeHeap_.empty()) {
        // Get smallest entry from heap (coord then readId)
        std::pop_heap(mergeHeap_.begin(), mergeHeap_.end(), heapLess_);
        HeapEntry top = std::move(mergeHeap_.back());  // Move instead of copy to avoid double-free
        mergeHeap_.pop_back();
        
        // Get the record for this entry
        BAMRecord* record = nullptr;
        if (top.sourceId == -1) {
            // In-memory source - use recordIdx to get the correct record
            if (top.recordIdx >= records_.size()) {
                continue;  // Invalid index, try next heap entry
            }
            record = &records_[top.recordIdx];
        } else {
            // Spill file source - move record into persistent buffer for this call
            returnedRecord_ = std::move(top.record);
            record = &returnedRecord_;
        }
        
        if (record == nullptr || record->data == nullptr) {
            continue;
        }
        
        // Advance the source for the record we're returning
        if (top.sourceId == -1) {
            // Advance to next in-memory record
            size_t nextIdx = top.recordIdx + 1;
            
            // Push next in-memory record if available
            if (nextIdx < records_.size()) {
                HeapEntry nextEntry;
                nextEntry.coordKey = SamtoolsSorterHelpers::computeCoordKey(records_[nextIdx].data);
                nextEntry.readId = records_[nextIdx].readId;
                nextEntry.sourceId = -1;
                nextEntry.recordIdx = nextIdx;
                nextEntry.bamPtr = records_[nextIdx].data;
                // Initialize record to empty for in-memory (we use recordIdx to index into records_)
                nextEntry.record.data = nullptr;
                nextEntry.record.size = 0;
                nextEntry.record.hasY = false;
                nextEntry.record.readId = 0;
                mergeHeap_.push_back(std::move(nextEntry));
                std::push_heap(mergeHeap_.begin(), mergeHeap_.end(), heapLess_);
            }
        } else {
            // Read next record from this spill file
            SpillFileReader* reader = spillReaders_[top.sourceId];
            if (reader->readNext()) {
                // Push next record from this source
                HeapEntry nextEntry;
                nextEntry.coordKey = SamtoolsSorterHelpers::computeCoordKey(reader->currentRecord.data);
                nextEntry.readId = reader->currentRecord.readId;
                nextEntry.sourceId = reader->sourceId;
                nextEntry.recordIdx = 0;  // Unused for spill files
                nextEntry.bamPtr = reader->currentRecord.data;
                // CRITICAL: Copy record data into heap entry so it persists after readNext() advances
                nextEntry.record.size = reader->currentRecord.size;
                nextEntry.record.hasY = reader->currentRecord.hasY;
                nextEntry.record.readId = reader->currentRecord.readId;
                nextEntry.record.data = new char[reader->currentRecord.size];
                memcpy(nextEntry.record.data, reader->currentRecord.data, reader->currentRecord.size);
                mergeHeap_.push_back(std::move(nextEntry));
                std::push_heap(mergeHeap_.begin(), mergeHeap_.end(), heapLess_);
            } else {
                // File exhausted, will be cleaned up later
            }
        }
        
        // Return the record (spill records are buffered until the next call)
        *bamData = record->data;
        *bamSize = record->size;
        *readId = record->readId;
        *hasY = record->hasY;
        
        return true;
    }
    
    return false;
}

void SamtoolsSorter::cleanupSpillFiles() {
    // Close and delete all spill file readers
    for (auto* reader : spillReaders_) {
        delete reader;
    }
    spillReaders_.clear();
    
    // Delete spill files
    for (const string& spillFile : spillFiles_) {
        remove(spillFile.c_str());
    }
    spillFiles_.clear();
}

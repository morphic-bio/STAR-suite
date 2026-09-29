#include "input/FastxMateReaders.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>

namespace star {
namespace input {

namespace {

constexpr int kNoEndLog = 1000;  // end of input caused by a mate other than mate 1

uint64_t nowNs() {
    return static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());
}

// Same as fastqHeaderExtraFromCurrentLine in ReadAlignChunk_processChunks.cpp.
std::string headerExtra(std::istream& in) {
    std::string extra;
    std::getline(in, extra);
    const size_t firstNonSpace = extra.find_first_not_of(" \t");
    if (firstNonSpace == std::string::npos) {
        extra.clear();
    } else if (firstNonSpace > 0) {
        extra.erase(0, firstNonSpace);
    }
    while (!extra.empty() && static_cast<unsigned char>(extra.back()) < 33) {
        extra.pop_back();
    }
    return extra;
}

// Same as illuminaFilterFlagFromHeaderExtra in ReadAlignChunk_processChunks.cpp.
char illuminaFilterFlag(const std::string& extra) {
    std::string field2;
    std::istringstream extraStream(extra);
    extraStream >> field2;
    if (field2.length() >= 4 && field2[1] == ':' && field2[2] == 'Y' && field2[3] == ':') {
        return 'Y';
    }
    return 'N';
}

const char* formatName(char format) {
    return format == '>' ? "FASTA" : "FASTQ";
}

}  // namespace

// ---------------------------------------------------------------- InputBuf

FastxMateReader::InputBuf::InputBuf(FastxMateReader* owner, std::streambuf* source)
    : owner_(owner), source_(source), buffer_(kFastxMateReadBufferBytes) {
    setg(buffer_.data(), buffer_.data(), buffer_.data());
}

// Refills from the mate's stream. The permit is given back for the wait, so a
// reader never holds one while its producer (which may itself need a permit)
// is behind.
FastxMateReader::InputBuf::int_type FastxMateReader::InputBuf::underflow() {
    if (gptr() < egptr()) {
        return traits_type::to_int_type(*gptr());
    }
    const bool held = owner_->holding_;
    if (held) {
        owner_->releasePermit();
    }
    const uint64_t start = nowNs();
    const std::streamsize got =
        source_->sgetn(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
    owner_->stats_.inputWaitNs += nowNs() - start;
    if (held) {
        owner_->acquirePermit();
    }
    if (got <= 0) {
        return traits_type::eof();
    }
    owner_->stats_.bytes += static_cast<uint64_t>(got);
    setg(buffer_.data(), buffer_.data(), buffer_.data() + got);
    return traits_type::to_int_type(*gptr());
}

// --------------------------------------------------------- FastxMateReader

FastxMateReader::FastxMateReader(uint32_t mate, std::istream* source, int initialLane,
                                 const BgzfWorkPermitHooks& hooks,
                                 const FastxMateLimits& limits)
    : mate_(mate), source_(source), lane_(initialLane), hooks_(hooks), limits_(limits) {}

FastxMateReader::~FastxMateReader() {
    requestStop();
    join();
}

void FastxMateReader::start() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (started_ || stop_.load()) {
        return;
    }
    const long long longest = std::max(limits_.nameSeqLineMax, limits_.seqLineMax);
    line_.assign(static_cast<size_t>(longest) + 2, '\0');
    owned_.reserve(kFastxMateBatchesInFlight);
    for (size_t i = 0; i < kFastxMateBatchesInFlight; ++i) {
        owned_.emplace_back(new FastxMateBatch());
        owned_.back()->records.reserve(kFastxMateBatchRecords);
        free_.push_back(owned_.back().get());
    }
    started_ = true;
    thread_ = std::thread(&FastxMateReader::run, this);
}

void FastxMateReader::requestStop() {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        stop_.store(true);
        if (!started_) {
            finished_ = true;
        }
    }
    freeCv_.notify_all();
    readyCv_.notify_all();
}

void FastxMateReader::join() {
    if (thread_.joinable()) {
        thread_.join();
    }
}

FastxMateReader::Stats FastxMateReader::stats() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return stats_;
}

void FastxMateReader::observeLocked() {
    if (hooks_.observe == nullptr) {
        return;
    }
    const uint64_t ready = ready_.size();
    hooks_.observe(hooks_.context, this, ready, ready + (filling_ ? 1 : 0),
                   kFastxMateBatchesInFlight, 1, consumerWaiting_ ? 1 : 0, finished_ ? 0 : 1);
}

FastxMateBatch* FastxMateReader::takeFree() {
    std::unique_lock<std::mutex> lock(mutex_);
    if (free_.empty() && !stop_.load()) {
        const uint64_t start = nowNs();
        freeCv_.wait(lock, [this] { return stop_.load() || !free_.empty(); });
        stats_.freeWaitNs += nowNs() - start;
    }
    if (stop_.load()) {
        return nullptr;
    }
    FastxMateBatch* batch = free_.back();
    free_.pop_back();
    filling_ = true;
    observeLocked();
    return batch;
}

void FastxMateReader::pushReady(FastxMateBatch* batch) {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        ready_.push_back(batch);
        filling_ = false;
        observeLocked();
    }
    readyCv_.notify_one();
}

FastxMateBatch* FastxMateReader::popReady() {
    std::unique_lock<std::mutex> lock(mutex_);
    if (ready_.empty() && !finished_) {
        consumerWaiting_ = true;
        observeLocked();
        const uint64_t start = nowNs();
        readyCv_.wait(lock, [this] { return !ready_.empty() || finished_; });
        stats_.consumerWaitNs += nowNs() - start;
        consumerWaiting_ = false;
    }
    if (ready_.empty()) {
        observeLocked();
        return nullptr;
    }
    FastxMateBatch* batch = ready_.front();
    ready_.pop_front();
    observeLocked();
    return batch;
}

void FastxMateReader::recycle(FastxMateBatch* batch) {
    if (batch == nullptr) {
        return;
    }
    {
        std::lock_guard<std::mutex> lock(mutex_);
        free_.push_back(batch);
    }
    freeCv_.notify_one();
}

void FastxMateReader::acquirePermit() {
    // A stopping reader parses at most to its next stop check; it does not
    // wait for a permit to do so.
    if (holding_ || !hooks_.enabled() || stop_.load()) {
        return;
    }
    const uint64_t wait = hooks_.acquire(hooks_.context);
    if (wait == UINT64_MAX) {
        return;  // the permit pool is not enabled
    }
    holding_ = true;
    holdWaitNs_ = wait;
    holdStartNs_ = nowNs();
    holdArenaStart_ = current_ != nullptr ? current_->arena.size() : 0;
    ++stats_.permitAcquires;
    stats_.permitWaitNs += wait;
}

void FastxMateReader::releasePermit() {
    if (!holding_) {
        return;
    }
    const size_t arena = current_ != nullptr ? current_->arena.size() : holdArenaStart_;
    const uint64_t bytes = arena > holdArenaStart_ ? arena - holdArenaStart_ : 0;
    hooks_.release(hooks_.context, holdWaitNs_, 1, bytes, nowNs() - holdStartNs_);
    holding_ = false;
}

void FastxMateReader::run() {
    InputBuf buffer(this, source_->rdbuf());
    std::istream in(&buffer);
    for (;;) {
        FastxMateBatch* batch = takeFree();
        if (batch == nullptr) {
            break;
        }
        batch->reset(lane_);
        current_ = batch;
        const uint64_t start = nowNs();
        const uint64_t inputWaitBefore = stats_.inputWaitNs;
        const uint64_t permitWaitBefore = stats_.permitWaitNs;
        acquirePermit();
        bool last = false;
        try {
            last = fillBatch(batch, in);
        } catch (const std::exception& e) {
            batch->end = FastxMateEnd::Error;
            batch->errorText = "EXITING because of FATAL ERROR in the FASTX reader for mate " +
                               std::to_string(mate_ + 1) + ": " + e.what() + "\n";
            last = true;
        }
        releasePermit();
        current_ = nullptr;
        const uint64_t elapsed = nowNs() - start;
        const uint64_t waited = (stats_.inputWaitNs - inputWaitBefore) +
                                (stats_.permitWaitNs - permitWaitBefore);
        stats_.parseNs += elapsed > waited ? elapsed - waited : 0;
        if (stop_.load()) {
            recycle(batch);
            break;
        }
        stats_.records += batch->records.size();
        ++stats_.batches;
        pushReady(batch);
        if (last) {
            break;
        }
    }
    {
        std::lock_guard<std::mutex> lock(mutex_);
        finished_ = true;
        filling_ = false;
        observeLocked();
    }
    readyCv_.notify_all();
}

void FastxMateReader::markLane(FastxMateBatch* batch, int nextLane) {
    batch->end = FastxMateEnd::Lane;
    batch->nextLane = nextLane;
    batch->laneRecords = laneRecords_;
    lane_ = nextLane;
    laneRecords_ = 0;
}

// The single-threaded loop's fastqReadOneLine: getline with its limit, drop
// one trailing byte below 33, end with a newline. For a one-byte result it
// wrote the newline over the preceding newline and added nothing.
bool FastxMateReader::appendLine(FastxMateBatch* batch, std::istream& in,
                                 uint32_t* offset, uint32_t* length) {
    in.getline(line_.data(), static_cast<std::streamsize>(limits_.nameSeqLineMax + 1));
    const std::streamsize got = in.gcount();
    *offset = static_cast<uint32_t>(batch->arena.size());
    if (got <= 0) {
        *length = 0;
        return false;  // end of input inside the record
    }
    if (got == 1) {
        *length = 0;
        return true;
    }
    std::streamsize kept = got;
    if (static_cast<int>(line_[static_cast<size_t>(got - 2)]) < 33) {
        --kept;
    }
    batch->arena.insert(batch->arena.end(), line_.data(), line_.data() + (kept - 1));
    batch->arena.push_back('\n');
    *length = static_cast<uint32_t>(kept);
    return true;
}

bool FastxMateReader::parseFastq(FastxMateBatch* batch, std::istream& in,
                                 FastxMateRecord* record, const std::string& token) {
    if (mate_ == 0) {
        std::string readId = token;
        // removeStringEndControl
        if (!readId.empty() && static_cast<int>(readId.back()) < 33) {
            readId.pop_back();
        }
        record->idOff = batch->append(readId.data(), readId.size());
        record->idLen = static_cast<uint32_t>(readId.size());
    }
    const std::string extra = headerExtra(in);
    if (mate_ == 0) {
        record->filter = illuminaFilterFlag(extra);
    }
    record->extraOff = batch->append(extra.data(), extra.size());
    record->extraLen = static_cast<uint32_t>(extra.size());
    if (!appendLine(batch, in, &record->seqOff, &record->seqLen)) {
        return false;
    }
    in.ignore(static_cast<std::streamsize>(limits_.nameSeqLineMax), '\n');
    return appendLine(batch, in, &record->qualOff, &record->qualLen);
}

// The single-threaded loop's FASTA branch: the header token, the rest of the
// header ignored, then sequence lines joined until the next record start.
void FastxMateReader::parseFasta(FastxMateBatch* batch, std::istream& in,
                                 FastxMateRecord* record, const std::string& token) {
    record->idOff = batch->append(token.data(), token.size());
    record->idLen = static_cast<uint32_t>(token.size());
    in.ignore(static_cast<std::streamsize>(limits_.nameSeqLineMax), '\n');
    const size_t start = batch->arena.size();
    record->seqOff = static_cast<uint32_t>(start);
    int next = in.peek();
    while (next != '@' && next != '>' && next != ' ' && next != '\n' && in.good()) {
        in.getline(line_.data(), static_cast<std::streamsize>(limits_.seqLineMax + 1));
        const std::streamsize got = in.gcount();
        if (got < 2) {
            break;
        }
        std::streamsize kept = got - 1;
        if (static_cast<int>(line_[static_cast<size_t>(kept - 1)]) < 33) {
            --kept;
        }
        batch->arena.insert(batch->arena.end(), line_.data(), line_.data() + kept);
        next = in.peek();
    }
    batch->arena.push_back('\n');
    record->seqLen = static_cast<uint32_t>(batch->arena.size() - start);
}

// Fills one batch. Returns true when this mate's input has ended (or failed).
bool FastxMateReader::fillBatch(FastxMateBatch* batch, std::istream& in) {
    std::string token;
    while (batch->records.size() < kFastxMateBatchRecords) {
        if (stop_.load(std::memory_order_relaxed)) {
            batch->end = FastxMateEnd::Input;
            batch->laneRecords = laneRecords_;
            return true;
        }
        FastxMateRecord record;
        token.clear();
        // The single-threaded loop tested mates 1 and 2 for good() before
        // each read and logged the end of input only when both passed.
        const bool wasGood = in.good();
        if (mate_ == 0) {
            // The single-threaded loop peeks mate 0 to decide what comes next.
            const int next = in.peek();
            if (next == '@' || next == '>') {
                in >> token;
                record.format = static_cast<char>(next);
            } else if (next == ' ' || next == '\n' || !in.good()) {
                batch->end = FastxMateEnd::Input;
                batch->endChar = static_cast<int>(static_cast<char>(next));
                batch->endSilent = !wasGood;
                batch->laneRecords = laneRecords_;
                return true;
            } else {
                std::string word;
                in >> word;
                if (word == "FILE") {
                    int nextLane = 0;
                    in >> nextLane;
                    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    markLane(batch, nextLane);
                    return false;
                }
                std::string rest;
                std::getline(in, rest);
                batch->end = FastxMateEnd::Error;
                batch->errorWord = word;
                batch->errorRest = rest;
                batch->laneRecords = laneRecords_;
                return true;
            }
        } else {
            // Other mates follow mate 0 in the single-threaded loop, which
            // reads their ID token without checking it.
            if (!(in >> token)) {
                batch->end = FastxMateEnd::Input;
                batch->endSilent = !wasGood;
                batch->laneRecords = laneRecords_;
                return true;
            }
            if (token == "FILE") {
                int nextLane = 0;
                in >> nextLane;
                in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                markLane(batch, nextLane);
                return false;
            }
            record.format = token[0] == '>' ? '>' : '@';
        }
        if (record.format == '@') {
            if (!parseFastq(batch, in, &record, token)) {
                // Incomplete last record: this mate's input ends here.
                batch->end = FastxMateEnd::Input;
                batch->laneRecords = laneRecords_;
                return true;
            }
        } else {
            parseFasta(batch, in, &record, token);
        }
        batch->records.push_back(record);
        ++laneRecords_;
    }
    batch->end = FastxMateEnd::None;
    batch->laneRecords = laneRecords_;
    return false;
}

// ---------------------------------------------------- FastxMateReaderGroup

FastxMateReaderGroup::FastxMateReaderGroup(const std::vector<std::istream*>& streams,
                                           int initialLane,
                                           const BgzfWorkPermitHooks& hooks,
                                           const FastxMateLimits& limits) {
    for (size_t m = 0; m < streams.size(); ++m) {
        readers_.emplace_back(new FastxMateReader(static_cast<uint32_t>(m), streams[m],
                                                  initialLane, hooks, limits));
    }
    cursors_.resize(streams.size());
}

FastxMateReaderGroup::~FastxMateReaderGroup() {
    stopAndJoin();
}

void FastxMateReaderGroup::ensureStarted() {
    if (started_ || stopped_) {
        return;
    }
    for (auto& reader : readers_) {
        reader->start();
    }
    started_ = true;
}

void FastxMateReaderGroup::stopAndJoin() {
    stopped_ = true;
    for (auto& reader : readers_) {
        reader->requestStop();
    }
    for (auto& reader : readers_) {
        reader->join();
    }
}

void FastxMateReaderGroup::normalize(uint32_t mate) {
    Cursor& cursor = cursors_[mate];
    for (;;) {
        if (cursor.dead) {
            return;
        }
        if (cursor.batch == nullptr) {
            cursor.batch = readers_[mate]->popReady();
            cursor.index = 0;
            if (cursor.batch == nullptr) {
                cursor.dead = true;
                return;
            }
        }
        if (cursor.index < cursor.batch->records.size() ||
            cursor.batch->end != FastxMateEnd::None) {
            return;
        }
        readers_[mate]->recycle(cursor.batch);
        cursor.batch = nullptr;
    }
}

bool FastxMateReaderGroup::hasRecord(uint32_t mate) const {
    const Cursor& cursor = cursors_[mate];
    return !cursor.dead && cursor.batch != nullptr &&
           cursor.index < cursor.batch->records.size();
}

FastxMateEnd FastxMateReaderGroup::endOf(uint32_t mate) const {
    const Cursor& cursor = cursors_[mate];
    if (cursor.dead || cursor.batch == nullptr) {
        return FastxMateEnd::Input;
    }
    return cursor.batch->end;
}

std::string FastxMateReaderGroup::laneCount(uint32_t mate) const {
    const Cursor& cursor = cursors_[mate];
    if (cursor.dead || cursor.batch == nullptr) {
        return "?";
    }
    return std::to_string(cursor.batch->laneRecords);
}

// Skips the rest of this mate's current lane, up to its lane or input end.
void FastxMateReaderGroup::drainLane(uint32_t mate) {
    for (;;) {
        normalize(mate);
        Cursor& cursor = cursors_[mate];
        if (cursor.dead) {
            return;
        }
        cursor.index = cursor.batch->records.size();
        if (cursor.batch->end != FastxMateEnd::None) {
            return;
        }
    }
}

FastxMateReaderGroup::PairStatus FastxMateReaderGroup::nextPair(FastxChunkFillContext& context,
                                                                int* laneOut, int* endCharOut) {
    const uint32_t n = mates();
    for (;;) {
        for (uint32_t m = 0; m < n; ++m) {
            normalize(m);
        }
        const int lane = cursors_[0].batch != nullptr ? cursors_[0].batch->lane
                                                      : *context.readFilesIndex;

        if (hasRecord(0)) {
            uint32_t short_ = n;
            for (uint32_t m = 1; m < n; ++m) {
                if (!hasRecord(m)) {
                    short_ = m;
                    break;
                }
            }
            if (short_ == n) {
                return PairStatus::Pair;
            }
            for (uint32_t m = 1; m < n; ++m) {
                if (hasRecord(m)) {
                    continue;
                }
                const FastxMateEnd end = endOf(m);
                if (end == FastxMateEnd::Error) {
                    context.onFatal(cursors_[m].batch->errorText);
                    return PairStatus::Error;
                }
                if (end == FastxMateEnd::Input) {
                    if (!extraWarned_ && context.onWarning) {
                        context.onWarning(" FASTX mate " + std::to_string(m + 1) + " ends after " +
                            laneCount(m) + " reads of input file # " + std::to_string(lane) +
                            ", before mate 1. STAR skips the remaining reads of the other mates.\n");
                    }
                    extraWarned_ = true;
                    *endCharOut = kNoEndLog;
                    return PairStatus::End;
                }
            }
            // A mate reached the end of this input file first: use the reads
            // all mates have and skip the rest of the file.
            for (uint32_t m = 0; m < n; ++m) {
                drainLane(m);
            }
            ++truncatedLanes_;
            if (context.onWarning) {
                std::string counts;
                for (uint32_t m = 0; m < n; ++m) {
                    counts += (m ? ", mate " : "mate ") + std::to_string(m + 1) + ": " + laneCount(m);
                }
                context.onWarning(" FASTX mates have different numbers of reads in input file # " +
                    std::to_string(lane) + " (" + counts + "). STAR mapped the reads all mates have "
                    "and skipped the rest of this file.\n");
            }
            continue;
        }

        const FastxMateEnd end0 = endOf(0);
        if (end0 == FastxMateEnd::Error) {
            const FastxMateBatch& batch = *cursors_[0].batch;
            if (!batch.errorText.empty()) {
                context.onFatal(batch.errorText);
            } else {
                context.onBadRecordStart(batch.errorWord, batch.errorRest);
            }
            return PairStatus::Error;
        }
        if (end0 == FastxMateEnd::Input) {
            if (!extraWarned_) {
                for (uint32_t m = 1; m < n; ++m) {
                    if (hasRecord(m) || endOf(m) == FastxMateEnd::Lane) {
                        if (context.onWarning) {
                            context.onWarning(" FASTX mate 1 ends after " + laneCount(0) +
                                " reads of input file # " + std::to_string(lane) + ", before mate " +
                                std::to_string(m + 1) + ". STAR skips the remaining reads of the other mates.\n");
                        }
                        extraWarned_ = true;
                        break;
                    }
                }
            }
            *endCharOut = cursors_[0].dead ? -1 : cursors_[0].batch->endChar;
            // The single-threaded loop logged the end of input only when mates
            // 1 and 2 were both still good() at the record start.
            const bool mate0Silent = !cursors_[0].dead && cursors_[0].batch->endSilent;
            const bool mate1Silent = n > 1 && !hasRecord(1) && endOf(1) == FastxMateEnd::Input &&
                                     !cursors_[1].dead && cursors_[1].batch->endSilent;
            if (mate0Silent || mate1Silent) {
                *endCharOut = kNoEndLog;
            }
            return PairStatus::End;
        }

        // Mate 0 is at a FILE marker.
        bool drained = false;
        for (uint32_t m = 1; m < n; ++m) {
            if (hasRecord(m)) {
                drainLane(m);
                drained = true;
            }
        }
        if (drained) {
            ++truncatedLanes_;
            if (context.onWarning) {
                std::string counts;
                for (uint32_t m = 0; m < n; ++m) {
                    counts += (m ? ", mate " : "mate ") + std::to_string(m + 1) + ": " + laneCount(m);
                }
                context.onWarning(" FASTX mates have different numbers of reads in input file # " +
                    std::to_string(lane) + " (" + counts + "). STAR mapped the reads all mates have "
                    "and skipped the rest of this file.\n");
            }
            continue;
        }
        const int next = cursors_[0].batch->nextLane;
        for (uint32_t m = 1; m < n; ++m) {
            const FastxMateEnd end = endOf(m);
            if (end == FastxMateEnd::Error) {
                context.onFatal(cursors_[m].batch->errorText);
                return PairStatus::Error;
            }
            if (end == FastxMateEnd::Input) {
                if (!extraWarned_ && context.onWarning) {
                    context.onWarning(" FASTX mate " + std::to_string(m + 1) + " has no input file # " +
                        std::to_string(next) + ". STAR skips the remaining files of the other mates.\n");
                }
                extraWarned_ = true;
                *endCharOut = kNoEndLog;
                return PairStatus::End;
            }
            if (cursors_[m].batch->nextLane != next) {
                context.onFatal("EXITING because of FATAL INPUT ERROR: the mates' read files are out of step: "
                    "mate 1 starts input file # " + std::to_string(next) + " where mate " +
                    std::to_string(m + 1) + " starts input file # " +
                    std::to_string(cursors_[m].batch->nextLane) + ".\n"
                    "SOLUTION: give the same number of files for every mate in --readFilesIn.\n");
                return PairStatus::Error;
            }
        }
        for (uint32_t m = 0; m < n; ++m) {
            readers_[m]->recycle(cursors_[m].batch);
            cursors_[m].batch = nullptr;
            cursors_[m].index = 0;
        }
        *laneOut = next;
        return PairStatus::Lane;
    }
}

bool FastxMateReaderGroup::appendPair(char* const* chunkIn, unsigned long long* totals,
                                      FastxChunkFillContext& context) {
    const uint32_t n = mates();
    const FastxMateBatch* batches[kFastxMateMaxMates];
    const FastxMateRecord* records[kFastxMateMaxMates];
    for (uint32_t m = 0; m < n; ++m) {
        batches[m] = cursors_[m].batch;
        records[m] = &cursors_[m].batch->records[cursors_[m].index];
    }
    const std::string readNumber = std::to_string(*context.iReadAll);
    for (uint32_t m = 1; m < n; ++m) {
        if (records[m]->format != records[0]->format) {
            context.onFatal("EXITING because of FATAL INPUT ERROR: read " + readNumber + " is " +
                formatName(records[0]->format) + " in mate 1 but " + formatName(records[m]->format) +
                " in mate " + std::to_string(m + 1) + ".\n"
                "SOLUTION: give all mates in the same format.\n");
            return false;
        }
    }
    const std::string lane = std::to_string(*context.readFilesIndex);
    auto fits = [&](uint32_t m, unsigned long long need) {
        if (totals[m] + need + 1 <= context.chunkArrayBytes) {
            return true;
        }
        context.onFatal("EXITING because of FATAL INPUT ERROR: read " + readNumber + " (mate " +
            std::to_string(m + 1) + ") does not fit in the input chunk buffer.\n"
            "SOLUTION: increase --limitIObufferSize.\n");
        return false;
    };
    auto put = [](char*& out, const char* data, size_t size) {
        if (size != 0) {
            std::memcpy(out, data, size);
            out += size;
        }
    };

    if (records[0]->format == '@') {
        // Every mate's header carries mate 1's ID and filter flag.
        std::string core = context.fastqReadIdNumber
            ? "@" + readNumber
            : std::string(batches[0]->at(records[0]->idOff), records[0]->idLen);
        core += ' ';
        core += readNumber;
        core += ' ';
        core += records[0]->filter;
        core += ' ';
        core += lane;
        for (uint32_t m = 0; m < n; ++m) {
            const FastxMateRecord& record = *records[m];
            const FastxMateBatch& batch = *batches[m];
            const unsigned long long need = core.size() +
                (record.extraLen != 0 ? 1ULL + record.extraLen : 0ULL) + 1 +
                record.seqLen + 2 + record.qualLen;
            if (!fits(m, need)) {
                return false;
            }
            char* out = chunkIn[m] + totals[m];
            put(out, core.data(), core.size());
            if (record.extraLen != 0) {
                *out++ = ' ';
                put(out, batch.at(record.extraOff), record.extraLen);
            }
            *out++ = '\n';
            put(out, batch.at(record.seqOff), record.seqLen);
            *out++ = '+';
            *out++ = '\n';
            put(out, batch.at(record.qualOff), record.qualLen);
            totals[m] += need;
        }
        return true;
    }

    const std::string tail = " " + readNumber + " N " + lane + " \n";
    for (uint32_t m = 0; m < n; ++m) {
        const FastxMateRecord& record = *records[m];
        const FastxMateBatch& batch = *batches[m];
        const std::string id = context.fastaReadIdNumber
            ? ">" + readNumber
            : std::string(batch.at(record.idOff), record.idLen);
        const unsigned long long need = id.size() + tail.size() + record.seqLen;
        if (!fits(m, need)) {
            return false;
        }
        char* out = chunkIn[m] + totals[m];
        put(out, id.data(), id.size());
        put(out, tail.data(), tail.size());
        put(out, batch.at(record.seqOff), record.seqLen);
        totals[m] += need;
    }
    return true;
}

void FastxMateReaderGroup::fillChunk(char* const* chunkIn, unsigned long long* totals,
                                     FastxChunkFillContext& context) {
    ensureStarted();
    const uint32_t n = mates();
    while (totals[0] < context.chunkInSizeBytes && totals[1] < context.chunkInSizeBytes) {
        if (*context.iReadAll == context.readMapNumber) {
            break;
        }
        int lane = 0;
        int endChar = kNoEndLog;
        const PairStatus status = nextPair(context, &lane, &endChar);
        if (status == PairStatus::Pair) {
            ++*context.iReadAll;
            if (!appendPair(chunkIn, totals, context)) {
                return;
            }
            for (uint32_t m = 0; m < n; ++m) {
                ++cursors_[m].index;
            }
            ++pairs_;
            continue;
        }
        if (status == PairStatus::Lane) {
            *context.readFilesIndex = lane;
            if (context.onLaneStart) {
                context.onLaneStart(lane);
            }
            continue;
        }
        if (status == PairStatus::End) {
            // The single-threaded loop logged this each time mate 1 showed a
            // blank line, and once at the end of its stream.
            if (endChar != kNoEndLog && (!endLogged_ || endChar != -1) && context.onLog) {
                context.onLog("Thread #" + std::to_string(context.thread) +
                              " end of input stream, nextChar=" + std::to_string(endChar));
            }
            if (endChar != kNoEndLog) {
                endLogged_ = true;
            }
            break;
        }
        return;  // PairStatus::Error: reported through the callbacks
    }
}

std::string FastxMateReaderGroup::summary() const {
    std::ostringstream out;
    out << "Fastx mate readers: " << pairs_ << " reads";
    if (truncatedLanes_ != 0) {
        out << "; " << truncatedLanes_ << " input file(s) truncated to the shorter mate";
    }
    out << "\n" << std::fixed << std::setprecision(3);
    for (uint32_t m = 0; m < mates(); ++m) {
        const FastxMateReader::Stats s = readers_[m]->stats();
        out << "  mate " << (m + 1)
            << ": records=" << s.records
            << " batches=" << s.batches
            << " inputMiB=" << (static_cast<double>(s.bytes) / 1048576.0)
            << " parseSeconds=" << (static_cast<double>(s.parseNs) / 1e9)
            << " inputWaitSeconds=" << (static_cast<double>(s.inputWaitNs) / 1e9)
            << " queueFullSeconds=" << (static_cast<double>(s.freeWaitNs) / 1e9)
            << " fillerWaitSeconds=" << (static_cast<double>(s.consumerWaitNs) / 1e9)
            << " permits=" << s.permitAcquires
            << " permitWaitSeconds=" << (static_cast<double>(s.permitWaitNs) / 1e9)
            << "\n";
    }
    return out.str();
}

}  // namespace input
}  // namespace star

// Harness for the per-mate FASTX reader threads (input/FastxMateReaders).
//
// Part 1 compares FastxMateReaderGroup::fillChunk with a copy of the
// single-threaded chunk fill in ReadAlignChunk_processChunks.cpp (FASTQ,
// FASTA, FILE markers, the bad-record error), on inputs where the two must
// agree byte for byte: chunk text, chunk boundaries, iReadAll, readFilesIndex
// and the "Starting to map file" and "end of input stream" Log.out lines. Each case runs with several chunk
// sizes, with and without --outSAMreadID Number and a --readMapNumber limit,
// and with and without a one-permit pool.
//
// Part 2 checks the intended behaviour where the single-threaded loop had
// none: mates with different read counts (truncate per file, warn), lane
// markers out of step (fatal), mates in different formats (fatal).
//
// Part 3 checks the deliberate change for blank lines where a read header is
// expected (the author's decisions of 29 and 30 Sep), in each of three mates,
// for FASTQ and FASTA: a single blank line before a header, and any number of
// blank lines at the end of a file, are skipped with one WARNING per file and
// counted per file; a second blank line in a row before the end of the file,
// or a blank line before anything but a header, is fatal, naming file and line.
//
// Usage: fastx_mate_reader_harness <scratch dir>

#include "input/FastxMateReaders.h"

#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using star::input::BgzfWorkPermitHooks;
using star::input::FastxChunkFillContext;
using star::input::FastxMateLimits;
using star::input::FastxMateReaderGroup;

namespace {

// STAR's defaults (IncludeDefine.h without COMPILE_FOR_LONG_READS).
constexpr long long kNameLengthMax = 50000;
constexpr long long kSeqLengthMax = 650;
constexpr long long kNameSeqLengthMax = kNameLengthMax > kSeqLengthMax ? kNameLengthMax : kSeqLengthMax;

struct Options {
    unsigned long long chunkInSizeBytes = 1ULL << 20;
    unsigned long long readMapNumber = static_cast<unsigned long long>(-1);
    bool readIdNumber = false;
    bool permits = false;
};

struct Chunk {
    std::string mate[3];
    unsigned long long iReadAllAfter = 0;
    int laneAfter = 0;
};

struct Result {
    std::vector<Chunk> chunks;
    std::vector<int> lanes;
    bool badRecord = false;
    std::string badWord, badRest;
    unsigned long long badRead = 0;
    bool fatal = false;
    std::string fatalText;
    std::vector<std::string> warnings;
    std::vector<std::string> logs;  // "end of input stream" lines
    std::string summary;            // FastxMateReaderGroup::summary(), new path only
    unsigned long long reads = 0;
};

unsigned long long arrayBytes(const Options& options) {
    // Parameters.cpp: chunkInSizeBytes = array - 2*(DEF_readSeqLengthMax+1) - 2*DEF_readNameLengthMax
    return options.chunkInSizeBytes + 2 * (kSeqLengthMax + 1) + 2 * kNameLengthMax;
}

// ---- copy of the single-threaded loop's helpers ----------------------------

std::string legacyHeaderExtra(std::istream& in) {
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

char legacyFilter(const std::string& extra) {
    std::string field2;
    std::istringstream extraStream(extra);
    extraStream >> field2;
    if (field2.length() >= 4 && field2[1] == ':' && field2[2] == 'Y' && field2[3] == ':') {
        return 'Y';
    }
    return 'N';
}

// fastqReadOneLine, with the index arithmetic done in signed form: arrIn
// always points past a header or "+" line, so the bytes before it exist.
long long legacyReadOneLine(std::istream& in, char* arrIn) {
    in.getline(arrIn, kNameSeqLengthMax + 1);
    long long lenIn = in.gcount();
    if (static_cast<int>(arrIn[lenIn - 2]) < 33) {
        --lenIn;
    }
    arrIn[lenIn - 1] = '\n';
    return lenIn;
}

// ---- copy of the single-threaded fill (ReadAlignChunk_processChunks.cpp) ----

Result runLegacy(const std::vector<std::string>& files, const Options& options) {
    Result result;
    const size_t mates = files.size();
    std::ifstream in[3];
    for (size_t m = 0; m < mates; ++m) {
        in[m].open(files[m].c_str());
    }
    unsigned long long iReadAll = 0;
    int readFilesIndex = 0;
    const unsigned long long array = arrayBytes(options);
    std::vector<std::vector<char>> buffers(3, std::vector<char>(array, '\n'));
    char* chunkIn[3] = {buffers[0].data(), buffers[1].data(), buffers[2].data()};
    for (;;) {
        unsigned long long total[3] = {0, 0, 0};
        bool newFile = false;
        while (total[0] < options.chunkInSizeBytes && total[1] < options.chunkInSizeBytes &&
               in[0].good() && in[1].good()) {
            char nextChar = in[0].peek();
            if (iReadAll == options.readMapNumber) {
                break;
            } else if (nextChar == '@') {
                iReadAll++;
                std::string readID;
                in[0] >> readID;
                if (static_cast<int>(readID.back()) < 33) {
                    readID.pop_back();
                }
                std::vector<std::string> extras(mates);
                extras[0] = legacyHeaderExtra(in[0]);
                if (options.readIdNumber) {
                    readID = "@" + std::to_string(iReadAll);
                }
                const char filter = legacyFilter(extras[0]);
                for (size_t m = 1; m < mates; ++m) {
                    std::string mateReadID;
                    in[m] >> mateReadID;
                    extras[m] = legacyHeaderExtra(in[m]);
                }
                for (size_t m = 0; m < mates; ++m) {
                    std::string header = readID + ' ' + std::to_string(iReadAll) + ' ' + filter + ' ' +
                                         std::to_string(readFilesIndex);
                    if (!extras[m].empty()) {
                        header += ' ' + extras[m];
                    }
                    total[m] += 1 + header.copy(chunkIn[m] + total[m], header.size(), 0);
                    chunkIn[m][total[m] - 1] = '\n';
                }
                for (size_t m = 0; m < mates; ++m) {
                    total[m] += legacyReadOneLine(in[m], chunkIn[m] + total[m]);
                    in[m].ignore(kNameSeqLengthMax, '\n');
                    chunkIn[m][total[m]] = '+';
                    chunkIn[m][total[m] + 1] = '\n';
                    total[m] += 2;
                    total[m] += legacyReadOneLine(in[m], chunkIn[m] + total[m]);
                }
            } else if (nextChar == '>') {
                iReadAll++;
                for (size_t m = 0; m < mates; ++m) {
                    std::string head;
                    if (options.readIdNumber) {
                        head = ">" + std::to_string(iReadAll);
                    } else {
                        in[m] >> head;
                    }
                    std::memcpy(chunkIn[m] + total[m], head.data(), head.size());
                    total[m] += head.size();
                    in[m].ignore(kNameSeqLengthMax, '\n');
                    total[m] += std::sprintf(chunkIn[m] + total[m], " %llu %c %i \n", iReadAll, 'N', readFilesIndex);
                    char c = in[m].peek();
                    while (c != '@' && c != '>' && c != ' ' && c != '\n' && in[m].good()) {
                        in[m].getline(chunkIn[m] + total[m], kSeqLengthMax + 1);
                        if (in[m].gcount() < 2) {
                            break;
                        }
                        total[m] += in[m].gcount() - 1;
                        if (static_cast<int>(chunkIn[m][total[m] - 1]) < 33) {
                            total[m]--;
                        }
                        c = in[m].peek();
                    }
                    chunkIn[m][total[m]] = '\n';
                    total[m]++;
                }
            } else if (nextChar == ' ' || nextChar == '\n' || !in[0].good()) {
                result.logs.push_back("Thread #0 end of input stream, nextChar=" +
                                      std::to_string(static_cast<int>(nextChar)));
                break;
            } else {
                std::string word1;
                in[0] >> word1;
                if (word1 == "FILE") {
                    newFile = true;
                } else {
                    std::string str1;
                    std::getline(in[0], str1);
                    result.badRecord = true;
                    result.badWord = word1;
                    result.badRest = str1;
                    result.badRead = iReadAll + 1;
                    result.reads = iReadAll;
                    return result;
                }
            }
            if (newFile) {
                in[0] >> readFilesIndex;
                result.lanes.push_back(readFilesIndex);
                for (size_t m = 0; m < mates; ++m) {
                    in[m].ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }
                newFile = false;
            }
        }
        if (mates == 2 && ((total[0] == 0) != (total[1] == 0))) {
            result.fatal = true;
            result.fatalText = "paired mates have unequal FASTQ buffering";
            break;
        }
        if (total[0] == 0) {
            break;
        }
        Chunk chunk;
        for (size_t m = 0; m < mates; ++m) {
            chunk.mate[m].assign(chunkIn[m], total[m]);
        }
        chunk.iReadAllAfter = iReadAll;
        chunk.laneAfter = readFilesIndex;
        result.chunks.push_back(chunk);
    }
    result.reads = iReadAll;
    return result;
}

// ---- the module under test ---------------------------------------------------

struct OnePermitPool {
    std::mutex mutex;
    std::condition_variable cv;
    int available = 1;
    unsigned long long acquires = 0;
};
OnePermitPool g_pool;

uint64_t poolAcquire(void*) {
    std::unique_lock<std::mutex> lock(g_pool.mutex);
    g_pool.cv.wait(lock, [] { return g_pool.available > 0; });
    --g_pool.available;
    ++g_pool.acquires;
    return 0;
}

void poolRelease(void*, uint64_t, uint64_t, uint64_t, uint64_t) {
    {
        std::lock_guard<std::mutex> lock(g_pool.mutex);
        ++g_pool.available;
    }
    g_pool.cv.notify_one();
}

// Names used in messages: "<path>#<input file index>".
std::string laneFileName(const std::string& path, int lane) {
    return path + "#" + std::to_string(lane);
}

Result runThreads(const std::vector<std::string>& files, const Options& options) {
    Result result;
    const size_t mates = files.size();
    std::ifstream in[3];
    std::vector<std::istream*> streams;
    std::vector<std::vector<std::string>> names(mates);
    for (size_t m = 0; m < mates; ++m) {
        in[m].open(files[m].c_str());
        streams.push_back(&in[m]);
        for (int lane = 0; lane < 8; ++lane) {
            names[m].push_back(laneFileName(files[m], lane));
        }
    }
    BgzfWorkPermitHooks hooks;
    if (options.permits) {
        hooks.acquire = poolAcquire;
        hooks.release = poolRelease;
    }
    FastxMateLimits limits;
    limits.nameSeqLineMax = kNameSeqLengthMax;
    limits.seqLineMax = kSeqLengthMax;
    unsigned long long iReadAll = 0;
    int readFilesIndex = 0;
    const unsigned long long array = arrayBytes(options);
    std::vector<std::vector<char>> buffers(3, std::vector<char>(array, '\n'));
    char* chunkIn[3] = {buffers[0].data(), buffers[1].data(), buffers[2].data()};
    bool stop = false;
    {
        FastxMateReaderGroup group(streams, 0, hooks, limits, names);
        FastxChunkFillContext context;
        context.chunkInSizeBytes = options.chunkInSizeBytes;
        context.chunkArrayBytes = array;
        context.readMapNumber = options.readMapNumber;
        context.fastqReadIdNumber = options.readIdNumber;
        context.fastaReadIdNumber = options.readIdNumber;
        context.iReadAll = &iReadAll;
        context.readFilesIndex = &readFilesIndex;
        context.onLaneStart = [&](int lane) { result.lanes.push_back(lane); };
        context.onLog = [&](const std::string& line) { result.logs.push_back(line); };
        context.onWarning = [&](const std::string& text) { result.warnings.push_back(text); };
        context.onBadRecordStart = [&](const std::string& word, const std::string& rest) {
            result.badRecord = true;
            result.badWord = word;
            result.badRest = rest;
            result.badRead = iReadAll + 1;
            stop = true;
        };
        context.onFatal = [&](const std::string& text) {
            result.fatal = true;
            result.fatalText = text;
            stop = true;
        };
        for (;;) {
            unsigned long long total[3] = {0, 0, 0};
            group.fillChunk(chunkIn, total, context);
            if (stop || total[0] == 0) {
                break;
            }
            Chunk chunk;
            for (size_t m = 0; m < mates; ++m) {
                chunk.mate[m].assign(chunkIn[m], total[m]);
            }
            chunk.iReadAllAfter = iReadAll;
            chunk.laneAfter = readFilesIndex;
            result.chunks.push_back(chunk);
        }
        group.stopAndJoin();
        result.summary = group.summary();
    }
    result.reads = iReadAll;
    return result;
}

// ---- inputs ------------------------------------------------------------------

std::string bases(std::mt19937& rng, size_t length) {
    static const char kBases[] = "ACGTN";
    std::uniform_int_distribution<int> pick(0, 3);
    std::string out(length, 'A');
    for (auto& c : out) {
        c = kBases[pick(rng)];
    }
    return out;
}

std::string fastq(const std::string& name, const std::string& extra, const std::string& seq,
                  const std::string& eol = "\n") {
    return "@" + name + (extra.empty() ? std::string() : " " + extra) + eol + seq + eol + "+" + eol +
           std::string(seq.size(), 'I') + eol;
}

std::string fasta(const std::string& name, const std::string& seq, size_t width) {
    std::string out = ">" + name + " description\n";
    for (size_t at = 0; at < seq.size(); at += width) {
        out += seq.substr(at, width) + "\n";
    }
    return out;
}

// Paired FASTQ text: mate 0 of length len0, mate 1 of length len1.
void pairedFastq(std::mt19937& rng, size_t count, size_t first, size_t len0, size_t len1,
                 std::string* mate0, std::string* mate1, const std::string& eol = "\n") {
    for (size_t i = 0; i < count; ++i) {
        const std::string name = "r" + std::to_string(first + i);
        const std::string flag = (i % 17 == 3) ? "Y" : "N";
        *mate0 += fastq(name, "1:" + flag + ":0:ACGTACGT", bases(rng, len0), eol);
        if (mate1 != nullptr) {
            *mate1 += fastq(name, "2:" + flag + ":0:ACGTACGT", bases(rng, len1), eol);
        }
    }
}

struct Case {
    Case(const std::string& nameIn, const std::vector<std::string>& matesIn)
        : name(nameIn), mates(matesIn) {}
    std::string name;
    std::vector<std::string> mates;
    // The new path warns here by design (D2: one mate's stream ends first).
    bool expectWarning = false;
};

std::vector<Case> identityCases() {
    std::mt19937 rng(20260929);
    std::vector<Case> cases;
    {
        Case c{"pe_basic", {"", ""}};
        pairedFastq(rng, 3000, 0, 90, 28, &c.mates[0], &c.mates[1]);
        cases.push_back(c);
    }
    {
        Case c{"pe_extras", {"", ""}};
        const char* extras[] = {"1:Y:0:AC", "", "\tBC:Z:AAAA", "  3:N:0:x  ", "1:N:0:1\tCB:Z:x", "a:Y:b"};
        for (size_t i = 0; i < 600; ++i) {
            const std::string name = "e" + std::to_string(i);
            c.mates[0] += fastq(name, extras[i % 6], bases(rng, 50));
            c.mates[1] += fastq(name, extras[(i + 1) % 6], bases(rng, 30));
        }
        cases.push_back(c);
    }
    {
        Case c{"pe_crlf", {"", ""}};
        pairedFastq(rng, 400, 0, 60, 28, &c.mates[0], &c.mates[1], "\r\n");
        cases.push_back(c);
    }
    {
        Case c{"pe_lanes", {"", ""}};
        for (size_t m = 0; m < 2; ++m) {
            c.mates[m] += "FILE 0\n";
        }
        pairedFastq(rng, 500, 0, 90, 28, &c.mates[0], &c.mates[1]);
        for (size_t m = 0; m < 2; ++m) {
            c.mates[m] += "FILE 1\nFILE 2\n";  // lane 1 is empty
        }
        pairedFastq(rng, 700, 500, 90, 28, &c.mates[0], &c.mates[1]);
        cases.push_back(c);
    }
    {
        Case c{"se", {""}};
        pairedFastq(rng, 1200, 0, 75, 0, &c.mates[0], nullptr);
        cases.push_back(c);
    }
    {
        Case c{"three_mates", {"", "", ""}};
        for (size_t i = 0; i < 900; ++i) {
            const std::string name = "t" + std::to_string(i);
            c.mates[0] += fastq(name, "1:N:0:A", bases(rng, 80));
            c.mates[1] += fastq(name, "2:N:0:A", bases(rng, 28));
            c.mates[2] += fastq(name, "3:N:0:A", bases(rng, 10));
        }
        cases.push_back(c);
    }
    {
        Case c{"fasta", {"", ""}};
        for (size_t i = 0; i < 800; ++i) {
            c.mates[0] += fasta("f" + std::to_string(i), bases(rng, 70), 1000);
            c.mates[1] += fasta("f" + std::to_string(i), bases(rng, 40), 1000);
        }
        cases.push_back(c);
    }
    {
        Case c{"fasta_multiline_lanes", {"FILE 0\n", "FILE 0\n"}};
        for (size_t i = 0; i < 300; ++i) {
            c.mates[0] += fasta("g" + std::to_string(i), bases(rng, 95), 30);
            c.mates[1] += fasta("g" + std::to_string(i), bases(rng, 61), 20);
        }
        c.mates[0] += "FILE 1\n";
        c.mates[1] += "FILE 1\n";
        for (size_t i = 0; i < 200; ++i) {
            c.mates[0] += fasta("h" + std::to_string(i), bases(rng, 33), 7);
            c.mates[1] += fasta("h" + std::to_string(i), bases(rng, 33), 100);
        }
        cases.push_back(c);
    }
    {
        Case c{"no_final_newline", {"", ""}};
        pairedFastq(rng, 300, 0, 90, 28, &c.mates[0], &c.mates[1]);
        c.mates[0].pop_back();
        c.mates[1].pop_back();
        cases.push_back(c);
    }
    {
        Case c{"control_and_high_bytes", {"", ""}};
        for (size_t i = 0; i < 300; ++i) {
            const std::string suffix = (i % 3 == 0) ? "\x80" : ((i % 3 == 1) ? "\x01" : "");
            const std::string name = "h" + std::to_string(i) + suffix;
            c.mates[0] += fastq(name, i % 2 ? "1:N:0:A\x7f" : "", bases(rng, 90));
            c.mates[1] += fastq(name, "2:N:0:A", bases(rng, 28));
        }
        cases.push_back(c);
    }
    {
        Case c{"long_extras", {"", ""}};
        for (size_t i = 0; i < 200; ++i) {
            const std::string name = "l" + std::to_string(i);
            c.mates[0] += fastq(name, "1:N:0:" + std::string(3000 + i, 'X'), bases(rng, 90));
            c.mates[1] += fastq(name, "2:N:0:" + std::string(10, 'Y'), bases(rng, 28));
        }
        cases.push_back(c);
    }
    {
        // Only mate 2 lacks a final newline.
        Case c{"mate2_no_final_newline", {"", ""}};
        pairedFastq(rng, 300, 0, 90, 28, &c.mates[0], &c.mates[1]);
        c.mates[1].pop_back();
        cases.push_back(c);
    }
    {
        // Lines at STAR's limits: FASTA lines of DEF_readSeqLengthMax bytes,
        // FASTQ reads of that length and header lines near
        // DEF_readNameLengthMax bytes.
        Case c{"line_limits_fasta", {"", ""}};
        for (size_t i = 0; i < 60; ++i) {
            const std::string name = "L" + std::to_string(i);
            c.mates[0] += fasta(name, bases(rng, 3 * kSeqLengthMax), kSeqLengthMax);
            c.mates[1] += fasta(name, bases(rng, kSeqLengthMax), kSeqLengthMax);
        }
        cases.push_back(c);
    }
    {
        Case c{"line_limits_fastq", {"", ""}};
        for (size_t i = 0; i < 60; ++i) {
            const std::string name = "M" + std::to_string(i);
            const std::string pad(static_cast<size_t>(kNameLengthMax) - 64 - name.size(), 'X');
            c.mates[0] += fastq(name, "1:N:0:" + pad, bases(rng, kSeqLengthMax));
            c.mates[1] += fastq(name, "2:N:0:A", bases(rng, kSeqLengthMax));
        }
        cases.push_back(c);
    }
    {
        // A FASTA line one byte over DEF_readSeqLengthMax in mate 1 stops
        // that stream after the record, in both paths.
        Case c{"fasta_line_over_limit", {"", ""}};
        c.expectWarning = true;
        for (size_t i = 0; i < 40; ++i) {
            const std::string name = "O" + std::to_string(i);
            c.mates[0] += fasta(name, bases(rng, i == 25 ? kSeqLengthMax + 1 : 80), 1000);
            c.mates[1] += fasta(name, bases(rng, 40), 1000);
        }
        cases.push_back(c);
    }
    {
        // Several refills of the reader's 1 MiB input buffer.
        Case c{"big", {"", ""}};
        pairedFastq(rng, 40000, 0, 90, 28, &c.mates[0], &c.mates[1]);
        cases.push_back(c);
    }
    {
        // Mate 0 has a line that is neither a read nor a FILE marker.
        Case c{"bad_record", {"", ""}};
        pairedFastq(rng, 250, 0, 90, 28, &c.mates[0], &c.mates[1]);
        c.mates[0] += "garbage line here\n";
        pairedFastq(rng, 50, 250, 90, 28, &c.mates[0], &c.mates[1]);
        cases.push_back(c);
    }
    return cases;
}

std::vector<std::string> writeCase(const std::string& dir, const Case& c) {
    std::vector<std::string> files;
    for (size_t m = 0; m < c.mates.size(); ++m) {
        const std::string path = dir + "/" + c.name + "_mate" + std::to_string(m + 1) + ".txt";
        std::ofstream out(path.c_str(), std::ios::binary);
        out << c.mates[m];
        files.push_back(path);
    }
    return files;
}

std::string describe(const Options& o) {
    std::ostringstream out;
    out << "chunk=" << o.chunkInSizeBytes << " readIdNumber=" << o.readIdNumber
        << " readMapNumber=" << (o.readMapNumber == static_cast<unsigned long long>(-1) ? std::string("all")
                                                                                          : std::to_string(o.readMapNumber))
        << " permits=" << o.permits;
    return out.str();
}

std::string compare(const Result& a, const Result& b, size_t mates) {
    if (a.chunks.size() != b.chunks.size()) {
        return "chunk count " + std::to_string(a.chunks.size()) + " vs " + std::to_string(b.chunks.size());
    }
    for (size_t i = 0; i < a.chunks.size(); ++i) {
        for (size_t m = 0; m < mates; ++m) {
            if (a.chunks[i].mate[m] != b.chunks[i].mate[m]) {
                return "chunk " + std::to_string(i) + " mate " + std::to_string(m + 1) + " text differs";
            }
        }
        if (a.chunks[i].iReadAllAfter != b.chunks[i].iReadAllAfter) {
            return "chunk " + std::to_string(i) + " iReadAll differs";
        }
        if (a.chunks[i].laneAfter != b.chunks[i].laneAfter) {
            return "chunk " + std::to_string(i) + " readFilesIndex differs";
        }
    }
    if (a.lanes != b.lanes) {
        return "lane events differ";
    }
    if (a.logs != b.logs) {
        return "end-of-input log lines differ (" + std::to_string(a.logs.size()) + " vs " +
               std::to_string(b.logs.size()) + ")";
    }
    if (a.badRecord != b.badRecord || a.badWord != b.badWord || a.badRest != b.badRest ||
        a.badRead != b.badRead) {
        return "bad-record error differs";
    }
    if (a.fatal != b.fatal) {
        return "fatal state differs: " + a.fatalText + " | " + b.fatalText;
    }
    if (a.reads != b.reads) {
        return "read count " + std::to_string(a.reads) + " vs " + std::to_string(b.reads);
    }
    return "";
}

int identityPart(const std::string& dir) {
    int failures = 0;
    int runs = 0;
    const unsigned long long chunkSizes[] = {1ULL, 997ULL, 65536ULL, 1ULL << 20};
    for (const Case& c : identityCases()) {
        const std::vector<std::string> files = writeCase(dir, c);
        int variant = 0;
        for (unsigned long long chunk : chunkSizes) {
            for (int readIdNumber = 0; readIdNumber < 2; ++readIdNumber) {
                for (int limited = 0; limited < 2; ++limited) {
                    Options options;
                    options.chunkInSizeBytes = chunk;
                    options.readIdNumber = readIdNumber != 0;
                    options.readMapNumber = limited ? 777ULL : static_cast<unsigned long long>(-1);
                    options.permits = (variant++ % 2) == 1;
                    const Result legacy = runLegacy(files, options);
                    const Result threads = runThreads(files, options);
                    std::string difference = compare(legacy, threads, files.size());
                    if (difference.empty() && !c.expectWarning && !threads.warnings.empty()) {
                        difference = "unexpected WARNING: " + threads.warnings.front();
                    }
                    ++runs;
                    if (!difference.empty()) {
                        ++failures;
                        std::cout << "FAIL identity " << c.name << " [" << describe(options) << "]: "
                                  << difference << "\n";
                    }
                }
            }
        }
        std::cout << "checked identity " << c.name << "\n";
    }
    std::cout << "identity: " << runs << " runs, " << failures << " failures\n";
    return failures;
}

// ---- intended behaviour where the single-threaded loop had none ---------------

int expect(bool condition, const std::string& what) {
    if (!condition) {
        std::cout << "FAIL behaviour " << what << "\n";
        return 1;
    }
    return 0;
}

int behaviourPart(const std::string& dir) {
    std::mt19937 rng(7);
    int failures = 0;
    Options options;
    {
        Case c{"mate2_short_in_lane0", {"FILE 0\n", "FILE 0\n"}};
        pairedFastq(rng, 300, 0, 90, 28, &c.mates[0], nullptr);
        pairedFastq(rng, 250, 0, 28, 0, &c.mates[1], nullptr);
        c.mates[0] += "FILE 1\n";
        c.mates[1] += "FILE 1\n";
        pairedFastq(rng, 200, 300, 90, 28, &c.mates[0], &c.mates[1]);
        for (int permits = 0; permits < 2; ++permits) {
            options.permits = permits != 0;
            const Result r = runThreads(writeCase(dir, c), options);
            failures += expect(r.reads == 450 && r.warnings.size() == 1 && !r.fatal &&
                               r.lanes == std::vector<int>({0, 1}), c.name);
        }
    }
    {
        Case c{"mate1_short_at_end", {"", ""}};
        pairedFastq(rng, 300, 0, 90, 28, &c.mates[0], nullptr);
        pairedFastq(rng, 350, 0, 28, 0, &c.mates[1], nullptr);
        const Result r = runThreads(writeCase(dir, c), options);
        failures += expect(r.reads == 300 && r.warnings.size() == 1 && !r.fatal, c.name);
    }
    {
        Case c{"mate2_short_at_end", {"", ""}};
        pairedFastq(rng, 350, 0, 90, 28, &c.mates[0], nullptr);
        pairedFastq(rng, 300, 0, 28, 0, &c.mates[1], nullptr);
        const Result r = runThreads(writeCase(dir, c), options);
        failures += expect(r.reads == 300 && r.warnings.size() == 1 && !r.fatal, c.name);
    }
    {
        Case c{"lanes_out_of_step", {"FILE 0\n", "FILE 0\n"}};
        pairedFastq(rng, 100, 0, 90, 28, &c.mates[0], &c.mates[1]);
        c.mates[0] += "FILE 1\n";
        c.mates[1] += "FILE 2\n";
        pairedFastq(rng, 100, 100, 90, 28, &c.mates[0], &c.mates[1]);
        const Result r = runThreads(writeCase(dir, c), options);
        failures += expect(r.fatal && r.reads == 100, c.name);
    }
    {
        Case c{"format_mismatch", {"", ""}};
        pairedFastq(rng, 10, 0, 90, 28, &c.mates[0], nullptr);
        for (size_t i = 0; i < 10; ++i) {
            c.mates[1] += fasta("r" + std::to_string(i), bases(rng, 28), 100);
        }
        const Result r = runThreads(writeCase(dir, c), options);
        failures += expect(r.fatal, c.name);
    }
    std::cout << "behaviour: " << failures << " failures\n";
    return failures;
}

// ---- blank lines where a read header is expected (deliberate change) ---------

bool contains(const std::string& text, const std::string& part) {
    return text.find(part) != std::string::npos;
}

// Record i of each of `mates` mates, in one format.
std::vector<std::vector<std::string>> blankRecords(std::mt19937& rng, char format, size_t count,
                                                   size_t mates) {
    static const size_t kLengths[3] = {90, 28, 10};
    std::vector<std::vector<std::string>> out(mates);
    for (size_t i = 0; i < count; ++i) {
        const std::string name = "k" + std::to_string(i);
        for (size_t m = 0; m < mates; ++m) {
            out[m].push_back(format == '@'
                ? fastq(name, std::to_string(m + 1) + ":N:0:A", bases(rng, kLengths[m]))
                : fasta(name, bases(rng, kLengths[m]), 1000));
        }
    }
    return out;
}

// Records [first, last) of one mate, with text inserted before the records
// named (index relative to first) and after the last one.
std::string joinRecords(const std::vector<std::string>& records, size_t first, size_t last,
                        const std::map<size_t, std::string>& before, const std::string& tail) {
    std::string out;
    for (size_t i = first; i < last; ++i) {
        const auto it = before.find(i - first);
        if (it != before.end()) {
            out += it->second;
        }
        out += records[i];
    }
    return out + tail;
}

std::string skippedLine(size_t mate, unsigned long long count, const std::string& file) {
    return "mate " + std::to_string(mate + 1) + ": skipped " + std::to_string(count) +
           " blank line(s) where a read header was expected in " + file + "\n";
}

int blankLinePart(const std::string& dir) {
    std::mt19937 rng(11);
    int failures = 0;
    Options options;
    const size_t kMates = 3;
    const size_t kCount = 300;
    const size_t kSplit = 150;  // two-file cases: file 0 holds reads [0, kSplit)
    const std::map<size_t, std::string> none;
    for (const char format : {'@', '>'}) {
        const std::string formatName = format == '@' ? "fastq" : "fasta";
        const unsigned long long linesPerRead = format == '@' ? 4 : 2;
        const std::vector<std::vector<std::string>> records =
            blankRecords(rng, format, kCount, kMates);

        // Two-file reference. FASTA: the sequence loop reads a FILE marker line
        // right after a read's sequence as sequence (v1.10.0 and here alike), so
        // every FASTA mate ends file 0 with one blank line, which ends the read.
        const std::string fileEnd = format == '@' ? "" : "\n";
        Case clean{"blank_clean_" + formatName, {}};
        Case cleanTwo{"blank_clean_two_files_" + formatName, {}};
        for (size_t m = 0; m < kMates; ++m) {
            clean.mates.push_back(joinRecords(records[m], 0, kCount, none, ""));
            cleanTwo.mates.push_back("FILE 0\n" + joinRecords(records[m], 0, kSplit, none, fileEnd) +
                                     "FILE 1\n" + joinRecords(records[m], kSplit, kCount, none, ""));
        }
        const Result reference = runThreads(writeCase(dir, clean), options);
        const Result referenceTwo = runThreads(writeCase(dir, cleanTwo), options);
        failures += expect(reference.reads == kCount && reference.warnings.empty() &&
                           !reference.fatal && !contains(reference.summary, "blank line"),
                           clean.name);
        failures += expect(referenceTwo.reads == kCount && !referenceTwo.fatal &&
                           referenceTwo.warnings.size() == (format == '@' ? 0U : kMates) &&
                           referenceTwo.lanes == std::vector<int>({0, 1}), cleanTwo.name);

        for (size_t mate = 0; mate < kMates; ++mate) {
            const std::string tag = formatName + "_mate" + std::to_string(mate + 1);
            auto variant = [&](const std::string& name, const std::map<size_t, std::string>& before,
                               const std::string& tail) {
                Case c{name + "_" + tag, clean.mates};
                c.mates[mate] = joinRecords(records[mate], 0, kCount, before, tail);
                return c;
            };

            // A blank line before a read header: skipped, one WARNING naming
            // the file and line, counted in the summary; the chunks are those
            // of the input without it.
            {
                const Case c = variant("blank_then_header", {{99, "\n"}}, "");
                const std::vector<std::string> files = writeCase(dir, c);
                for (int permits = 0; permits < 2; ++permits) {
                    options.permits = permits != 0;
                    const Result r = runThreads(files, options);
                    failures += expect(
                        compare(reference, r, kMates).empty() && r.warnings.size() == 1 &&
                        contains(r.warnings[0], laneFileName(files[mate], 0) + " (mate " +
                                 std::to_string(mate + 1) + "): blank line at line " +
                                 std::to_string(linesPerRead * 99 + 1) + " ") &&
                        contains(r.summary, skippedLine(mate, 1, laneFileName(files[mate], 0))),
                        c.name + (permits ? " [permits]" : ""));
                }
                options.permits = false;
            }

            // A blank line before anything else is fatal, naming the file and
            // the line; the reads before it are mapped.
            const std::string otherFormat = format == '@' ? ">x\nACGT\n" : "@x\nACGT\n+\nIIII\n";
            const std::pair<std::string, std::string> bad[] = {
                {"blank_then_garbage", "garbage line\n"},
                {"blank_then_indented_header", "  " + std::string(1, format) + "x 1:N:0:A\n"},
                {"blank_then_other_format", otherFormat},
            };
            for (const auto& item : bad) {
                const Case c = variant(item.first, {{99, "\n" + item.second}}, "");
                const std::vector<std::string> files = writeCase(dir, c);
                const Result r = runThreads(files, options);
                failures += expect(
                    r.fatal && r.reads == 99 && r.warnings.size() == 1 &&
                    contains(r.fatalText, "malformed input in read file " +
                             laneFileName(files[mate], 0) + " (mate " + std::to_string(mate + 1) +
                             "): line " + std::to_string(linesPerRead * 99 + 2) +
                             " follows a blank line (line " + std::to_string(linesPerRead * 99 + 1) +
                             ") but is not " + (format == '@' ? "a FASTQ" : "a FASTA")),
                    c.name);
            }

            // A second blank line in a row before the end of the file is fatal,
            // reported at the second blank line, whatever follows it.
            const std::pair<std::string, std::string> doubles[] = {
                {"double_blank_then_header", "\n\n"},
                {"double_blank_spaces_then_header", "\n \t\n"},
                {"double_blank_then_garbage", "\n\ngarbage line\n"},
            };
            for (const auto& item : doubles) {
                const Case c = variant(item.first, {{99, item.second}}, "");
                const std::vector<std::string> files = writeCase(dir, c);
                const Result r = runThreads(files, options);
                failures += expect(
                    r.fatal && r.reads == 99 && r.warnings.size() == 1 &&
                    contains(r.fatalText, "malformed input in read file " +
                             laneFileName(files[mate], 0) + " (mate " + std::to_string(mate + 1) +
                             "): line " + std::to_string(linesPerRead * 99 + 2) +
                             " is a second blank line in a row (after line " +
                             std::to_string(linesPerRead * 99 + 1) + ")"),
                    c.name);
            }

            // Several single blank lines in one file (one of spaces and tabs;
            // FASTQ also a CRLF one, while in FASTA a CRLF line after a read is
            // an empty sequence line, as in v1.10.0) and three at the end of the
            // file: one WARNING, every line counted, chunks as without them.
            {
                std::map<size_t, std::string> before = {
                    {10, "\n"}, {20, "\n"}, {30, " \t\n"}, {40, "\n"}};
                unsigned long long expected = 4 + 3;
                if (format == '@') {
                    before[50] = "\r\n";
                    ++expected;
                }
                const Case c = variant("several_blank_lines", before, "\n\n\n");
                const std::vector<std::string> files = writeCase(dir, c);
                const Result r = runThreads(files, options);
                failures += expect(
                    compare(reference, r, kMates).empty() && r.warnings.size() == 1 &&
                    contains(r.warnings[0], "blank line at line " +
                             std::to_string(linesPerRead * 10 + 1) + " ") &&
                    contains(r.summary, skippedLine(mate, expected, laneFileName(files[mate], 0))),
                    c.name);
            }

            // Two input files: a WARNING and a count per file; three blank
            // lines at the end of file 0 (before the next file's marker) and
            // two at the end of file 1 (end of input) are skipped.
            {
                Case c{"two_files_blank_lines_" + tag, cleanTwo.mates};
                c.mates[mate] = "FILE 0\n" + joinRecords(records[mate], 0, kSplit, {{5, "\n"}}, "\n\n\n") +
                                "FILE 1\n" + joinRecords(records[mate], kSplit, kCount,
                                                         {{7, "\n"}, {8, "\n"}}, "\n\n");
                const std::vector<std::string> files = writeCase(dir, c);
                const Result r = runThreads(files, options);
                // FASTA: the other mates warn once each for their file-0 blank line.
                const size_t otherWarnings = format == '@' ? 0 : kMates - 1;
                std::vector<std::string> mine;
                for (const std::string& w : r.warnings) {
                    if (contains(w, " (mate " + std::to_string(mate + 1) + "): ")) {
                        mine.push_back(w);
                    }
                }
                failures += expect(
                    compare(referenceTwo, r, kMates).empty() &&
                    r.warnings.size() == 2 + otherWarnings && mine.size() == 2 &&
                    contains(mine[0], laneFileName(files[mate], 0) + " (mate " +
                             std::to_string(mate + 1) + "): blank line at line " +
                             std::to_string(linesPerRead * 5 + 1) + " ") &&
                    contains(mine[1], laneFileName(files[mate], 1) + " (mate " +
                             std::to_string(mate + 1) + "): blank line at line " +
                             std::to_string(linesPerRead * 7 + 1) + " ") &&
                    contains(r.summary, skippedLine(mate, 4, laneFileName(files[mate], 0))) &&
                    contains(r.summary, skippedLine(mate, 4, laneFileName(files[mate], 1))),
                    c.name);
            }
        }
    }
    std::cout << "blank lines: " << failures << " failures\n";
    return failures;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: fastx_mate_reader_harness <scratch dir>\n";
        return 2;
    }
    const std::string dir = argv[1];
    const int failures = identityPart(dir) + behaviourPart(dir) + blankLinePart(dir);
    std::cout << (failures == 0 ? "PASS" : "FAIL") << " fastx_mate_reader_harness"
              << " (one-permit pool acquires: " << g_pool.acquires << ")\n";
    return failures == 0 ? 0 : 1;
}

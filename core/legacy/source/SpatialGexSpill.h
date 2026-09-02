#ifndef STAR_SPATIAL_GEX_SPILL_H
#define STAR_SPATIAL_GEX_SPILL_H

#include "SpatialGex.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace spatial_gex {
namespace spill {

// A completed run is immutable. It is retained until all final matrices and
// their completion marker have committed successfully.
struct Run {
    std::string path;
    std::uint64_t reads = 0;
    std::uint64_t candidates = 0;
    std::uint64_t bytes = 0;
};

// Sorts a thread-local segment by (gene, raw UMI, source ordinal), writes a
// compact binary run through a .partial file, and atomically renames it.
bool writeRun(const std::string &directory, const std::string &sourceRevision,
              std::uint32_t threadIndex, std::uint64_t runIndex,
              const std::vector<ReadEvidence> &reads,
              const std::vector<CandidateEvidence> &candidates,
              Run &run, std::string &error);

class Cursor {
  public:
    static std::unique_ptr<Cursor> open(const Run &run,
                                        const std::string &sourceRevision,
                                        std::string &error);
    ~Cursor();
    Cursor(const Cursor &) = delete;
    Cursor &operator=(const Cursor &) = delete;

    // Returns one complete read/candidate bundle. False with an empty error is
    // clean EOF after validating the checksum and completion trailer.
    bool next(ReadEvidence &read, std::vector<CandidateEvidence> &candidates,
              std::string &error);

  private:
    struct Impl;
    explicit Cursor(std::unique_ptr<Impl> impl);
    std::unique_ptr<Impl> impl_;
};

} // namespace spill
} // namespace spatial_gex

#endif

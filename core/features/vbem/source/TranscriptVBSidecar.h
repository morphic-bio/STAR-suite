#ifndef STAR_SUITE_TRANSCRIPT_VB_SIDECAR_H
#define STAR_SUITE_TRANSCRIPT_VB_SIDECAR_H

#include "libem/em_types.h"
#include "libem/fld_accumulator.h"
#include "libem/gc_bias.h"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace transcript_vb_sidecar {

struct Metadata {
    std::string source_revision;
    std::string reference_id;
    std::string transcript_dictionary_sha256;
    std::string transcriptome_fasta_fingerprint;
    std::string quant_parameters_sha256;
    std::string calibration_sha256;
    std::string sample_id;
    std::string input_manifest_sha256;

    std::uint32_t transcript_count = 0;
    std::uint32_t shard_ordinal = 0;
    std::uint32_t shard_count = 1;
    std::uint32_t library_format_id = 0;
    std::uint64_t first_pair = 0;
    std::uint64_t pair_count = 0;
};

struct Evidence {
    ECTable ec_table;
    std::array<double, GCFragModel::GC_BINS> gc_counts{};
    FLDAccumulator::State fld_state;
    std::uint64_t processed_fragments = 0;
    std::uint64_t dropped_incompat = 0;
    std::uint64_t dropped_missing_mate_fields = 0;
    std::uint64_t dropped_unknown_obs_fmt = 0;
};

struct Sidecar {
    Metadata metadata;
    Evidence evidence;
};

std::string sha256Hex(const std::string& value);

bool writeAtomic(const std::string& path, const Sidecar& sidecar,
                 std::string& error);
bool read(const std::string& path, Sidecar& sidecar, std::string& error);

bool compatible(const Metadata& expected, const Metadata& observed,
                std::string& error);
bool merge(const std::vector<Sidecar>& inputs, Sidecar& output,
           std::string& error);

}  // namespace transcript_vb_sidecar

#endif

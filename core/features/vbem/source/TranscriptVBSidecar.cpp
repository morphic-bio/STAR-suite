#include "TranscriptVBSidecar.h"

#include <openssl/evp.h>
#include <openssl/sha.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace transcript_vb_sidecar {
namespace {

const char kHeaderMagic[8] = {'S', 'T', 'V', 'B', 'S', 'C', '0', '1'};
const char kFooterMagic[8] = {'S', 'T', 'V', 'B', 'E', 'N', 'D', '1'};
const std::uint16_t kMajorVersion = 1;
const std::uint16_t kMinorVersion = 0;
const std::uint8_t kLittleEndian = 1;
const std::size_t kFooterBytes = 8 + SHA256_DIGEST_LENGTH + 8;
const std::uint64_t kMaxSidecarBytes = (std::uint64_t{1} << 40);
const std::uint64_t kMaxECs = (std::uint64_t{1} << 32);
const std::uint64_t kMaxVectorValues = (std::uint64_t{1} << 34);

enum TlvType : std::uint16_t {
    kSourceRevision = 1,
    kReferenceId = 2,
    kTranscriptDictionarySha256 = 3,
    kQuantParametersSha256 = 4,
    kCalibrationSha256 = 5,
    kSampleId = 6,
    kInputManifestSha256 = 7,
    kTranscriptomeFastaFingerprint = 8
};

std::string errnoText(const std::string& prefix) {
    return prefix + ": " + std::strerror(errno);
}

void appendU8(std::vector<std::uint8_t>& out, std::uint8_t value) {
    out.push_back(value);
}

void appendU16(std::vector<std::uint8_t>& out, std::uint16_t value) {
    out.push_back(static_cast<std::uint8_t>(value & 0xff));
    out.push_back(static_cast<std::uint8_t>((value >> 8) & 0xff));
}

void appendU32(std::vector<std::uint8_t>& out, std::uint32_t value) {
    for (unsigned shift = 0; shift < 32; shift += 8) {
        out.push_back(static_cast<std::uint8_t>((value >> shift) & 0xff));
    }
}

void appendU64(std::vector<std::uint8_t>& out, std::uint64_t value) {
    for (unsigned shift = 0; shift < 64; shift += 8) {
        out.push_back(static_cast<std::uint8_t>((value >> shift) & 0xff));
    }
}

void appendDouble(std::vector<std::uint8_t>& out, double value) {
    static_assert(sizeof(double) == sizeof(std::uint64_t),
                  "TranscriptVB sidecars require IEEE-754 binary64 doubles");
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    appendU64(out, bits);
}

void appendBytes(std::vector<std::uint8_t>& out, const void* data,
                 std::size_t size) {
    const std::uint8_t* begin = static_cast<const std::uint8_t*>(data);
    out.insert(out.end(), begin, begin + size);
}

void patchU32(std::vector<std::uint8_t>& out, std::size_t offset,
              std::uint32_t value) {
    for (unsigned shift = 0; shift < 32; shift += 8) {
        out[offset++] = static_cast<std::uint8_t>((value >> shift) & 0xff);
    }
}

void appendTlv(std::vector<std::uint8_t>& out, std::uint16_t type,
               const std::string& value) {
    appendU16(out, type);
    appendU16(out, 0);
    appendU32(out, static_cast<std::uint32_t>(value.size()));
    appendBytes(out, value.data(), value.size());
}

std::string digestHex(const unsigned char* digest) {
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (std::size_t i = 0; i < SHA256_DIGEST_LENGTH; ++i) {
        out << std::setw(2) << static_cast<unsigned>(digest[i]);
    }
    return out.str();
}

struct Cursor {
    const std::vector<std::uint8_t>& bytes;
    std::size_t pos;
    std::size_t limit;
    std::string* error;

    Cursor(const std::vector<std::uint8_t>& value, std::size_t begin,
           std::size_t end, std::string& error_out)
        : bytes(value), pos(begin), limit(end), error(&error_out) {}

    bool take(std::size_t size, const std::uint8_t*& data) {
        if (pos > limit || size > limit - pos) {
            *error = "sidecar is truncated or contains an invalid section length";
            return false;
        }
        data = bytes.data() + pos;
        pos += size;
        return true;
    }

    bool u8(std::uint8_t& value) {
        const std::uint8_t* data = nullptr;
        if (!take(1, data)) return false;
        value = data[0];
        return true;
    }

    bool u16(std::uint16_t& value) {
        const std::uint8_t* data = nullptr;
        if (!take(2, data)) return false;
        value = static_cast<std::uint16_t>(data[0]) |
                (static_cast<std::uint16_t>(data[1]) << 8);
        return true;
    }

    bool u32(std::uint32_t& value) {
        const std::uint8_t* data = nullptr;
        if (!take(4, data)) return false;
        value = 0;
        for (unsigned i = 0; i < 4; ++i) {
            value |= static_cast<std::uint32_t>(data[i]) << (8 * i);
        }
        return true;
    }

    bool u64(std::uint64_t& value) {
        const std::uint8_t* data = nullptr;
        if (!take(8, data)) return false;
        value = 0;
        for (unsigned i = 0; i < 8; ++i) {
            value |= static_cast<std::uint64_t>(data[i]) << (8 * i);
        }
        return true;
    }

    bool f64(double& value) {
        std::uint64_t bits = 0;
        if (!u64(bits)) return false;
        std::memcpy(&value, &bits, sizeof(value));
        return true;
    }
};

bool validateEvidence(const Sidecar& sidecar, std::string& error) {
    const Evidence& evidence = sidecar.evidence;
    const std::uint64_t n_transcripts = sidecar.metadata.transcript_count;
    if (sidecar.metadata.shard_count == 0 ||
        sidecar.metadata.shard_ordinal >= sidecar.metadata.shard_count) {
        error = "invalid shard ordinal/count in TranscriptVB sidecar metadata";
        return false;
    }
    if (evidence.ec_table.n_transcripts != n_transcripts) {
        error = "sidecar EC transcript count disagrees with its metadata";
        return false;
    }
    if (evidence.ec_table.n_ecs != evidence.ec_table.ecs.size()) {
        error = "sidecar EC count disagrees with the EC vector";
        return false;
    }
    for (const EC& ec : evidence.ec_table.ecs) {
        if (ec.transcript_ids.empty()) {
            error = "sidecar contains an empty equivalence class";
            return false;
        }
        if (!ec.weights.empty() && ec.weights.size() != ec.transcript_ids.size()) {
            error = "sidecar EC weight and transcript arrays have different lengths";
            return false;
        }
        if (!std::isfinite(ec.count) || ec.count < 0.0) {
            error = "sidecar contains an invalid EC fragment count";
            return false;
        }
        for (std::uint32_t id : ec.transcript_ids) {
            if (id >= n_transcripts) {
                error = "sidecar contains a transcript ID outside its dictionary";
                return false;
            }
        }
        for (double weight : ec.weights) {
            if (!std::isfinite(weight) || weight < 0.0) {
                error = "sidecar contains an invalid EC weight";
                return false;
            }
        }
    }
    for (double value : evidence.gc_counts) {
        if (!std::isfinite(value) || value < 0.0) {
            error = "sidecar contains invalid GC observations";
            return false;
        }
    }
    if (evidence.fld_state.histogram.size() !=
        static_cast<std::size_t>(FLDAccumulator::MAX_FRAG_LEN + 1)) {
        error = "sidecar contains an invalid FLD histogram length";
        return false;
    }
    return true;
}

bool writeAll(int fd, const std::uint8_t* data, std::size_t size,
              std::string& error) {
    while (size > 0) {
        const ssize_t result = ::write(fd, data, size);
        if (result < 0) {
            if (errno == EINTR) continue;
            error = errnoText("could not write TranscriptVB sidecar");
            return false;
        }
        if (result == 0) {
            error = "short write while creating TranscriptVB sidecar";
            return false;
        }
        data += result;
        size -= static_cast<std::size_t>(result);
    }
    return true;
}

bool sameField(const char* name, const std::string& expected,
               const std::string& observed, std::string& error) {
    if (expected == observed) return true;
    error = std::string("incompatible TranscriptVB sidecars: ") + name +
            " differs (expected '" + expected + "', observed '" + observed + "')";
    return false;
}

}  // namespace

std::string sha256Hex(const std::string& value) {
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256(reinterpret_cast<const unsigned char*>(value.data()), value.size(), digest);
    return digestHex(digest);
}

bool writeAtomic(const std::string& path, const Sidecar& sidecar,
                 std::string& error) {
    if (path.empty() || path == "-" || path == "None") {
        error = "TranscriptVB sidecar output path is empty";
        return false;
    }
    if (!validateEvidence(sidecar, error)) return false;

    std::uint64_t transcript_id_count = 0;
    std::uint64_t signature_id_count = 0;
    std::uint64_t weight_count = 0;
    for (const EC& ec : sidecar.evidence.ec_table.ecs) {
        transcript_id_count += ec.transcript_ids.size();
        signature_id_count += ec.signature_ids.size();
        weight_count += ec.weights.size();
    }

    std::vector<std::uint8_t> header;
    appendBytes(header, kHeaderMagic, sizeof(kHeaderMagic));
    appendU16(header, kMajorVersion);
    appendU16(header, kMinorVersion);
    appendU8(header, kLittleEndian);
    appendU8(header, 0x03);  // GC and FLD sections are present.
    appendU16(header, 0);
    const std::size_t header_size_offset = header.size();
    appendU32(header, 0);
    const std::size_t payload_size_offset = header.size();
    appendU64(header, 0);
    appendU64(header, sidecar.evidence.ec_table.ecs.size());
    appendU64(header, transcript_id_count);
    appendU64(header, signature_id_count);
    appendU64(header, weight_count);
    appendU32(header, sidecar.metadata.transcript_count);
    appendU32(header, sidecar.metadata.shard_ordinal);
    appendU32(header, sidecar.metadata.shard_count);
    appendU32(header, sidecar.metadata.library_format_id);
    appendU64(header, sidecar.metadata.first_pair);
    appendU64(header, sidecar.metadata.pair_count);
    appendU64(header, sidecar.evidence.processed_fragments);
    appendU64(header, sidecar.evidence.dropped_incompat);
    appendU64(header, sidecar.evidence.dropped_missing_mate_fields);
    appendU64(header, sidecar.evidence.dropped_unknown_obs_fmt);
    appendU32(header, 8);
    appendU32(header, 0);
    appendTlv(header, kSourceRevision, sidecar.metadata.source_revision);
    appendTlv(header, kReferenceId, sidecar.metadata.reference_id);
    appendTlv(header, kTranscriptDictionarySha256,
              sidecar.metadata.transcript_dictionary_sha256);
    appendTlv(header, kQuantParametersSha256,
              sidecar.metadata.quant_parameters_sha256);
    appendTlv(header, kCalibrationSha256,
              sidecar.metadata.calibration_sha256);
    appendTlv(header, kSampleId, sidecar.metadata.sample_id);
    appendTlv(header, kInputManifestSha256,
              sidecar.metadata.input_manifest_sha256);
    appendTlv(header, kTranscriptomeFastaFingerprint,
              sidecar.metadata.transcriptome_fasta_fingerprint);
    while ((header.size() % 64) != 0) header.push_back(0);

    std::vector<std::uint8_t> payload;
    payload.reserve(sidecar.evidence.ec_table.ecs.size() * 48 + 16384);
    for (const EC& ec : sidecar.evidence.ec_table.ecs) {
        appendU32(payload, static_cast<std::uint32_t>(ec.transcript_ids.size()));
        appendU32(payload, static_cast<std::uint32_t>(ec.signature_ids.size()));
        appendU32(payload, static_cast<std::uint32_t>(ec.weights.size()));
        appendU32(payload, 0);
        appendDouble(payload, ec.count);
        for (std::uint32_t id : ec.transcript_ids) appendU32(payload, id);
        for (std::uint32_t id : ec.signature_ids) appendU32(payload, id);
        for (double weight : ec.weights) appendDouble(payload, weight);
    }
    for (double value : sidecar.evidence.gc_counts) appendDouble(payload, value);
    appendU32(payload, static_cast<std::uint32_t>(
        sidecar.evidence.fld_state.histogram.size()));
    appendU32(payload, 0);
    appendDouble(payload, sidecar.evidence.fld_state.total_mass);
    appendDouble(payload, sidecar.evidence.fld_state.weighted_sum);
    appendU64(payload, sidecar.evidence.fld_state.minimum);
    for (double value : sidecar.evidence.fld_state.histogram) {
        appendDouble(payload, value);
    }

    if (header.size() > std::numeric_limits<std::uint32_t>::max()) {
        error = "TranscriptVB sidecar metadata header is too large";
        return false;
    }
    patchU32(header, header_size_offset,
             static_cast<std::uint32_t>(header.size()));
    for (unsigned shift = 0; shift < 64; shift += 8) {
        header[payload_size_offset + shift / 8] =
            static_cast<std::uint8_t>((static_cast<std::uint64_t>(payload.size()) >> shift) & 0xff);
    }

    unsigned char digest[SHA256_DIGEST_LENGTH];
    unsigned int digest_size = 0;
    EVP_MD_CTX* context = EVP_MD_CTX_new();
    if (context == nullptr ||
        EVP_DigestInit_ex(context, EVP_sha256(), nullptr) != 1 ||
        EVP_DigestUpdate(context, header.data(), header.size()) != 1 ||
        EVP_DigestUpdate(context, payload.data(), payload.size()) != 1 ||
        EVP_DigestFinal_ex(context, digest, &digest_size) != 1 ||
        digest_size != SHA256_DIGEST_LENGTH) {
        if (context != nullptr) EVP_MD_CTX_free(context);
        error = "could not compute the TranscriptVB sidecar SHA-256 checksum";
        return false;
    }
    EVP_MD_CTX_free(context);

    std::vector<std::uint8_t> footer;
    appendBytes(footer, kFooterMagic, sizeof(kFooterMagic));
    appendBytes(footer, digest, sizeof(digest));
    appendU64(footer, header.size() + payload.size() + kFooterBytes);

    const std::string partial = path + ".partial";
    const int fd = ::open(partial.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
    if (fd < 0) {
        error = errnoText("could not create TranscriptVB sidecar " + partial);
        return false;
    }
    bool ok = writeAll(fd, header.data(), header.size(), error) &&
              writeAll(fd, payload.data(), payload.size(), error) &&
              writeAll(fd, footer.data(), footer.size(), error);
    if (ok && ::fsync(fd) != 0) {
        error = errnoText("could not fsync TranscriptVB sidecar " + partial);
        ok = false;
    }
    if (::close(fd) != 0 && ok) {
        error = errnoText("could not close TranscriptVB sidecar " + partial);
        ok = false;
    }
    if (!ok) {
        ::unlink(partial.c_str());
        return false;
    }
    if (::rename(partial.c_str(), path.c_str()) != 0) {
        error = errnoText("could not publish TranscriptVB sidecar " + path);
        ::unlink(partial.c_str());
        return false;
    }
    return true;
}

bool read(const std::string& path, Sidecar& sidecar, std::string& error) {
    struct stat info;
    if (::stat(path.c_str(), &info) != 0) {
        error = errnoText("could not stat TranscriptVB sidecar " + path);
        return false;
    }
    if (info.st_size < 0 || static_cast<std::uint64_t>(info.st_size) > kMaxSidecarBytes ||
        static_cast<std::size_t>(info.st_size) < kFooterBytes + 128) {
        error = "TranscriptVB sidecar has an invalid file size";
        return false;
    }
    std::vector<std::uint8_t> bytes(static_cast<std::size_t>(info.st_size));
    std::ifstream input(path.c_str(), std::ios::binary);
    if (!input.is_open()) {
        error = "could not open TranscriptVB sidecar " + path;
        return false;
    }
    input.read(reinterpret_cast<char*>(bytes.data()), bytes.size());
    if (!input || static_cast<std::size_t>(input.gcount()) != bytes.size()) {
        error = "could not read the complete TranscriptVB sidecar " + path;
        return false;
    }

    const std::size_t footer_offset = bytes.size() - kFooterBytes;
    if (std::memcmp(bytes.data() + footer_offset, kFooterMagic, sizeof(kFooterMagic)) != 0) {
        error = "TranscriptVB sidecar is incomplete or has an invalid footer";
        return false;
    }
    const unsigned char* recorded_digest = bytes.data() + footer_offset + 8;
    unsigned char observed_digest[SHA256_DIGEST_LENGTH];
    SHA256(bytes.data(), footer_offset, observed_digest);
    if (std::memcmp(recorded_digest, observed_digest, SHA256_DIGEST_LENGTH) != 0) {
        error = "TranscriptVB sidecar content checksum mismatch";
        return false;
    }
    std::string footer_error;
    Cursor footer(bytes, footer_offset + 8 + SHA256_DIGEST_LENGTH,
                  bytes.size(), footer_error);
    std::uint64_t recorded_size = 0;
    if (!footer.u64(recorded_size) || recorded_size != bytes.size()) {
        error = "TranscriptVB sidecar footer size does not match the file";
        return false;
    }

    Cursor cursor(bytes, 0, footer_offset, error);
    const std::uint8_t* magic = nullptr;
    if (!cursor.take(8, magic) ||
        std::memcmp(magic, kHeaderMagic, sizeof(kHeaderMagic)) != 0) {
        error = "not a TranscriptVB evidence sidecar";
        return false;
    }
    std::uint16_t major = 0, minor = 0, reserved16 = 0;
    std::uint8_t endian = 0, flags = 0;
    std::uint32_t header_bytes = 0;
    std::uint64_t payload_bytes = 0, ec_count = 0;
    std::uint64_t transcript_id_count = 0, signature_id_count = 0, weight_count = 0;
    if (!cursor.u16(major) || !cursor.u16(minor) || !cursor.u8(endian) ||
        !cursor.u8(flags) || !cursor.u16(reserved16) || !cursor.u32(header_bytes) ||
        !cursor.u64(payload_bytes) || !cursor.u64(ec_count) ||
        !cursor.u64(transcript_id_count) || !cursor.u64(signature_id_count) ||
        !cursor.u64(weight_count)) return false;
    if (major != kMajorVersion || minor > kMinorVersion || endian != kLittleEndian) {
        error = "unsupported TranscriptVB sidecar schema or byte order";
        return false;
    }
    if ((flags & 0x03) != 0x03 || header_bytes < cursor.pos ||
        header_bytes > footer_offset || payload_bytes != footer_offset - header_bytes) {
        error = "TranscriptVB sidecar header/payload boundaries are invalid";
        return false;
    }
    if (ec_count > kMaxECs || transcript_id_count > kMaxVectorValues ||
        signature_id_count > kMaxVectorValues || weight_count > kMaxVectorValues) {
        error = "TranscriptVB sidecar declares unreasonable array sizes";
        return false;
    }

    Sidecar parsed;
    std::uint32_t tlv_count = 0, reserved32 = 0;
    if (!cursor.u32(parsed.metadata.transcript_count) ||
        !cursor.u32(parsed.metadata.shard_ordinal) ||
        !cursor.u32(parsed.metadata.shard_count) ||
        !cursor.u32(parsed.metadata.library_format_id) ||
        !cursor.u64(parsed.metadata.first_pair) ||
        !cursor.u64(parsed.metadata.pair_count) ||
        !cursor.u64(parsed.evidence.processed_fragments) ||
        !cursor.u64(parsed.evidence.dropped_incompat) ||
        !cursor.u64(parsed.evidence.dropped_missing_mate_fields) ||
        !cursor.u64(parsed.evidence.dropped_unknown_obs_fmt) ||
        !cursor.u32(tlv_count) || !cursor.u32(reserved32)) return false;
    if (cursor.pos > header_bytes) {
        error = "TranscriptVB sidecar fixed metadata exceeds the declared header";
        return false;
    }
    if (tlv_count > 1024) {
        error = "TranscriptVB sidecar has an unreasonable TLV count";
        return false;
    }
    for (std::uint32_t i = 0; i < tlv_count; ++i) {
        std::uint16_t type = 0, tlv_flags = 0;
        std::uint32_t length = 0;
        if (!cursor.u16(type) || !cursor.u16(tlv_flags) || !cursor.u32(length)) return false;
        if (cursor.pos > header_bytes || length > header_bytes - cursor.pos) {
            error = "TranscriptVB sidecar TLV exceeds the metadata header";
            return false;
        }
        const std::uint8_t* value_bytes = nullptr;
        if (!cursor.take(length, value_bytes)) return false;
        const std::string value(reinterpret_cast<const char*>(value_bytes), length);
        switch (type) {
            case kSourceRevision: parsed.metadata.source_revision = value; break;
            case kReferenceId: parsed.metadata.reference_id = value; break;
            case kTranscriptDictionarySha256:
                parsed.metadata.transcript_dictionary_sha256 = value; break;
            case kQuantParametersSha256:
                parsed.metadata.quant_parameters_sha256 = value; break;
            case kCalibrationSha256:
                parsed.metadata.calibration_sha256 = value; break;
            case kSampleId: parsed.metadata.sample_id = value; break;
            case kInputManifestSha256:
                parsed.metadata.input_manifest_sha256 = value; break;
            case kTranscriptomeFastaFingerprint:
                parsed.metadata.transcriptome_fasta_fingerprint = value; break;
            default: break;
        }
    }
    cursor.pos = header_bytes;
    cursor.limit = footer_offset;

    parsed.evidence.ec_table.n_transcripts = parsed.metadata.transcript_count;
    parsed.evidence.ec_table.n_ecs = static_cast<std::size_t>(ec_count);
    parsed.evidence.ec_table.ecs.reserve(static_cast<std::size_t>(ec_count));
    std::uint64_t observed_transcript_ids = 0;
    std::uint64_t observed_signature_ids = 0;
    std::uint64_t observed_weights = 0;
    for (std::uint64_t i = 0; i < ec_count; ++i) {
        std::uint32_t n_tx = 0, n_sig = 0, n_weights = 0;
        double count = 0.0;
        if (!cursor.u32(n_tx) || !cursor.u32(n_sig) || !cursor.u32(n_weights) ||
            !cursor.u32(reserved32) || !cursor.f64(count)) return false;
        if (n_tx == 0 || n_tx > parsed.metadata.transcript_count ||
            (n_weights != 0 && n_weights != n_tx)) {
            error = "TranscriptVB sidecar contains an invalid EC row shape";
            return false;
        }
        observed_transcript_ids += n_tx;
        observed_signature_ids += n_sig;
        observed_weights += n_weights;
        if (observed_transcript_ids > transcript_id_count ||
            observed_signature_ids > signature_id_count ||
            observed_weights > weight_count) {
            error = "TranscriptVB sidecar EC arrays exceed their header totals";
            return false;
        }
        EC ec;
        ec.count = count;
        ec.transcript_ids.resize(n_tx);
        ec.signature_ids.resize(n_sig);
        ec.weights.resize(n_weights);
        for (std::uint32_t& id : ec.transcript_ids) if (!cursor.u32(id)) return false;
        for (std::uint32_t& id : ec.signature_ids) if (!cursor.u32(id)) return false;
        for (double& value : ec.weights) if (!cursor.f64(value)) return false;
        parsed.evidence.ec_table.ecs.push_back(std::move(ec));
    }
    if (observed_transcript_ids != transcript_id_count ||
        observed_signature_ids != signature_id_count || observed_weights != weight_count) {
        error = "TranscriptVB sidecar EC array totals do not match the header";
        return false;
    }
    for (double& value : parsed.evidence.gc_counts) if (!cursor.f64(value)) return false;
    std::uint32_t histogram_size = 0;
    if (!cursor.u32(histogram_size) || !cursor.u32(reserved32) ||
        !cursor.f64(parsed.evidence.fld_state.total_mass) ||
        !cursor.f64(parsed.evidence.fld_state.weighted_sum) ||
        !cursor.u64(parsed.evidence.fld_state.minimum)) return false;
    if (histogram_size != static_cast<std::uint32_t>(FLDAccumulator::MAX_FRAG_LEN + 1)) {
        error = "TranscriptVB sidecar FLD histogram has the wrong length";
        return false;
    }
    parsed.evidence.fld_state.histogram.resize(histogram_size);
    for (double& value : parsed.evidence.fld_state.histogram) {
        if (!cursor.f64(value)) return false;
    }
    if (cursor.pos != footer_offset) {
        error = "TranscriptVB sidecar payload has trailing or missing bytes";
        return false;
    }
    if (!validateEvidence(parsed, error)) return false;
    sidecar = std::move(parsed);
    return true;
}

bool compatible(const Metadata& expected, const Metadata& observed,
                std::string& error) {
    if (expected.transcript_count != observed.transcript_count) {
        error = "incompatible TranscriptVB sidecars: transcript counts differ";
        return false;
    }
    if (expected.library_format_id != observed.library_format_id) {
        error = "incompatible TranscriptVB sidecars: library formats differ";
        return false;
    }
    return sameField("source revision", expected.source_revision,
                     observed.source_revision, error) &&
           sameField("reference", expected.reference_id,
                     observed.reference_id, error) &&
           sameField("transcript dictionary", expected.transcript_dictionary_sha256,
                     observed.transcript_dictionary_sha256, error) &&
           sameField("transcriptome FASTA", expected.transcriptome_fasta_fingerprint,
                     observed.transcriptome_fasta_fingerprint, error) &&
           sameField("quantification parameters", expected.quant_parameters_sha256,
                     observed.quant_parameters_sha256, error) &&
           sameField("calibration", expected.calibration_sha256,
                     observed.calibration_sha256, error) &&
           sameField("sample", expected.sample_id, observed.sample_id, error) &&
           sameField("input manifest", expected.input_manifest_sha256,
                     observed.input_manifest_sha256, error);
}

bool merge(const std::vector<Sidecar>& inputs, Sidecar& output,
           std::string& error) {
    if (inputs.empty()) {
        error = "cannot merge an empty TranscriptVB sidecar set";
        return false;
    }
    std::vector<const Sidecar*> ordered;
    ordered.reserve(inputs.size());
    for (const Sidecar& input : inputs) ordered.push_back(&input);
    std::sort(ordered.begin(), ordered.end(), [](const Sidecar* lhs, const Sidecar* rhs) {
        return lhs->metadata.shard_ordinal < rhs->metadata.shard_ordinal;
    });
    for (std::size_t i = 0; i < ordered.size(); ++i) {
        if (!compatible(ordered[0]->metadata, ordered[i]->metadata, error)) return false;
        if (ordered[i]->metadata.shard_count !=
            ordered[0]->metadata.shard_count) {
            error = "incompatible TranscriptVB sidecars: shard counts differ";
            return false;
        }
        if (i > 0 && ordered[i - 1]->metadata.shard_ordinal ==
                         ordered[i]->metadata.shard_ordinal) {
            error = "duplicate shard ordinal in TranscriptVB sidecar set";
            return false;
        }
    }
    if (ordered.size() != ordered[0]->metadata.shard_count) {
        error = "incomplete TranscriptVB sidecar set: input count differs from declared shard count";
        return false;
    }
    for (std::size_t i = 0; i < ordered.size(); ++i) {
        if (ordered[i]->metadata.shard_ordinal != i) {
            error = "incomplete TranscriptVB sidecar set: shard ordinals are not contiguous from zero";
            return false;
        }
    }

    Sidecar result;
    result.metadata = ordered[0]->metadata;
    result.metadata.shard_ordinal = 0;
    result.metadata.shard_count = 1;
    result.metadata.first_pair = ordered.front()->metadata.first_pair;
    result.metadata.pair_count = 0;
    result.evidence.ec_table.n_transcripts = result.metadata.transcript_count;
    result.evidence.gc_counts.fill(0.0);

    std::map<std::vector<std::uint32_t>, std::size_t> ec_index;
    bool have_fld = false;
    FLDAccumulator combined_fld;
    for (const Sidecar* input : ordered) {
        result.metadata.pair_count += input->metadata.pair_count;
        result.evidence.processed_fragments += input->evidence.processed_fragments;
        result.evidence.dropped_incompat += input->evidence.dropped_incompat;
        result.evidence.dropped_missing_mate_fields +=
            input->evidence.dropped_missing_mate_fields;
        result.evidence.dropped_unknown_obs_fmt +=
            input->evidence.dropped_unknown_obs_fmt;
        for (std::size_t i = 0; i < result.evidence.gc_counts.size(); ++i) {
            result.evidence.gc_counts[i] += input->evidence.gc_counts[i];
        }
        FLDAccumulator shard_fld;
        if (!shard_fld.restore(input->evidence.fld_state, error)) return false;
        if (!have_fld) {
            if (!combined_fld.restore(input->evidence.fld_state, error)) return false;
            have_fld = true;
        } else {
            combined_fld.combine(shard_fld);
        }

        for (const EC& input_ec : input->evidence.ec_table.ecs) {
            const std::vector<std::uint32_t>& signature = input_ec.signature_ids.empty()
                ? input_ec.transcript_ids : input_ec.signature_ids;
            auto found = ec_index.find(signature);
            if (found == ec_index.end()) {
                const std::size_t index = result.evidence.ec_table.ecs.size();
                result.evidence.ec_table.ecs.push_back(input_ec);
                ec_index[signature] = index;
            } else {
                EC& target = result.evidence.ec_table.ecs[found->second];
                if (target.transcript_ids != input_ec.transcript_ids ||
                    target.weights.size() != input_ec.weights.size()) {
                    error = "EC signature collision with incompatible TranscriptVB payloads";
                    return false;
                }
                target.count += input_ec.count;
                for (std::size_t i = 0; i < target.weights.size(); ++i) {
                    target.weights[i] += input_ec.weights[i];
                }
            }
        }
    }
    result.evidence.ec_table.n_ecs = result.evidence.ec_table.ecs.size();
    result.evidence.fld_state = combined_fld.snapshot();
    if (!validateEvidence(result, error)) return false;
    output = std::move(result);
    return true;
}

}  // namespace transcript_vb_sidecar

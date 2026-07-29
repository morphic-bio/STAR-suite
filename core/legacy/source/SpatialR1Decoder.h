#ifndef STAR_SPATIAL_R1_DECODER_H
#define STAR_SPATIAL_R1_DECODER_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace spatial_r1_decoder {

struct Candidate {
    std::uint32_t coordinateIndex = 0;
    std::uint32_t auditBits = 0;
    double logSequenceLikelihood = 0.0;
};

struct Result {
    std::uint32_t rawUmi = 0;
    bool rawUmiValid = false;
    bool rawUmiHadN = false;
    bool decoderAssigned = false;
    bool barcodeHadN = false;
    bool barcodeNHashChecked = false;
    bool barcodeHadUnsupportedBase = false;
    std::uint8_t barcodeNCount = 0;
    std::uint8_t bc1Edit = 0;
    std::uint8_t bc2Edit = 0;
    std::vector<Candidate> candidates;
};

struct ExactH0Counts {
    std::vector<std::uint64_t> bc1;
    std::vector<std::uint64_t> bc2;
    std::uint64_t reads = 0;

    void reset(std::size_t bc1Size, std::size_t bc2Size);
    void add(const ExactH0Counts &other);
};

struct Config {
    std::string bc1OligosPath;
    std::string bc2OligosPath;
    std::uint32_t gridRows = 3350;
    std::uint32_t gridColumns = 3350;
    int fullStartMin = 8;
    int fullStartMax = 12;
};

class Decoder {
  public:
    explicit Decoder(const Config &config);
    ~Decoder();
    Decoder(const Decoder &) = delete;
    Decoder &operator=(const Decoder &) = delete;

    bool decode(const char *sequence, std::size_t sequenceLength,
                const char *quality, std::size_t qualityLength,
                Result &result, ExactH0Counts *h0, std::string &error) const;

    std::size_t bc1Count() const;
    std::size_t bc2Count() const;
    std::uint32_t gridRows() const;
    std::uint32_t gridColumns() const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

bool encodeRawUmi9(const char *sequence, std::size_t length,
                   std::uint32_t &packed);
std::string decodeRawUmi9(std::uint32_t packed);

std::uint32_t packAuditBits(std::uint8_t tier, std::uint8_t bc1Edit,
                            std::uint8_t bc2Edit, std::uint8_t bc1ObservedLength,
                            std::uint8_t bc2ObservedLength,
                            std::uint8_t fullStart);

} // namespace spatial_r1_decoder

#endif

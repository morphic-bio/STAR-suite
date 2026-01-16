#ifndef CELLRANGER_FORMATTER_H
#define CELLRANGER_FORMATTER_H

/**
 * CellRangerFormatter: Format FASTA and GTF files using CellRanger-style rules.
 * 
 * Matches the Perl script format-fa-gtf.pl exactly:
 * - FASTA header transformation (chr renaming)
 * - GTF biotype filtering (two-pass)
 * - Version stripping (ENSG/ENST/ENSE)
 * 
 * Source of truth:
 * /mnt/pikachu/Bulk-RNA-seq/star_salmon_aws/widgets/star_salmon_aws/10x_format_fa_gtf/Dockerfiles/format-fa-gtf.pl
 */

#include <string>
#include <set>
#include <unordered_set>

namespace CellRangerFormatter {

// Configuration struct
struct Config {
    std::string inputFastaPath;    // Input FASTA file
    std::string outputFastaPath;   // Output FASTA file
    std::string inputGtfPath;      // Input GTF file
    std::string outputGtfPath;      // Output GTF file
    
    Config() {}
};

// Result struct
struct Result {
    bool success;
    std::string errorMessage;
    
    Result() : success(false) {}
};

// Main entry point: format both FASTA and GTF
Result format(const Config& config);

// Individual formatters (exposed for testing)
Result formatFasta(const std::string& inputPath, const std::string& outputPath);
Result formatGtf(const std::string& inputPath, const std::string& outputPath);

// Download reference files
// Uses FTP/HTTP with optional cksum verification for integrity
// expectedCksum: cksum value to verify (0 if not available)
// expectedSize: file size to verify (0 if not available)
// allowUntrusted: if true and URL is not in trusted table, proceed without cksum verification (with warning)
// cacheDir: directory for cksum cache file (empty string to disable cache)
// autoCksumUpdate: if true, attempt to auto-fill missing cksum from CHECKSUMS file for trusted URLs
bool downloadReference(const std::string& url, const std::string& outputPath, 
                       uint32_t expectedCksum, uint64_t expectedSize, bool allowUntrusted,
                       const std::string& cacheDir, bool autoCksumUpdate, bool replaceUnverifiableFiles,
                       std::string& errorMsg);

// Check if URL is in trusted Ensembl FTP URL table
bool isTrustedUrl(const std::string& url);

// Check if URL has a known cksum entry (from embedded map or cache file)
bool hasKnownCksum(const std::string& url);

// Check if URL is trusted OR has known cksum
bool isTrustedOrHasKnownCksum(const std::string& url);

// Load cksum cache file from disk
// Returns true on success, false on error
bool loadCksumCache(const std::string& cacheDir);

// Save cksum entry to cache file
// allowOverwrite: if true, update existing entry; if false, skip if already exists
// Returns true on success, false on error
bool saveCksumToCache(const std::string& cacheDir, const std::string& url, uint32_t crc, uint64_t size, bool allowOverwrite = false);

// Compute POSIX cksum of a file (CRC-32 XOR file size)
// Returns true on success, false on error
// On success, crc and size are filled with computed values
bool computeCksumFile(const std::string& filePath, uint32_t& crc, uint64_t& size);

// Auto-fill cksum from CHECKSUMS file
// Detects provider (Ensembl/GENCODE) from URL, downloads CHECKSUMS file, extracts cksum for target filename
// Writes to cache and updates runtime cache
// Returns true on success, false on error
// On success, crc and size are filled with extracted values
bool autoFillCksumFromChecksums(const std::string& url, const std::string& cacheDir, 
                                  uint32_t& crc, uint64_t& size, std::string& errorMsg);

// Get decompressed cksum for URL (from cache)
// Returns true if found and valid, false otherwise
bool getDecompressedCksumForUrl(const std::string& url, uint32_t& crc, uint64_t& size);

// Save decompressed cksum entry to cache file
// Returns true on success, false on error
bool saveDecompressedCksumToCache(const std::string& cacheDir, const std::string& url, uint32_t crc, uint64_t size);

} // namespace CellRangerFormatter

#endif // CELLRANGER_FORMATTER_H


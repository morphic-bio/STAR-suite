#ifndef LIBTXIMPORT_H
#define LIBTXIMPORT_H

#include <string>
#include <vector>
#include <unordered_map>

namespace tximport {

/// Transcript-level record from quant.sf
struct TxRecord {
    std::string name;        // transcript ID
    double length;           // raw length (bp)
    double eff_len;          // effective length
    double tpm;              // TPM abundance
    double num_reads;        // estimated counts (fractional)
};

/// Gene-level summary output
struct GeneSummary {
    std::string name;        // gene ID
    double length;           // gene length (avg of tx lengths, for fallback)
    double eff_len;          // TPM-weighted avg effective length
    double tpm;              // gene TPM (sum of tx TPMs)
    double num_reads;        // lengthScaledTPM counts
};

/// Summarization mode
enum class CountsFromAbundance {
    No,              // raw counts (sum)
    ScaledTPM,       // TPM-proportional, preserving total
    LengthScaledTPM  // TPM × length, preserving total (default for Salmon)
};

/// Which length to use
enum class LengthMode {
    Effective,  // Use EffectiveLength (default for Salmon)
    Raw         // Use Length
};

/// Statistics from summarization
struct SummarizeStats {
    size_t tx_total;           // total transcripts in input
    size_t tx_mapped;          // transcripts found in tx2gene
    size_t tx_missing;         // transcripts not in tx2gene (filtered)
    size_t genes_output;       // genes in output
    size_t genes_zero_tpm;     // genes with zero TPM (used fallback length)
};

/**
 * Summarize transcript-level quantification to gene-level.
 *
 * Implements tximport's summarizeToGene with lengthScaledTPM mode.
 *
 * @param transcripts   Vector of transcript records from quant.sf
 * @param tx2gene       Map: transcript_id -> gene_id
 * @param mode          Count transformation mode (default: LengthScaledTPM)
 * @param length_mode   Which length to use (default: Effective)
 * @param gene_order    Output: ordered list of gene IDs (first-seen order from tx2gene)
 * @param stats         Output: summarization statistics
 * @return              Vector of gene summaries in gene_order order
 */
std::vector<GeneSummary> summarize_to_gene(
    const std::vector<TxRecord>& transcripts,
    const std::unordered_map<std::string, std::string>& tx2gene,
    CountsFromAbundance mode = CountsFromAbundance::LengthScaledTPM,
    LengthMode length_mode = LengthMode::Effective,
    std::vector<std::string>* gene_order = nullptr,
    SummarizeStats* stats = nullptr
);

/**
 * Parse quant.sf file into transcript records.
 *
 * @param path    Path to quant.sf
 * @param records Output: vector of transcript records
 * @param error   Output: error message if failed
 * @return        true on success
 */
bool parse_quant_sf(
    const std::string& path,
    std::vector<TxRecord>& records,
    std::string& error
);

/**
 * Parse tx2gene file (TSV: transcript_id<TAB>gene_id).
 *
 * @param path    Path to tx2gene.tsv
 * @param tx2gene Output: map from transcript to gene
 * @param gene_order Output: ordered list of genes (first occurrence)
 * @param error   Output: error message if failed
 * @return        true on success
 */
bool parse_tx2gene(
    const std::string& path,
    std::unordered_map<std::string, std::string>& tx2gene,
    std::vector<std::string>& gene_order,
    std::string& error
);

/**
 * Write gene summaries to file in quant.sf format.
 *
 * @param path      Output path
 * @param summaries Gene summaries to write
 * @param error     Output: error message if failed
 * @param precision Decimal precision for floating point (default: 6)
 * @return          true on success
 */
bool write_gene_sf(
    const std::string& path,
    const std::vector<GeneSummary>& summaries,
    std::string& error,
    int precision = 6
);

} // namespace tximport

#endif // LIBTXIMPORT_H


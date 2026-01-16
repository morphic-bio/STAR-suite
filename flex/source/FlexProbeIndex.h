#ifndef FLEX_PROBE_INDEX_H
#define FLEX_PROBE_INDEX_H

/**
 * FlexProbeIndex: Filter 50bp gene probes to match a target GTF and build hybrid reference.
 * 
 * Filtering rules:
 *   1. 50bp A/C/G/T only (fail on invalid)
 *   2. Skip DEPRECATED probes (case-insensitive)
 *   3. Keep only probes whose gene_id exists in the target GTF
 *   4. Stable sort by gene_id then probe_id
 * 
 * Outputs:
 *   - filtered_probe_set.csv        Filtered probes (sorted)
 *   - probes_only.fa                Probe-only FASTA
 *   - probes_only.gtf               Probe-only GTF (no header)
 *   - genome.filtered.fa            Hybrid FASTA (base + probes)
 *   - genes.filtered.gtf            Hybrid GTF (base + probes)
 *   - probe_genes_exons.bed         BED file for probes
 *   - probe_list.txt                Unique gene IDs
 *   - metadata/reference_manifest.json
 */

#include <string>
#include <vector>
#include <set>
#include <unordered_set>

namespace FlexProbeIndex {

// Configuration struct for probe filtering
struct Config {
    std::string probeCSVPath;       // Input probe CSV (required)
    std::string gtfPath;            // Target GTF for gene matching (required, .gz supported)
    std::string baseFastaPath;      // Base genome FASTA (required)
    std::string outputDir;          // Output directory (required)
    uint32_t enforceProbeLength;    // Expected probe length (default: 50)
    bool removeDeprecated;          // Remove deprecated gene IDs from probe lists (default: false)
    
    Config() : enforceProbeLength(50), removeDeprecated(false) {}
};

// Single probe entry
struct ProbeEntry {
    std::string gene_id;
    std::string probe_seq;
    std::string probe_id;
    std::string included;
    std::string region;
    std::string gene_name;
    std::string original_line;  // Full CSV line for output
    
    // For stable sorting
    bool operator<(const ProbeEntry& other) const {
        if (gene_id != other.gene_id) return gene_id < other.gene_id;
        return probe_id < other.probe_id;
    }
};

// Statistics from filtering
struct FilterStats {
    uint32_t totalInput;
    uint32_t droppedDeprecated;
    uint32_t droppedIncludedFalse;
    uint32_t droppedNoMatch;
    uint32_t droppedInvalidSeq;
    uint32_t totalOutput;
    uint32_t uniqueGenes;
    
    FilterStats() : totalInput(0), droppedDeprecated(0), droppedIncludedFalse(0), droppedNoMatch(0),
                    droppedInvalidSeq(0), totalOutput(0), uniqueGenes(0) {}
};

// Result struct
struct Result {
    bool success;
    std::string errorMessage;
    FilterStats stats;
    
    Result() : success(false) {}
};

// Main entry point
Result run(const Config& config);

// Helper functions (exposed for testing)
bool parseGTFGeneIds(const std::string& gtfPath, std::unordered_set<std::string>& geneIds, std::string& errorMsg);
bool parseFastaContigNames(const std::string& fastaPath, std::unordered_set<std::string>& contigNames, std::string& errorMsg);
bool parseProbeCSV(const std::string& csvPath, const std::unordered_set<std::string>& gtfGeneIds,
                   uint32_t enforceLength, std::vector<ProbeEntry>& probes,
                   std::vector<std::string>& headerLines, std::string& headerRow,
                   FilterStats& stats, std::string& errorMsg);
std::string computeSHA256(const std::string& filePath);

} // namespace FlexProbeIndex

#endif // FLEX_PROBE_INDEX_H


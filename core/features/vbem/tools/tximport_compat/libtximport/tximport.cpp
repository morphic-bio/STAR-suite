#include "tximport.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <unordered_set>

namespace tximport {

// Internal helper: accumulator for gene-level aggregation
struct GeneAccum {
    double tpm_sum = 0.0;           // Σ TPM
    double counts_sum = 0.0;        // Σ NumReads
    double weighted_len_sum = 0.0;  // Σ (TPM × EffLen)
    double raw_len_sum = 0.0;       // Σ Length (for output)
    double selected_len_sum = 0.0;  // Σ selected length (eff or raw, for fallback)
    int tx_count = 0;               // number of transcripts
};

std::vector<GeneSummary> summarize_to_gene(
    const std::vector<TxRecord>& transcripts,
    const std::unordered_map<std::string, std::string>& tx2gene,
    CountsFromAbundance mode,
    LengthMode length_mode,
    std::vector<std::string>* gene_order_out,
    SummarizeStats* stats_out
) {
    // Track statistics
    SummarizeStats stats = {0, 0, 0, 0, 0};
    stats.tx_total = transcripts.size();

    // Accumulate by gene, tracking first-seen order
    std::unordered_map<std::string, GeneAccum> gene_accum;
    std::vector<std::string> gene_order;  // first-seen order
    std::unordered_set<std::string> seen_genes;

    for (const auto& tx : transcripts) {
        auto it = tx2gene.find(tx.name);
        if (it == tx2gene.end()) {
            stats.tx_missing++;
            continue;
        }
        stats.tx_mapped++;

        const std::string& gene_id = it->second;
        
        // Track gene order (first occurrence)
        if (seen_genes.find(gene_id) == seen_genes.end()) {
            seen_genes.insert(gene_id);
            gene_order.push_back(gene_id);
        }

        // Choose which length to use
        double len = (length_mode == LengthMode::Effective) ? tx.eff_len : tx.length;
        
        // Accumulate
        GeneAccum& acc = gene_accum[gene_id];
        acc.tpm_sum += tx.tpm;
        acc.counts_sum += tx.num_reads;
        acc.weighted_len_sum += tx.tpm * len;
        acc.raw_len_sum += tx.length;        // always use raw for output
        acc.selected_len_sum += len;         // selected length for fallback
        acc.tx_count++;
    }

    stats.genes_output = gene_order.size();

    // Sort gene order alphabetically with case-insensitive comparison
    // (matches R's locale-aware rowsum() behavior)
    std::sort(gene_order.begin(), gene_order.end(), [](const std::string& a, const std::string& b) {
        // Case-insensitive comparison
        size_t len = std::min(a.size(), b.size());
        for (size_t i = 0; i < len; ++i) {
            char ca = std::tolower(static_cast<unsigned char>(a[i]));
            char cb = std::tolower(static_cast<unsigned char>(b[i]));
            if (ca != cb) return ca < cb;
        }
        return a.size() < b.size();
    });

    // Build gene summaries
    std::vector<GeneSummary> summaries;
    summaries.reserve(gene_order.size());

    double total_raw_counts = 0.0;
    double total_scaled = 0.0;

    for (const auto& gene_id : gene_order) {
        const GeneAccum& acc = gene_accum[gene_id];
        
        GeneSummary gs;
        gs.name = gene_id;
        gs.tpm = acc.tpm_sum;

        // Compute gene effective length (TPM-weighted average)
        if (acc.tpm_sum > 0) {
            gs.eff_len = acc.weighted_len_sum / acc.tpm_sum;
        } else {
            // Fallback: unweighted average of selected lengths (matches tximport)
            gs.eff_len = acc.selected_len_sum / acc.tx_count;
            stats.genes_zero_tpm++;
        }

        // Length column: same as EffectiveLength (matches tximport behavior)
        // tximport's output$length matrix contains TPM-weighted effective length
        gs.length = gs.eff_len;

        // Initial counts (will be transformed below)
        gs.num_reads = acc.counts_sum;

        total_raw_counts += acc.counts_sum;

        summaries.push_back(gs);
    }

    // Apply counts transformation
    if (mode == CountsFromAbundance::LengthScaledTPM) {
        // lengthScaledTPM: counts proportional to TPM × length, scaled to preserve total
        // newCounts = TPM × avgLength
        // finalCounts = newCounts × (totalRaw / totalNew)
        
        for (auto& gs : summaries) {
            gs.num_reads = gs.tpm * gs.eff_len;  // length-scaled
            total_scaled += gs.num_reads;
        }

        if (total_scaled > 0) {
            double scale_factor = total_raw_counts / total_scaled;
            for (auto& gs : summaries) {
                gs.num_reads *= scale_factor;
            }
        }
    } else if (mode == CountsFromAbundance::ScaledTPM) {
        // scaledTPM: counts proportional to TPM, scaled to preserve total
        total_scaled = 0.0;
        for (auto& gs : summaries) {
            gs.num_reads = gs.tpm;  // TPM only
            total_scaled += gs.num_reads;
        }

        if (total_scaled > 0) {
            double scale_factor = total_raw_counts / total_scaled;
            for (auto& gs : summaries) {
                gs.num_reads *= scale_factor;
            }
        }
    }
    // CountsFromAbundance::No - keep raw counts (already set)

    // Output gene order if requested
    if (gene_order_out) {
        *gene_order_out = gene_order;
    }

    // Output stats if requested
    if (stats_out) {
        *stats_out = stats;
    }

    return summaries;
}

bool parse_quant_sf(
    const std::string& path,
    std::vector<TxRecord>& records,
    std::string& error
) {
    std::ifstream file(path);
    if (!file.is_open()) {
        error = "Cannot open file: " + path;
        return false;
    }

    std::string line;
    int line_num = 0;

    // Read header
    if (!std::getline(file, line)) {
        error = "Empty file: " + path;
        return false;
    }
    line_num++;

    // Validate header columns
    std::istringstream hss(line);
    std::string col;
    std::vector<std::string> headers;
    while (hss >> col) {
        headers.push_back(col);
    }

    // Find column indices
    int name_idx = -1, len_idx = -1, eff_idx = -1, tpm_idx = -1, reads_idx = -1;
    for (size_t i = 0; i < headers.size(); i++) {
        if (headers[i] == "Name") name_idx = i;
        else if (headers[i] == "Length") len_idx = i;
        else if (headers[i] == "EffectiveLength") eff_idx = i;
        else if (headers[i] == "TPM") tpm_idx = i;
        else if (headers[i] == "NumReads") reads_idx = i;
    }

    if (name_idx < 0 || len_idx < 0 || eff_idx < 0 || tpm_idx < 0 || reads_idx < 0) {
        error = "Missing required columns in " + path + 
                ". Expected: Name, Length, EffectiveLength, TPM, NumReads";
        return false;
    }

    records.clear();

    // Read data lines
    while (std::getline(file, line)) {
        line_num++;
        if (line.empty()) continue;

        std::istringstream ss(line);
        std::vector<std::string> fields;
        std::string field;
        while (ss >> field) {
            fields.push_back(field);
        }

        if (fields.size() < headers.size()) {
            error = "Line " + std::to_string(line_num) + ": insufficient fields";
            return false;
        }

        TxRecord rec;
        rec.name = fields[name_idx];
        
        try {
            rec.length = std::stod(fields[len_idx]);
            rec.eff_len = std::stod(fields[eff_idx]);
            rec.tpm = std::stod(fields[tpm_idx]);
            rec.num_reads = std::stod(fields[reads_idx]);
        } catch (const std::exception& e) {
            error = "Line " + std::to_string(line_num) + ": invalid number format";
            return false;
        }

        records.push_back(rec);
    }

    return true;
}

bool parse_tx2gene(
    const std::string& path,
    std::unordered_map<std::string, std::string>& tx2gene,
    std::vector<std::string>& gene_order,
    std::string& error
) {
    std::ifstream file(path);
    if (!file.is_open()) {
        error = "Cannot open file: " + path;
        return false;
    }

    tx2gene.clear();
    gene_order.clear();
    std::unordered_set<std::string> seen_genes;

    std::string line;
    int line_num = 0;

    while (std::getline(file, line)) {
        line_num++;
        if (line.empty()) continue;

        // Skip comment lines
        if (line[0] == '#') continue;

        // Parse tab-separated: transcript_id<TAB>gene_id
        size_t tab_pos = line.find('\t');
        if (tab_pos == std::string::npos) {
            error = "Line " + std::to_string(line_num) + ": missing tab separator";
            return false;
        }

        std::string tx_id = line.substr(0, tab_pos);
        std::string gene_id = line.substr(tab_pos + 1);

        // Trim trailing whitespace from gene_id
        while (!gene_id.empty() && (gene_id.back() == '\r' || gene_id.back() == ' ' || gene_id.back() == '\t')) {
            gene_id.pop_back();
        }

        if (tx_id.empty() || gene_id.empty()) {
            error = "Line " + std::to_string(line_num) + ": empty transcript or gene ID";
            return false;
        }

        tx2gene[tx_id] = gene_id;

        // Track gene order (first occurrence)
        if (seen_genes.find(gene_id) == seen_genes.end()) {
            seen_genes.insert(gene_id);
            gene_order.push_back(gene_id);
        }
    }

    if (tx2gene.empty()) {
        error = "No mappings found in: " + path;
        return false;
    }

    return true;
}

bool write_gene_sf(
    const std::string& path,
    const std::vector<GeneSummary>& summaries,
    std::string& error,
    int precision
) {
    std::ofstream file(path);
    if (!file.is_open()) {
        error = "Cannot create file: " + path;
        return false;
    }

    // Write header
    file << "Name\tLength\tEffectiveLength\tTPM\tNumReads\n";

    // Write data
    file << std::fixed << std::setprecision(precision);
    for (const auto& gs : summaries) {
        file << gs.name << "\t"
             << std::setprecision(3) << gs.length << "\t"
             << std::setprecision(3) << gs.eff_len << "\t"
             << std::setprecision(precision) << gs.tpm << "\t"
             << std::setprecision(3) << gs.num_reads << "\n";
    }

    if (!file.good()) {
        error = "Error writing to: " + path;
        return false;
    }

    return true;
}

} // namespace tximport


/**
 * tximport_compat - CLI tool for tximport-compatible gene-level summarization
 *
 * Usage:
 *   tximport_compat --quant quant.sf --tx2gene tx2gene.tsv --output quant.genes.sf
 *
 * Options:
 *   --quant FILE       Input transcript-level quant.sf
 *   --tx2gene FILE     Transcript-to-gene mapping (TSV: tx_id<TAB>gene_id)
 *   --output FILE      Output gene-level summary
 *   --mode MODE        Count mode: lengthScaledTPM (default), scaledTPM, no
 *   --length-mode MODE Length: effective (default), raw
 *   --sort             Sort output by gene ID (default: first-seen order)
 *   --stats            Print summarization statistics
 *   -h, --help         Show help
 */

#include <iostream>
#include <string>
#include <algorithm>
#include "libtximport/tximport.h"

using namespace tximport;

void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " --quant FILE --tx2gene FILE --output FILE [options]\n"
              << "\n"
              << "tximport-compatible gene-level summarization\n"
              << "\n"
              << "Required:\n"
              << "  --quant FILE       Input transcript-level quant.sf\n"
              << "  --tx2gene FILE     Transcript-to-gene mapping (tx_id<TAB>gene_id)\n"
              << "  --output FILE      Output gene-level summary file\n"
              << "\n"
              << "Options:\n"
              << "  --mode MODE        Count transformation mode (default: lengthScaledTPM)\n"
              << "                     lengthScaledTPM - TPM × length, normalized (tximport default)\n"
              << "                     scaledTPM       - TPM only, normalized\n"
              << "                     no              - raw counts (sum of transcripts)\n"
              << "  --length-mode MODE Which length to use (default: effective)\n"
              << "                     effective - EffectiveLength column\n"
              << "                     raw       - Length column\n"
              << "  --sort             Sort output by gene ID (default: first-seen order)\n"
              << "  --stats            Print summarization statistics to stderr\n"
              << "  -h, --help         Show this help message\n"
              << "\n"
              << "Output format matches quant.sf: Name, Length, EffectiveLength, TPM, NumReads\n";
}

int main(int argc, char* argv[]) {
    std::string quant_path;
    std::string tx2gene_path;
    std::string output_path;
    CountsFromAbundance mode = CountsFromAbundance::LengthScaledTPM;
    LengthMode length_mode = LengthMode::Effective;
    bool sort_output = false;
    bool print_stats = false;

    // Parse arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--quant" && i + 1 < argc) {
            quant_path = argv[++i];
        } else if (arg == "--tx2gene" && i + 1 < argc) {
            tx2gene_path = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            output_path = argv[++i];
        } else if (arg == "--mode" && i + 1 < argc) {
            std::string m = argv[++i];
            if (m == "lengthScaledTPM") {
                mode = CountsFromAbundance::LengthScaledTPM;
            } else if (m == "scaledTPM") {
                mode = CountsFromAbundance::ScaledTPM;
            } else if (m == "no") {
                mode = CountsFromAbundance::No;
            } else {
                std::cerr << "Error: Unknown mode '" << m << "'\n";
                return 1;
            }
        } else if (arg == "--length-mode" && i + 1 < argc) {
            std::string m = argv[++i];
            if (m == "effective") {
                length_mode = LengthMode::Effective;
            } else if (m == "raw") {
                length_mode = LengthMode::Raw;
            } else {
                std::cerr << "Error: Unknown length-mode '" << m << "'\n";
                return 1;
            }
        } else if (arg == "--sort") {
            sort_output = true;
        } else if (arg == "--stats") {
            print_stats = true;
        } else {
            std::cerr << "Error: Unknown argument '" << arg << "'\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    // Validate required arguments
    if (quant_path.empty() || tx2gene_path.empty() || output_path.empty()) {
        std::cerr << "Error: --quant, --tx2gene, and --output are required\n";
        print_usage(argv[0]);
        return 1;
    }

    std::string error;

    // Parse tx2gene
    std::unordered_map<std::string, std::string> tx2gene;
    std::vector<std::string> tx2gene_order;
    if (!parse_tx2gene(tx2gene_path, tx2gene, tx2gene_order, error)) {
        std::cerr << "Error parsing tx2gene: " << error << "\n";
        return 1;
    }

    // Parse quant.sf
    std::vector<TxRecord> transcripts;
    if (!parse_quant_sf(quant_path, transcripts, error)) {
        std::cerr << "Error parsing quant.sf: " << error << "\n";
        return 1;
    }

    // Summarize to gene level
    std::vector<std::string> gene_order;
    SummarizeStats stats;
    auto summaries = summarize_to_gene(transcripts, tx2gene, mode, length_mode,
                                       &gene_order, &stats);

    // Sort if requested
    if (sort_output) {
        std::sort(summaries.begin(), summaries.end(),
                  [](const GeneSummary& a, const GeneSummary& b) {
                      return a.name < b.name;
                  });
    }

    // Write output
    if (!write_gene_sf(output_path, summaries, error)) {
        std::cerr << "Error writing output: " << error << "\n";
        return 1;
    }

    // Print stats if requested
    if (print_stats) {
        std::cerr << "=== Summarization Statistics ===\n";
        std::cerr << "Transcripts total:    " << stats.tx_total << "\n";
        std::cerr << "Transcripts mapped:   " << stats.tx_mapped << "\n";
        std::cerr << "Transcripts missing:  " << stats.tx_missing << "\n";
        std::cerr << "Genes output:         " << stats.genes_output << "\n";
        std::cerr << "Genes with zero TPM:  " << stats.genes_zero_tpm << "\n";
    }

    return 0;
}



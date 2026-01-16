/*
 * compute_gc_bias - Compute GC-corrected effective lengths
 * 
 * Combines expected GC, observed GC, and FLD to compute effective lengths
 */

#include "../../source/libem/gc_bias.h"
#include "../../source/libem/effective_length.h"
#include "../../source/libem/alignment_model.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <iomanip>

void print_usage(const char* progname) {
    std::cerr << "Usage: " << progname << " [options]\n\n";
    std::cerr << "Required:\n";
    std::cerr << "  --expected-gc <file>     Expected GC distribution (TSV: gc_bin <tab> prob)\n";
    std::cerr << "  --observed-gc <file>     Observed GC distribution (TSV: gc_bin <tab> prob)\n";
    std::cerr << "  --fld <file>             Fragment length distribution (TSV: len <tab> count <tab> prob)\n";
    std::cerr << "  --transcriptome <file>   Transcriptome FASTA file\n";
    std::cerr << "  --output <file>          Output effective lengths (TSV: transcript_id <tab> eff_length)\n\n";
    std::cerr << "Options:\n";
    std::cerr << "  --max-bias-ratio <float> Maximum bias ratio (default: 1000.0)\n";
    std::cerr << "  -h, --help               Show this help\n";
}

int main(int argc, char* argv[]) {
    std::string expected_gc_file;
    std::string observed_gc_file;
    std::string fld_file;
    std::string transcriptome_file;
    std::string output_file;
    double max_bias_ratio = 1000.0;
    
    // Parse command line
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--expected-gc") == 0 && i + 1 < argc) {
            expected_gc_file = argv[++i];
        } else if (strcmp(argv[i], "--observed-gc") == 0 && i + 1 < argc) {
            observed_gc_file = argv[++i];
        } else if (strcmp(argv[i], "--fld") == 0 && i + 1 < argc) {
            fld_file = argv[++i];
        } else if (strcmp(argv[i], "--transcriptome") == 0 && i + 1 < argc) {
            transcriptome_file = argv[++i];
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "--max-bias-ratio") == 0 && i + 1 < argc) {
            max_bias_ratio = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
    }
    
    // Validate required arguments
    if (expected_gc_file.empty() || observed_gc_file.empty() || 
        fld_file.empty() || transcriptome_file.empty() || output_file.empty()) {
        std::cerr << "Error: All required arguments must be provided\n";
        print_usage(argv[0]);
        return 1;
    }
    
    // Load expected GC
    GCFragModel expected_gc;
    std::cerr << "Loading expected GC from: " << expected_gc_file << "\n";
    if (!expected_gc.loadFromFile(expected_gc_file)) {
        std::cerr << "Error: Failed to load expected GC\n";
        return 1;
    }
    
    // Load observed GC
    GCFragModel observed_gc;
    std::cerr << "Loading observed GC from: " << observed_gc_file << "\n";
    if (!observed_gc.loadFromFile(observed_gc_file)) {
        std::cerr << "Error: Failed to load observed GC\n";
        return 1;
    }
    
    // Compute bias ratio
    std::cerr << "Computing GC bias ratio...\n";
    std::vector<double> gc_bias = observed_gc.computeBiasRatio(expected_gc, max_bias_ratio);
    
    // Load transcriptome
    libem::Transcriptome transcriptome;
    std::cerr << "Loading transcriptome from: " << transcriptome_file << "\n";
    if (!transcriptome.loadFromFasta(transcriptome_file)) {
        std::cerr << "Error: Failed to load transcriptome\n";
        return 1;
    }
    
    // Load FLD and setup effective length calculator
    EffectiveLengthCalculator calc;
    std::cerr << "Loading FLD from: " << fld_file << "\n";
    if (!calc.loadFLD(fld_file)) {
        std::cerr << "Error: Failed to load FLD\n";
        return 1;
    }
    
    calc.loadGCBias(gc_bias);
    
    // Compute effective lengths for all transcripts
    std::cerr << "Computing effective lengths for " << transcriptome.size() << " transcripts...\n";
    std::vector<double> raw_lengths;
    raw_lengths.reserve(transcriptome.size());
    
    for (size_t i = 0; i < transcriptome.size(); ++i) {
        const libem::TranscriptSequence* txp = transcriptome.getTranscript(static_cast<uint32_t>(i));
        if (txp) {
            raw_lengths.push_back(static_cast<double>(txp->length()));
        } else {
            raw_lengths.push_back(0.0);
        }
    }
    
    std::vector<double> eff_lengths = calc.computeAllEffectiveLengths(transcriptome, raw_lengths);
    
    // Write output
    std::cerr << "Writing effective lengths to: " << output_file << "\n";
    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open output file: " << output_file << "\n";
        return 1;
    }
    
    out << std::setprecision(17);
    for (size_t i = 0; i < eff_lengths.size(); ++i) {
        out << i << "\t" << eff_lengths[i] << "\n";
    }
    
    std::cerr << "Done. Computed effective lengths for " << eff_lengths.size() << " transcripts.\n";
    
    return 0;
}

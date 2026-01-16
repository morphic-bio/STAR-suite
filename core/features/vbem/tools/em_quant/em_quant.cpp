#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>
#include "../../source/libem/em_types.h"
#include "../../source/libem/ec_loader.h"
#include "../../source/libem/em_engine.h"
#include "../../source/libem/vb_engine.h"

void print_usage(const char* prog_name) {
    std::cerr << "Usage: " << prog_name << " [options] -e <eq_classes.txt> -o <output.tsv>\n\n";
    std::cerr << "Options:\n";
    std::cerr << "  -e, --ec <file>        Equivalence classes file (required)\n";
    std::cerr << "  -l, --lengths <file>   Transcript lengths file (quant.sf format, optional)\n";
    std::cerr << "  --effective-lengths <file> GC-corrected effective lengths (TSV: id <tab> eff_len, optional)\n";
    std::cerr << "  -o, --output <file>    Output TSV file (required)\n";
    std::cerr << "  --vb                   Use Variational Bayes instead of EM\n";
    std::cerr << "  --vb-prior <float>     Dirichlet prior concentration (default: 0.01)\n";
    std::cerr << "  --per-nucleotide-prior Use per-nucleotide prior (prior * eff_length)\n";
    std::cerr << "                         Default: per-transcript prior (Salmon default)\n";
    std::cerr << "  --max-iters <int>      Maximum iterations (default: 200 for VB, 1000 for EM)\n";
    std::cerr << "  --tolerance <float>    Convergence tolerance (default: 1e-8 for VB, 1e-6 for EM)\n";
    std::cerr << "  --zero-threshold <float> Zero out NumReads below threshold (default: 0.001)\n";
    std::cerr << "                         Set to 0 to disable zeroing (for exact EM/VB)\n";
    std::cerr << "  --init-by-length       Initialize abundances by transcript length\n";
    std::cerr << "  --uniform-init         Use uniform initialization for VB (alpha = prior only)\n";
    std::cerr << "                         Default: prior + unique_counts (backward compatible)\n";
    std::cerr << "                         Salmon uses txp.projectedCounts + uniform mixing\n";
    std::cerr << "  --threads <int>        Number of threads (default: OMP_NUM_THREADS)\n";
    std::cerr << "  --debug-trace <file>   Enable debug tracing for selected transcripts\n";
    std::cerr << "                         Use EM_QUANT_DEBUG_TXPS env var (comma-separated IDs)\n";
    std::cerr << "  -v, --verbose          Verbose output\n";
    std::cerr << "  -h, --help             Show this help message\n";
}

int main(int argc, char* argv[]) {
    std::string ec_file;
    std::string lengths_file;
    std::string effective_lengths_file;
    std::string output_file;
    EMParams params;
    bool verbose = false;
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--ec") == 0) {
            if (i + 1 < argc) {
                ec_file = argv[++i];
            } else {
                std::cerr << "Error: -e/--ec requires a file argument\n";
                return 1;
            }
        } else if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--lengths") == 0) {
            if (i + 1 < argc) {
                lengths_file = argv[++i];
            } else {
                std::cerr << "Error: -l/--lengths requires a file argument\n";
                return 1;
            }
        } else if (strcmp(argv[i], "--effective-lengths") == 0) {
            if (i + 1 < argc) {
                effective_lengths_file = argv[++i];
            } else {
                std::cerr << "Error: --effective-lengths requires a file argument\n";
                return 1;
            }
        } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            if (i + 1 < argc) {
                output_file = argv[++i];
            } else {
                std::cerr << "Error: -o/--output requires a file argument\n";
                return 1;
            }
        } else if (strcmp(argv[i], "--vb") == 0) {
            params.use_vb = true;
        } else if (strcmp(argv[i], "--vb-prior") == 0) {
            if (i + 1 < argc) {
                params.vb_prior = std::stod(argv[++i]);
            } else {
                std::cerr << "Error: --vb-prior requires a float argument\n";
                return 1;
            }
        } else if (strcmp(argv[i], "--per-nucleotide-prior") == 0) {
            params.per_transcript_prior = false;
        } else if (strcmp(argv[i], "--max-iters") == 0) {
            if (i + 1 < argc) {
                params.max_iters = std::stoul(argv[++i]);
            } else {
                std::cerr << "Error: --max-iters requires an int argument\n";
                return 1;
            }
        } else if (strcmp(argv[i], "--tolerance") == 0) {
            if (i + 1 < argc) {
                params.tolerance = std::stod(argv[++i]);
            } else {
                std::cerr << "Error: --tolerance requires a float argument\n";
                return 1;
            }
        } else if (strcmp(argv[i], "--zero-threshold") == 0) {
            if (i + 1 < argc) {
                params.zero_threshold = std::stod(argv[++i]);
            } else {
                std::cerr << "Error: --zero-threshold requires a float argument\n";
                return 1;
            }
        } else if (strcmp(argv[i], "--init-by-length") == 0) {
            params.init_by_length = true;
        } else if (strcmp(argv[i], "--uniform-init") == 0) {
            params.use_uniform_init = true;
        } else if (strcmp(argv[i], "--threads") == 0) {
            if (i + 1 < argc) {
                params.threads = std::stoi(argv[++i]);
            } else {
                std::cerr << "Error: --threads requires an int argument\n";
                return 1;
            }
        } else if (strcmp(argv[i], "--debug-trace") == 0) {
            if (i + 1 < argc) {
                params.debug_file = argv[++i];
                params.debug_trace = true;
                // Parse transcript IDs from environment variable
                const char* debug_txps = std::getenv("EM_QUANT_DEBUG_TXPS");
                if (debug_txps) {
                    std::string txps_str(debug_txps);
                    size_t pos = 0;
                    while ((pos = txps_str.find(',')) != std::string::npos) {
                        params.debug_transcripts.push_back(txps_str.substr(0, pos));
                        txps_str.erase(0, pos + 1);
                    }
                    if (!txps_str.empty()) {
                        params.debug_transcripts.push_back(txps_str);
                    }
                }
            } else {
                std::cerr << "Error: --debug-trace requires a file argument\n";
                return 1;
            }
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Error: Unknown option: " << argv[i] << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }
    
    // Check required arguments
    if (ec_file.empty()) {
        std::cerr << "Error: -e/--ec is required\n";
        print_usage(argv[0]);
        return 1;
    }
    if (output_file.empty()) {
        std::cerr << "Error: -o/--output is required\n";
        print_usage(argv[0]);
        return 1;
    }
    
    try {
        // Load equivalence classes
        if (verbose) {
            std::cerr << "Loading equivalence classes from: " << ec_file << "\n";
        }
        std::vector<std::string> transcript_names;
        ECTable ecs = load_ec_file(ec_file, transcript_names);
        
        if (verbose) {
            std::cerr << "Loaded " << ecs.n_transcripts << " transcripts and " 
                      << ecs.n_ecs << " equivalence classes\n";
        }
        
        // Initialize transcript state
        TranscriptState state;
        state.resize(ecs.n_transcripts);
        
        // Copy transcript names
        for (size_t i = 0; i < ecs.n_transcripts; ++i) {
            state.names[i] = transcript_names[i];
        }
        
        // Load transcript lengths
        if (!lengths_file.empty()) {
            if (verbose) {
                std::cerr << "Loading transcript lengths from: " << lengths_file << "\n";
            }
            load_transcript_lengths(state, lengths_file);
        } else {
            if (verbose) {
                std::cerr << "No lengths file provided, using placeholder lengths\n";
            }
            load_transcript_lengths(state, "");
        }
        
        // Override effective lengths if GC-corrected file provided
        if (!effective_lengths_file.empty()) {
            if (verbose) {
                std::cerr << "Loading GC-corrected effective lengths from: " << effective_lengths_file << "\n";
            }
            load_effective_lengths(state, effective_lengths_file);
        }
        
        // Run EM or VB
        EMResult result;
        if (params.use_vb) {
            if (verbose) {
                std::cerr << "Running Variational Bayes with prior=" << params.vb_prior << "\n";
            }
            result = run_vb(ecs, state, params);
        } else {
            if (verbose) {
                std::cerr << "Running EM algorithm\n";
            }
            result = run_em(ecs, state, params);
        }
        
        if (verbose) {
            std::cerr << "Converged: " << (result.converged ? "yes" : "no") << "\n";
            std::cerr << "Iterations: " << result.iterations << "\n";
            std::cerr << "Final " << (params.use_vb ? "ELBO" : "log-likelihood") 
                      << ": " << result.final_ll << "\n";
        }
        
        // Write output TSV
        std::ofstream out(output_file);
        if (!out.is_open()) {
            std::cerr << "Error: Failed to open output file: " << output_file << "\n";
            return 1;
        }
        
        // Write header
        out << "Name\tLength\tEffectiveLength\tTPM\tNumReads\n";
        
        // Write transcript data
        for (size_t i = 0; i < state.n; ++i) {
            out << state.names[i] << "\t"
                << state.lengths[i] << "\t"
                << state.eff_lengths[i] << "\t"
                << result.tpm[i] << "\t"
                << result.counts[i] << "\n";
        }
        
        out.close();
        
        if (verbose) {
            std::cerr << "Output written to: " << output_file << "\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

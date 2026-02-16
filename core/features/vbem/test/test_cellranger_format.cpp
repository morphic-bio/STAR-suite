/**
 * Test driver for CellRangerFormatter
 * 
 * Usage: ./test_cellranger_format <input_fasta> <input_gtf> <output_fasta> <output_gtf> [legacy|updated]
 * 
 * Example:
 *   ./test_cellranger_format \
 *     test_genome.fa \
 *     test_genes.gtf \
 *     output_genome.fa \
 *     output_genes.gtf
 */

#include <iostream>
#include <cstring>
#include "../source/CellRangerFormatter.h"

void printUsage(const char* progName) {
    std::cerr << "Usage: " << progName << " <input_fasta> <input_gtf> <output_fasta> <output_gtf> [legacy|updated]\n";
    std::cerr << "\nExample:\n";
    std::cerr << "  " << progName << " test_genome.fa test_genes.gtf output_genome.fa output_genes.gtf\n";
    std::cerr << "  " << progName << " test_genome.fa test_genes.gtf output_genome.fa output_genes.gtf updated\n";
}

int main(int argc, char** argv) {
    if (argc < 5) {
        printUsage(argv[0]);
        return 1;
    }
    
    CellRangerFormatter::Config config;
    config.inputFastaPath = argv[1];
    config.inputGtfPath = argv[2];
    config.outputFastaPath = argv[3];
    config.outputGtfPath = argv[4];
    config.useLegacyGtfFilter = true;

    if (argc >= 6) {
        std::string mode = argv[5];
        if (mode == "legacy") {
            config.useLegacyGtfFilter = true;
        } else if (mode == "updated") {
            config.useLegacyGtfFilter = false;
        } else {
            std::cerr << "ERROR: invalid mode '" << mode << "'. Use legacy or updated.\n";
            printUsage(argv[0]);
            return 1;
        }
    }
    
    std::cout << "CellRangerFormatter Test Driver\n";
    std::cout << "===============================\n";
    std::cout << "Input FASTA: " << config.inputFastaPath << "\n";
    std::cout << "Input GTF: " << config.inputGtfPath << "\n";
    std::cout << "Output FASTA: " << config.outputFastaPath << "\n";
    std::cout << "Output GTF: " << config.outputGtfPath << "\n";
    std::cout << "GTF Filter Mode: " << (config.useLegacyGtfFilter ? "legacy" : "updated") << "\n";
    std::cout << "\n";
    
    CellRangerFormatter::Result result = CellRangerFormatter::format(config);
    
    if (!result.success) {
        std::cerr << "ERROR: " << result.errorMessage << "\n";
        return 1;
    }
    
    std::cout << "SUCCESS: Files formatted successfully\n";
    return 0;
}

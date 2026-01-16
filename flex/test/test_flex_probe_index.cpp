/**
 * Test driver for FlexProbeIndex
 * 
 * Usage: ./test_flex_probe_index <probe_csv> <gtf> <base_fasta> <output_dir>
 * 
 * Example:
 *   ./test_flex_probe_index \
 *     /mnt/pikachu/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv \
 *     /home/lhhung/cellranger-9.0.1/external/cellranger_tiny_ref/genes/genes.gtf.gz \
 *     /home/lhhung/cellranger-9.0.1/external/cellranger_tiny_ref/fasta/genome.fa \
 *     ./test_flex_probe_output
 */

#include <iostream>
#include <cstring>
#include "../source/FlexProbeIndex.h"

void printUsage(const char* progName) {
    std::cerr << "Usage: " << progName << " <probe_csv> <gtf> <base_fasta> <output_dir> [--enforce-length N]\n";
    std::cerr << "\nExample:\n";
    std::cerr << "  " << progName << " probes.csv genes.gtf.gz genome.fa ./output\n";
}

int main(int argc, char** argv) {
    if (argc < 5) {
        printUsage(argv[0]);
        return 1;
    }
    
    FlexProbeIndex::Config config;
    config.probeCSVPath = argv[1];
    config.gtfPath = argv[2];
    config.baseFastaPath = argv[3];
    config.outputDir = argv[4];
    config.enforceProbeLength = 50;
    
    // Parse optional arguments
    for (int i = 5; i < argc; i++) {
        if (std::strcmp(argv[i], "--enforce-length") == 0 && i + 1 < argc) {
            config.enforceProbeLength = std::atoi(argv[++i]);
        }
    }
    
    std::cout << "FlexProbeIndex Test Driver\n";
    std::cout << "==========================\n";
    std::cout << "Probe CSV: " << config.probeCSVPath << "\n";
    std::cout << "GTF: " << config.gtfPath << "\n";
    std::cout << "Base FASTA: " << config.baseFastaPath << "\n";
    std::cout << "Output dir: " << config.outputDir << "\n";
    std::cout << "Enforce length: " << config.enforceProbeLength << "bp\n";
    std::cout << "\n";
    
    FlexProbeIndex::Result result = FlexProbeIndex::run(config);
    
    if (result.success) {
        std::cout << "SUCCESS!\n\n";
        std::cout << "Statistics:\n";
        std::cout << "  Total input probes:    " << result.stats.totalInput << "\n";
        std::cout << "  Dropped (DEPRECATED):  " << result.stats.droppedDeprecated << "\n";
        std::cout << "  Dropped (no GTF match):" << result.stats.droppedNoMatch << "\n";
        std::cout << "  Dropped (invalid seq): " << result.stats.droppedInvalidSeq << "\n";
        std::cout << "  Total output probes:   " << result.stats.totalOutput << "\n";
        std::cout << "  Unique genes:          " << result.stats.uniqueGenes << "\n";
        return 0;
    } else {
        std::cerr << "ERROR: " << result.errorMessage << "\n";
        return 1;
    }
}


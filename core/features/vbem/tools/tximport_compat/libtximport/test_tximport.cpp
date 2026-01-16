#include "tximport.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>

using namespace tximport;

// Helper for floating point comparison
bool approx_equal(double a, double b, double tol = 1e-6) {
    if (std::abs(a) < tol && std::abs(b) < tol) return true;
    return std::abs(a - b) / std::max(std::abs(a), std::abs(b)) < tol;
}

void test_basic_summarization() {
    std::cout << "Test: basic_summarization... ";

    // Create test data: 4 transcripts mapping to 2 genes
    std::vector<TxRecord> transcripts = {
        {"tx1", 1000, 800, 100.0, 80.0},   // gene1
        {"tx2", 500,  400, 50.0,  20.0},   // gene1
        {"tx3", 2000, 1800, 200.0, 360.0}, // gene2
        {"tx4", 1500, 1300, 100.0, 130.0}, // gene2
    };

    std::unordered_map<std::string, std::string> tx2gene = {
        {"tx1", "gene1"},
        {"tx2", "gene1"},
        {"tx3", "gene2"},
        {"tx4", "gene2"},
    };

    SummarizeStats stats;
    auto result = summarize_to_gene(transcripts, tx2gene, 
                                    CountsFromAbundance::LengthScaledTPM,
                                    LengthMode::Effective, nullptr, &stats);

    assert(result.size() == 2);
    assert(stats.tx_total == 4);
    assert(stats.tx_mapped == 4);
    assert(stats.tx_missing == 0);
    assert(stats.genes_output == 2);

    // Gene1: TPM = 100 + 50 = 150
    // Gene1: EffLen = (100*800 + 50*400) / 150 = 100000/150 = 666.67
    assert(result[0].name == "gene1");
    assert(approx_equal(result[0].tpm, 150.0));
    assert(approx_equal(result[0].eff_len, 666.666667, 0.01));

    // Gene2: TPM = 200 + 100 = 300
    // Gene2: EffLen = (200*1800 + 100*1300) / 300 = 490000/300 = 1633.33
    assert(result[1].name == "gene2");
    assert(approx_equal(result[1].tpm, 300.0));
    assert(approx_equal(result[1].eff_len, 1633.333333, 0.01));

    // lengthScaledTPM counts:
    // Total raw counts = 80 + 20 + 360 + 130 = 590
    // gene1_scaled = 150 * 666.67 = 100000
    // gene2_scaled = 300 * 1633.33 = 490000
    // total_scaled = 590000
    // scale = 590 / 590000 = 0.001
    // gene1_final = 100000 * 0.001 = 100
    // gene2_final = 490000 * 0.001 = 490
    assert(approx_equal(result[0].num_reads, 100.0, 0.01));
    assert(approx_equal(result[1].num_reads, 490.0, 0.01));

    // Verify total counts preserved
    double total = result[0].num_reads + result[1].num_reads;
    assert(approx_equal(total, 590.0, 0.01));

    std::cout << "PASS" << std::endl;
}

void test_missing_transcripts() {
    std::cout << "Test: missing_transcripts... ";

    std::vector<TxRecord> transcripts = {
        {"tx1", 1000, 800, 100.0, 80.0},
        {"tx2", 500,  400, 50.0,  20.0},
        {"tx_unknown", 2000, 1800, 200.0, 360.0},  // not in tx2gene
    };

    std::unordered_map<std::string, std::string> tx2gene = {
        {"tx1", "gene1"},
        {"tx2", "gene1"},
    };

    SummarizeStats stats;
    auto result = summarize_to_gene(transcripts, tx2gene,
                                    CountsFromAbundance::LengthScaledTPM,
                                    LengthMode::Effective, nullptr, &stats);

    assert(result.size() == 1);
    assert(stats.tx_missing == 1);
    assert(stats.tx_mapped == 2);
    assert(result[0].name == "gene1");

    std::cout << "PASS" << std::endl;
}

void test_zero_tpm_gene() {
    std::cout << "Test: zero_tpm_gene... ";

    std::vector<TxRecord> transcripts = {
        {"tx1", 1000, 800, 100.0, 80.0},   // gene1
        {"tx2", 500,  400, 0.0,   0.0},    // gene2 - zero TPM
    };

    std::unordered_map<std::string, std::string> tx2gene = {
        {"tx1", "gene1"},
        {"tx2", "gene2"},
    };

    SummarizeStats stats;
    auto result = summarize_to_gene(transcripts, tx2gene,
                                    CountsFromAbundance::LengthScaledTPM,
                                    LengthMode::Effective, nullptr, &stats);

    assert(result.size() == 2);
    assert(stats.genes_zero_tpm == 1);

    // Gene2 should use fallback length (selected length average = EffectiveLength)
    assert(result[1].name == "gene2");
    assert(approx_equal(result[1].tpm, 0.0));
    assert(approx_equal(result[1].eff_len, 400.0));  // fallback to effective length (400)
    assert(approx_equal(result[1].num_reads, 0.0));

    std::cout << "PASS" << std::endl;
}

void test_gene_order() {
    std::cout << "Test: gene_order (alphabetical, case-insensitive)... ";

    std::vector<TxRecord> transcripts = {
        {"tx3", 1000, 800, 100.0, 80.0},   // GeneB
        {"tx1", 500,  400, 50.0,  20.0},   // geneA (lowercase)
        {"tx2", 2000, 1800, 200.0, 360.0}, // geneA again
    };

    std::unordered_map<std::string, std::string> tx2gene = {
        {"tx1", "geneA"},   // lowercase 'g'
        {"tx2", "geneA"},
        {"tx3", "GeneB"},   // uppercase 'G'
    };

    std::vector<std::string> gene_order;
    auto result = summarize_to_gene(transcripts, tx2gene,
                                    CountsFromAbundance::LengthScaledTPM,
                                    LengthMode::Effective, &gene_order, nullptr);

    // Gene order should be alphabetical (case-insensitive)
    // geneA < GeneB (case-insensitive: 'a' < 'b')
    assert(gene_order.size() == 2);
    assert(gene_order[0] == "geneA");  // 'a' before 'b'
    assert(gene_order[1] == "GeneB");

    std::cout << "PASS" << std::endl;
}

void test_parse_quant_sf() {
    std::cout << "Test: parse_quant_sf... ";

    // Create temp file
    const char* test_file = "/tmp/test_quant.sf";
    {
        std::ofstream f(test_file);
        f << "Name\tLength\tEffectiveLength\tTPM\tNumReads\n";
        f << "tx1\t1000\t800\t100.5\t80.25\n";
        f << "tx2\t500\t400.5\t50.123456\t20.0\n";
    }

    std::vector<TxRecord> records;
    std::string error;
    bool ok = parse_quant_sf(test_file, records, error);

    assert(ok);
    assert(records.size() == 2);
    assert(records[0].name == "tx1");
    assert(approx_equal(records[0].length, 1000.0));
    assert(approx_equal(records[0].eff_len, 800.0));
    assert(approx_equal(records[0].tpm, 100.5));
    assert(approx_equal(records[0].num_reads, 80.25));

    assert(records[1].name == "tx2");
    assert(approx_equal(records[1].eff_len, 400.5));
    assert(approx_equal(records[1].tpm, 50.123456));

    std::remove(test_file);
    std::cout << "PASS" << std::endl;
}

void test_parse_tx2gene() {
    std::cout << "Test: parse_tx2gene... ";

    const char* test_file = "/tmp/test_tx2gene.tsv";
    {
        std::ofstream f(test_file);
        f << "tx1\tgene1\n";
        f << "tx2\tgene1\n";
        f << "tx3\tgene2\n";
    }

    std::unordered_map<std::string, std::string> tx2gene;
    std::vector<std::string> gene_order;
    std::string error;
    bool ok = parse_tx2gene(test_file, tx2gene, gene_order, error);

    assert(ok);
    assert(tx2gene.size() == 3);
    assert(tx2gene["tx1"] == "gene1");
    assert(tx2gene["tx3"] == "gene2");
    assert(gene_order.size() == 2);
    assert(gene_order[0] == "gene1");  // first occurrence
    assert(gene_order[1] == "gene2");

    std::remove(test_file);
    std::cout << "PASS" << std::endl;
}

void test_write_gene_sf() {
    std::cout << "Test: write_gene_sf... ";

    std::vector<GeneSummary> summaries = {
        {"gene1", 750.0, 600.0, 150.0, 100.0},
        {"gene2", 1750.0, 1550.0, 300.0, 490.0},
    };

    const char* test_file = "/tmp/test_genes.sf";
    std::string error;
    bool ok = write_gene_sf(test_file, summaries, error, 6);
    assert(ok);

    // Read back and verify
    std::ifstream f(test_file);
    std::string line;
    std::getline(f, line);  // header
    assert(line == "Name\tLength\tEffectiveLength\tTPM\tNumReads");

    std::getline(f, line);
    assert(line.substr(0, 5) == "gene1");

    std::remove(test_file);
    std::cout << "PASS" << std::endl;
}

void test_scaled_tpm_mode() {
    std::cout << "Test: scaledTPM mode... ";

    std::vector<TxRecord> transcripts = {
        {"tx1", 1000, 800, 100.0, 80.0},
        {"tx2", 500,  400, 100.0, 20.0},  // same TPM, different counts
    };

    std::unordered_map<std::string, std::string> tx2gene = {
        {"tx1", "gene1"},
        {"tx2", "gene2"},
    };

    auto result = summarize_to_gene(transcripts, tx2gene,
                                    CountsFromAbundance::ScaledTPM,
                                    LengthMode::Effective, nullptr, nullptr);

    // scaledTPM: counts proportional to TPM
    // Both genes have TPM=100, so they get equal counts
    // Total raw = 100, so each gets 50
    assert(approx_equal(result[0].num_reads, 50.0));
    assert(approx_equal(result[1].num_reads, 50.0));

    std::cout << "PASS" << std::endl;
}

void test_no_counts_transformation() {
    std::cout << "Test: no counts transformation... ";

    std::vector<TxRecord> transcripts = {
        {"tx1", 1000, 800, 100.0, 80.0},
        {"tx2", 500,  400, 100.0, 20.0},
    };

    std::unordered_map<std::string, std::string> tx2gene = {
        {"tx1", "gene1"},
        {"tx2", "gene2"},
    };

    auto result = summarize_to_gene(transcripts, tx2gene,
                                    CountsFromAbundance::No,
                                    LengthMode::Effective, nullptr, nullptr);

    // No transformation: raw counts preserved
    assert(approx_equal(result[0].num_reads, 80.0));
    assert(approx_equal(result[1].num_reads, 20.0));

    std::cout << "PASS" << std::endl;
}

int main() {
    std::cout << "\n=== libtximport Unit Tests ===" << std::endl;

    test_basic_summarization();
    test_missing_transcripts();
    test_zero_tpm_gene();
    test_gene_order();
    test_parse_quant_sf();
    test_parse_tx2gene();
    test_write_gene_sf();
    test_scaled_tpm_mode();
    test_no_counts_transformation();

    std::cout << "\nAll tests PASSED" << std::endl;
    return 0;
}


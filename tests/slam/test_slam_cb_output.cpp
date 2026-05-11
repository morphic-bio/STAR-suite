#include "SlamQuant.h"
#include "libem/slam_vb_overdisp.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <unistd.h>

static bool approx(double a, double b, double tol = 1e-9) {
    return std::fabs(a - b) <= tol;
}

static int check(bool ok, const std::string& label) {
    if (!ok) {
        std::cerr << "FAIL: " << label << "\n";
        return 1;
    }
    return 0;
}

static std::string readFile(const std::string& path) {
    std::ifstream in(path.c_str());
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

static std::vector<std::string> splitTab(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    std::istringstream ss(line);
    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }
    return fields;
}

struct Sums {
    double readCount = 0.0;
    double conversions = 0.0;
    double coverage = 0.0;
};

static std::unordered_map<std::string, Sums> reconstructStarCb(const std::string& path) {
    std::unordered_map<std::string, Sums> sums;
    std::ifstream in(path.c_str());
    std::string line;
    bool header = true;
    while (std::getline(in, line)) {
        if (header) {
            header = false;
            continue;
        }
        std::vector<std::string> f = splitTab(line);
        if (f.size() != 6) {
            continue;
        }
        const std::string& gene = f[1];
        const double nT = std::stod(f[3]);
        const double tc = std::stod(f[4]);
        const double n = std::stod(f[5]);
        sums[gene].readCount += n;
        sums[gene].conversions += tc * n;
        sums[gene].coverage += nT * n;
    }
    return sums;
}

int main() {
    int failed = 0;

    const MismatchHistogramKey wideKey = slamPackMismatchKey(300, 260);
    failed += check(slamKeyNT(wideKey) == 300, "wide key nT round-trip");
    failed += check(slamKeyTC(wideKey) == 260, "wide key TC round-trip");
    failed += check(wideKey > 0xFFFFu, "wide key exceeds old uint16_t key space");

    SlamQuant quant(2);
    quant.addRead(0, 300, 260, 2.5);
    quant.addRead(0, 5, 0, 1.0);
    quant.addRead(1, 7, 1, 3.0);

    const auto& genes = quant.genes();
    failed += check(genes[0].histogram.count(wideKey) == 1, "SlamQuant stores TC>255 bin");
    failed += check(approx(genes[0].histogram.at(wideKey), 2.5), "SlamQuant stores weighted bin count");
    failed += check(approx(genes[0].readCount, 3.5), "gene0 read count");
    failed += check(approx(genes[0].conversions, 650.0), "gene0 conversions");
    failed += check(approx(genes[0].coverage, 755.0), "gene0 coverage");

    SlamSolver solver(0.001, 0.90);
    MismatchHistogram hist;
    hist[wideKey] = 2.0;
    hist[slamPackMismatchKey(300, 0)] = 2.0;
    SlamResult solved = solver.solve(hist);
    failed += check(std::isfinite(solved.ntr) && std::isfinite(solved.log_likelihood),
                    "SlamSolver accepts TC>255 histogram bin");

    SlamVbOverdispSolver vb(0.001, 0.90, 50.0, 1.0, 1.0);
    VbOverdispResult vbSolved = vb.solve(hist);
    failed += check(std::isfinite(vbSolved.ntr_map) && std::isfinite(vbSolved.log_likelihood),
                    "SlamVbOverdispSolver accepts TC>255 histogram bin");

    const char* tmpEnv = std::getenv("TMP_DIR");
    const std::string tmpDir = (tmpEnv && tmpEnv[0] != '\0') ? tmpEnv : "/tmp";
    const std::string base = tmpDir + "/test_slam_cb_output_" + std::to_string(static_cast<long long>(getpid()));
    const std::vector<std::string> geneIds{"ENSG1", "ENSG2"};
    const std::vector<std::string> geneNames{"GeneA", ""};

    const std::string starPath = base + ".star.tsv";
    failed += check(quant.writeCountBinomial(geneIds, geneNames, starPath, "S1_", "star"),
                    "write STAR cB format");
    const std::string star = readFile(starPath);
    failed += check(star.find("sample\tfeature_id\tfeature_name\tnT\tTC\tn\n") == 0,
                    "STAR cB header");
    failed += check(star.find("S1\tENSG1\tGeneA\t300\t260\t2.5\n") != std::string::npos,
                    "STAR cB includes TC>255 row");
    failed += check(star.find("S1\tENSG2\tENSG2\t7\t1\t3\n") != std::string::npos,
                    "STAR cB falls back to gene ID for empty name");

    const auto sums = reconstructStarCb(starPath);
    failed += check(approx(sums.at("ENSG1").readCount, genes[0].readCount), "cB reconstructs read count");
    failed += check(approx(sums.at("ENSG1").conversions, genes[0].conversions), "cB reconstructs conversions");
    failed += check(approx(sums.at("ENSG1").coverage, genes[0].coverage), "cB reconstructs coverage");

    const std::string ezPath = base + ".ezbakr.tsv";
    failed += check(quant.writeCountBinomial(geneIds, geneNames, ezPath, "S1_", "ezbakr"),
                    "write EZbakR-like cB format");
    const std::string ez = readFile(ezPath);
    failed += check(ez.find("sample\trname\tsj\tGF\tXF\tnT\tTC\tn\n") == 0,
                    "EZbakR-like cB header");
    failed += check(ez.find("S1\t.\t.\tENSG1\t.\t300\t260\t2.5\n") != std::string::npos,
                    "EZbakR-like cB includes TC>255 row");

    const std::string gzPath = base + ".star.tsv.gz";
    failed += check(quant.writeCountBinomial(geneIds, geneNames, gzPath, "S1_", "star"),
                    "write gzip cB format");
    std::ifstream gz(gzPath.c_str(), std::ios::binary);
    unsigned char magic[2] = {0, 0};
    gz.read(reinterpret_cast<char*>(magic), 2);
    failed += check(magic[0] == 0x1f && magic[1] == 0x8b, "gzip cB magic");

    std::remove(starPath.c_str());
    std::remove(ezPath.c_str());
    std::remove(gzPath.c_str());

    if (failed == 0) {
        std::cout << "PASS: SLAM cB output checks\n";
    }
    return failed == 0 ? 0 : 1;
}

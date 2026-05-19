#include "VelocytoMexWriter.h"
#include "PfMultiMerge.h"
#include "SoloFeature.h"
#include "TimeFunctions.h"

#include <algorithm>
#include <cctype>
#include <climits>
#include <cerrno>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unordered_map>
#include <utility>
#include <zlib.h>

namespace VelocytoMexWriter {

namespace {

static const int kGzLevel = 3;

static bool fileExists(const string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}

static string normalizeDir(string path) {
    while (!path.empty() && path.back() == '/') {
        path.pop_back();
    }
    return path.empty() ? "." : path;
}

static bool dirExists(const string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

static bool ensureDir(const string& path) {
    const string dir = normalizeDir(path);
    if (dir.empty() || dir == ".") {
        return true;
    }
    if (dirExists(dir)) {
        return true;
    }

    const size_t slash = dir.find_last_of('/');
    if (slash != string::npos && slash > 0) {
        if (!ensureDir(dir.substr(0, slash))) {
            return false;
        }
    }
    if (mkdir(dir.c_str(), 0775) == 0) {
        return true;
    }
    return errno == EEXIST && dirExists(dir);
}

static string resolveOptionalMexFile(const string& mexDir, const string& basename) {
    const string plain = mexDir + "/" + basename;
    const string gz = plain + ".gz";
    if (fileExists(gz)) return gz;
    if (fileExists(plain)) return plain;
    return string();
}

static bool gzGetLine(gzFile file, string& lineOut) {
    lineOut.clear();
    if (file == nullptr) {
        return false;
    }

    char buffer[8192];
    while (true) {
        char* ret = gzgets(file, buffer, static_cast<int>(sizeof(buffer)));
        if (ret == nullptr) {
            return !lineOut.empty();
        }

        lineOut.append(buffer);
        if (!lineOut.empty() && lineOut.back() == '\n') {
            while (!lineOut.empty() &&
                   (lineOut.back() == '\n' || lineOut.back() == '\r')) {
                lineOut.pop_back();
            }
            return true;
        }
    }
}

static string stripBarcodeSuffix(const string& barcode) {
    const size_t dashPos = barcode.find_last_of('-');
    if (dashPos == string::npos || dashPos == barcode.size() - 1) {
        return barcode;
    }
    for (size_t i = dashPos + 1; i < barcode.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(barcode[i]))) {
            return barcode;
        }
    }
    return barcode.substr(0, dashPos);
}

static void stripBarcodeSuffixesInPlace(vector<string>& barcodes) {
    for (auto& barcode : barcodes) {
        barcode = stripBarcodeSuffix(barcode);
    }
}

static vector<MexWriter::Triplet> readMatrixTriplets(const string& matrixPath,
                                                     size_t nFeatures,
                                                     size_t nBarcodes) {
    vector<MexWriter::Triplet> triplets;
    uint32_t declaredRows = 0;
    uint32_t declaredCols = 0;
    uint64_t declaredNnz = 0;
    bool dimsSeen = false;

    auto parseLine = [&](const string& line) {
        if (line.empty() || line[0] == '%') {
            return;
        }
        if (!dimsSeen) {
            istringstream ss(line);
            if (!(ss >> declaredRows >> declaredCols >> declaredNnz)) {
                throw std::runtime_error("Invalid Matrix Market dimensions in " + matrixPath);
            }
            if (declaredRows != nFeatures || declaredCols != nBarcodes) {
                ostringstream err;
                err << "Velocyto matrix dimensions mismatch for " << matrixPath
                    << ": matrix=" << declaredRows << "x" << declaredCols
                    << " axis=" << nFeatures << "x" << nBarcodes;
                throw std::runtime_error(err.str());
            }
            dimsSeen = true;
            return;
        }

        istringstream ss(line);
        uint32_t row = 0;
        uint32_t col = 0;
        double val = 0.0;
        if (!(ss >> row >> col >> val)) {
            return;
        }
        if (row > 0 && col > 0 && row <= nFeatures && col <= nBarcodes) {
            MexWriter::Triplet t;
            t.gene_idx = row - 1;
            t.cell_idx = col - 1;
            t.count = static_cast<uint32_t>(val);
            triplets.push_back(t);
        }
    };

    const bool isGz = matrixPath.size() > 3 && matrixPath.substr(matrixPath.size() - 3) == ".gz";
    if (isGz) {
        gzFile file = gzopen(matrixPath.c_str(), "rb");
        if (file == nullptr) {
            throw std::runtime_error("Failed to open gz matrix: " + matrixPath);
        }
        string line;
        while (gzGetLine(file, line)) {
            parseLine(line);
        }
        const int rc = gzclose(file);
        if (rc != Z_OK) {
            throw std::runtime_error("Failed while reading gz matrix: " + matrixPath);
        }
    } else {
        ifstream file(matrixPath.c_str());
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open matrix: " + matrixPath);
        }
        string line;
        while (getline(file, line)) {
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
                line.pop_back();
            }
            parseLine(line);
        }
    }

    if (!dimsSeen) {
        throw std::runtime_error("Missing Matrix Market dimensions in " + matrixPath);
    }
    if (declaredNnz != triplets.size()) {
        // STAR outputs should be exact; make this a hard failure so corruption is not silent.
        ostringstream err;
        err << "Velocyto matrix nnz mismatch for " << matrixPath
            << ": declared=" << declaredNnz << " parsed=" << triplets.size();
        throw std::runtime_error(err.str());
    }
    return triplets;
}

static vector<MexWriter::Triplet> subsetVelocytoTriplets(const vector<MexWriter::Triplet>& triplets,
                                                       const vector<uint32_t>& oldToNew) {
    vector<MexWriter::Triplet> out;
    out.reserve(triplets.size());
    for (const auto& t : triplets) {
        if (t.cell_idx < oldToNew.size() && oldToNew[t.cell_idx] != UINT32_MAX) {
            MexWriter::Triplet nt;
            nt.cell_idx = oldToNew[t.cell_idx];
            nt.gene_idx = t.gene_idx;
            nt.count = t.count;
            out.push_back(nt);
        }
    }
    return out;
}

static bool writeGzLines(const string& path, const vector<string>& lines) {
    gzFile gz = gzopen(path.c_str(), "wb");
    if (gz == nullptr) {
        return false;
    }
    gzbuffer(gz, 1 << 20);
    gzsetparams(gz, kGzLevel, Z_DEFAULT_STRATEGY);
    for (const auto& line : lines) {
        if (gzwrite(gz, line.data(), line.size()) <= 0) {
            gzclose(gz);
            return false;
        }
        if (gzwrite(gz, "\n", 1) <= 0) {
            gzclose(gz);
            return false;
        }
    }
    return gzclose(gz) == Z_OK;
}

static bool writeVelocytoMatrixGz(const string& path,
                                  size_t nRows,
                                  size_t nCols,
                                  const vector<MexWriter::Triplet>& triplets) {
    gzFile gz = gzopen(path.c_str(), "wb");
    if (gz == nullptr) {
        return false;
    }
    gzbuffer(gz, 1 << 20);
    gzsetparams(gz, kGzLevel, Z_DEFAULT_STRATEGY);
    ostringstream hdr;
    hdr << "%%MatrixMarket matrix coordinate integer general\n"
        << "%\n"
        << nRows << " " << nCols << " " << triplets.size() << "\n";
    const string hdrStr = hdr.str();
    if (gzwrite(gz, hdrStr.data(), hdrStr.size()) <= 0) {
        gzclose(gz);
        return false;
    }
    for (const auto& t : triplets) {
        char buf[64];
        const int n = snprintf(buf, sizeof(buf), "%u %u %u\n",
                               t.gene_idx + 1, t.cell_idx + 1, t.count);
        if (n <= 0 || gzwrite(gz, buf, n) <= 0) {
            gzclose(gz);
            return false;
        }
    }
    return gzclose(gz) == Z_OK;
}

static vector<MexWriter::Triplet> sumVelocytoLayers(const vector<MexWriter::Triplet>& a,
                                                    const vector<MexWriter::Triplet>& b,
                                                    const vector<MexWriter::Triplet>& c) {
    map<pair<uint32_t, uint32_t>, uint32_t> sums;
    auto addLayer = [&](const vector<MexWriter::Triplet>& layer) {
        for (const auto& t : layer) {
            sums[{t.cell_idx, t.gene_idx}] += t.count;
        }
    };
    addLayer(a);
    addLayer(b);
    addLayer(c);

    vector<MexWriter::Triplet> out;
    out.reserve(sums.size());
    for (const auto& kv : sums) {
        MexWriter::Triplet t;
        t.cell_idx = kv.first.first;
        t.gene_idx = kv.first.second;
        t.count = kv.second;
        out.push_back(t);
    }
    std::sort(out.begin(), out.end(), [](const MexWriter::Triplet& x, const MexWriter::Triplet& y) {
        if (x.cell_idx != y.cell_idx) {
            return x.cell_idx < y.cell_idx;
        }
        return x.gene_idx < y.gene_idx;
    });
    return out;
}

struct MatrixSummary {
    size_t features = 0;
    size_t barcodes = 0;
    size_t nnzTotal = 0;
    size_t nnzSpliced = 0;
    size_t nnzUnspliced = 0;
    size_t nnzAmbiguous = 0;
};

static MatrixSummary summarize(const VelocytoMexData& data) {
    MatrixSummary s;
    s.features = data.features.size();
    s.barcodes = data.barcodes.size();
    s.nnzSpliced = data.spliced.size();
    s.nnzUnspliced = data.unspliced.size();
    s.nnzAmbiguous = data.ambiguous.size();
    s.nnzTotal = sumVelocytoLayers(data.spliced, data.unspliced, data.ambiguous).size();
    return s;
}

static string jsonEscape(const string& value) {
    string out;
    out.reserve(value.size() + 8);
    for (char c : value) {
        if (c == '\\' || c == '"') {
            out.push_back('\\');
        }
        out.push_back(c);
    }
    return out;
}

static void writeSummaryObject(ofstream& out, const MatrixSummary& s) {
    out << "{\n"
        << "    \"barcodes\": " << s.barcodes << ",\n"
        << "    \"features\": " << s.features << ",\n"
        << "    \"nnz_ambiguous\": " << s.nnzAmbiguous << ",\n"
        << "    \"nnz_spliced\": " << s.nnzSpliced << ",\n"
        << "    \"nnz_total\": " << s.nnzTotal << ",\n"
        << "    \"nnz_unspliced\": " << s.nnzUnspliced << "\n"
        << "  }";
}

static void writeManifest(const string& path,
                          const string& runDir,
                          const string& outputRoot,
                          const string& rawSource,
                          const string& filteredSource,
                          const string& filteredBarcodesSource,
                          const MatrixSummary& raw,
                          const MatrixSummary& filtered) {
    ofstream out(path.c_str());
    if (!out.is_open()) {
        throw std::runtime_error("Failed to write Velocyto manifest: " + path);
    }
    out << "{\n";
    out << "  \"filtered\": ";
    writeSummaryObject(out, filtered);
    out << ",\n";
    out << "  \"output_root\": \"" << jsonEscape(outputRoot) << "\",\n";
    out << "  \"raw\": ";
    writeSummaryObject(out, raw);
    out << ",\n";
    out << "  \"run_dir\": \"" << jsonEscape(runDir) << "\",\n";
    out << "  \"source\": {\n"
        << "    \"filtered_barcodes_source\": \"" << jsonEscape(filteredBarcodesSource) << "\",\n"
        << "    \"velocyto_filtered_dir\": \"" << jsonEscape(filteredSource) << "\",\n"
        << "    \"velocyto_raw_dir\": \"" << jsonEscape(rawSource) << "\"\n"
        << "  }\n";
    out << "}\n";
}

static VelocytoMexData subsetByBarcodeList(const VelocytoMexData& rawData,
                                           const vector<string>& requestedBarcodes,
                                           bool stripSuffixes) {
    unordered_map<string, uint32_t> sourceIndex;
    sourceIndex.reserve(rawData.barcodes.size() * 2);
    for (size_t i = 0; i < rawData.barcodes.size(); ++i) {
        const string key = stripSuffixes ? stripBarcodeSuffix(rawData.barcodes[i]) : rawData.barcodes[i];
        auto inserted = sourceIndex.emplace(key, static_cast<uint32_t>(i));
        if (!inserted.second) {
            throw std::runtime_error("Duplicate Velocyto barcode after normalization: " + key);
        }
    }

    vector<uint32_t> colIndices;
    colIndices.reserve(requestedBarcodes.size());
    vector<string> outputBarcodes;
    outputBarcodes.reserve(requestedBarcodes.size());
    for (const auto& barcode : requestedBarcodes) {
        const string key = stripSuffixes ? stripBarcodeSuffix(barcode) : barcode;
        auto it = sourceIndex.find(key);
        if (it == sourceIndex.end()) {
            throw std::runtime_error("Filtered Velocyto barcode missing from raw barcode list: " + key);
        }
        colIndices.push_back(it->second);
        outputBarcodes.push_back(key);
    }

    VelocytoMexData out = subsetVelocytoColumns(rawData, colIndices);
    out.barcodes = std::move(outputBarcodes);
    return out;
}

} // namespace

bool soloVelocytoRawReady(const string& soloOut) {
    const string rawSource = soloOut + "/Velocyto/raw";
    return !resolveOptionalMexFile(rawSource, "barcodes.tsv").empty() &&
           !resolveOptionalMexFile(rawSource, "features.tsv").empty() &&
           !resolveOptionalMexFile(rawSource, "spliced.mtx").empty() &&
           !resolveOptionalMexFile(rawSource, "unspliced.mtx").empty() &&
           !resolveOptionalMexFile(rawSource, "ambiguous.mtx").empty();
}

VelocytoRunMex loadSoloVelocytoMex(const string& soloOut) {
    const string rawSource = soloOut + "/Velocyto/raw";
    const string filteredSource = soloOut + "/Velocyto/filtered";
    if (!soloVelocytoRawReady(soloOut)) {
        throw std::runtime_error("OCM Velocyto loader: missing Solo.out/Velocyto/raw inputs under " + rawSource);
    }

    string filteredBarcodesPath = resolveOptionalMexFile(filteredSource, "barcodes.tsv");
    if (filteredBarcodesPath.empty()) {
        const string geneFullFiltered = soloOut + "/GeneFull/filtered";
        filteredBarcodesPath = resolveOptionalMexFile(geneFullFiltered, "barcodes.tsv");
    }
    if (filteredBarcodesPath.empty()) {
        const string geneFiltered = soloOut + "/Gene/filtered";
        filteredBarcodesPath = resolveOptionalMexFile(geneFiltered, "barcodes.tsv");
    }
    if (filteredBarcodesPath.empty()) {
        throw std::runtime_error("OCM Velocyto loader: no filtered barcode source found under " + soloOut);
    }

    VelocytoRunMex out;
    out.raw = readVelocytoMex(rawSource);
    if (!resolveOptionalMexFile(filteredSource, "features.tsv").empty() &&
        !resolveOptionalMexFile(filteredSource, "spliced.mtx").empty() &&
        !resolveOptionalMexFile(filteredSource, "unspliced.mtx").empty() &&
        !resolveOptionalMexFile(filteredSource, "ambiguous.mtx").empty()) {
        out.filtered = readVelocytoMex(filteredSource);
    } else {
        out.filtered = subsetByBarcodeList(out.raw, PfMultiMerge::readLines(filteredBarcodesPath), true);
    }
    return out;
}

VelocytoMexData readVelocytoMex(const string& mexDir) {
    VelocytoMexData data;
    const vector<string> featureLines =
        PfMultiMerge::readLines(PfMultiMerge::resolveMexFile(mexDir, "features.tsv"));
    for (const auto& line : featureLines) {
        istringstream ss(line);
        string id;
        string name;
        string type;
        if (getline(ss, id, '\t')) {
            data.features.push_back(id);
            if (getline(ss, name, '\t')) {
                data.featureNames.push_back(name);
                if (getline(ss, type, '\t')) {
                    data.featureTypes.push_back(type);
                } else {
                    data.featureTypes.push_back("Gene Expression");
                }
            } else {
                data.featureNames.push_back(id);
                data.featureTypes.push_back("Gene Expression");
            }
        }
    }
    data.barcodes = PfMultiMerge::readLines(PfMultiMerge::resolveMexFile(mexDir, "barcodes.tsv"));
    data.spliced = readMatrixTriplets(PfMultiMerge::resolveMexFile(mexDir, "spliced.mtx"),
                                      data.features.size(), data.barcodes.size());
    data.unspliced = readMatrixTriplets(PfMultiMerge::resolveMexFile(mexDir, "unspliced.mtx"),
                                        data.features.size(), data.barcodes.size());
    data.ambiguous = readMatrixTriplets(PfMultiMerge::resolveMexFile(mexDir, "ambiguous.mtx"),
                                        data.features.size(), data.barcodes.size());
    return data;
}

VelocytoMexData subsetVelocytoColumns(const VelocytoMexData& data,
                                      const vector<uint32_t>& colIndices) {
    VelocytoMexData out;
    out.features = data.features;
    out.featureNames = data.featureNames;
    out.featureTypes = data.featureTypes;

    vector<uint32_t> oldToNew(data.barcodes.size(), UINT32_MAX);
    out.barcodes.reserve(colIndices.size());
    for (uint32_t oldIdx : colIndices) {
        if (oldIdx >= data.barcodes.size()) {
            continue;
        }
        oldToNew[oldIdx] = static_cast<uint32_t>(out.barcodes.size());
        out.barcodes.push_back(data.barcodes[oldIdx]);
    }
    out.spliced = subsetVelocytoTriplets(data.spliced, oldToNew);
    out.unspliced = subsetVelocytoTriplets(data.unspliced, oldToNew);
    out.ambiguous = subsetVelocytoTriplets(data.ambiguous, oldToNew);
    return out;
}

int writeVelocytoGzDir(const string& outputDir, const VelocytoMexData& data) {
    if (!ensureDir(outputDir)) {
        return -1;
    }

    vector<string> featureLines;
    featureLines.reserve(data.features.size());
    for (size_t i = 0; i < data.features.size(); ++i) {
        const string name = (i < data.featureNames.size()) ? data.featureNames[i] : data.features[i];
        const string type = (i < data.featureTypes.size()) ? data.featureTypes[i] : "Gene Expression";
        featureLines.push_back(data.features[i] + "\t" + name + "\t" + type);
    }
    if (!writeGzLines(outputDir + "/features.tsv.gz", featureLines)) {
        return -1;
    }
    if (!writeGzLines(outputDir + "/barcodes.tsv.gz", data.barcodes)) {
        return -1;
    }

    const vector<MexWriter::Triplet> total =
        sumVelocytoLayers(data.spliced, data.unspliced, data.ambiguous);
    if (!writeVelocytoMatrixGz(outputDir + "/matrix.mtx.gz", data.features.size(), data.barcodes.size(), total)) {
        return -1;
    }
    if (!writeVelocytoMatrixGz(outputDir + "/spliced.mtx.gz", data.features.size(), data.barcodes.size(),
                               data.spliced)) {
        return -1;
    }
    if (!writeVelocytoMatrixGz(outputDir + "/unspliced.mtx.gz", data.features.size(), data.barcodes.size(),
                               data.unspliced)) {
        return -1;
    }
    if (!writeVelocytoMatrixGz(outputDir + "/ambiguous.mtx.gz", data.features.size(), data.barcodes.size(),
                               data.ambiguous)) {
        return -1;
    }
    return 0;
}

int materializeRunVelocytoMex(const string& runDirInput,
                              const string& outputRootInput,
                              ofstream& logStream) {
    const string runDir = normalizeDir(runDirInput);
    const string outputRoot = normalizeDir(outputRootInput);
    const string soloOut = runDir + "/Solo.out";
    const string rawSource = soloOut + "/Velocyto/raw";
    const string filteredSource = soloOut + "/Velocyto/filtered";

    if (resolveOptionalMexFile(rawSource, "barcodes.tsv").empty() ||
        resolveOptionalMexFile(rawSource, "features.tsv").empty() ||
        resolveOptionalMexFile(rawSource, "spliced.mtx").empty() ||
        resolveOptionalMexFile(rawSource, "unspliced.mtx").empty() ||
        resolveOptionalMexFile(rawSource, "ambiguous.mtx").empty()) {
        throw std::runtime_error("Velocyto MEX materializer: missing Solo.out/Velocyto/raw inputs under " + rawSource);
    }

    string filteredBarcodesPath = resolveOptionalMexFile(filteredSource, "barcodes.tsv");
    string filteredBarcodesSource = "Solo.out/Velocyto/filtered/barcodes.tsv";
    if (filteredBarcodesPath.empty()) {
        const string geneFullFiltered = soloOut + "/GeneFull/filtered";
        filteredBarcodesPath = resolveOptionalMexFile(geneFullFiltered, "barcodes.tsv");
        filteredBarcodesSource = "Solo.out/GeneFull/filtered/barcodes.tsv";
    }
    if (filteredBarcodesPath.empty()) {
        const string geneFiltered = soloOut + "/Gene/filtered";
        filteredBarcodesPath = resolveOptionalMexFile(geneFiltered, "barcodes.tsv");
        filteredBarcodesSource = "Solo.out/Gene/filtered/barcodes.tsv";
    }
    if (filteredBarcodesPath.empty()) {
        throw std::runtime_error("Velocyto MEX materializer: no filtered barcode source found");
    }

    logStream << timeMonthDayTime() << " ..... started Velocyto MEX materialization\n";
    VelocytoMexData rawData = readVelocytoMex(rawSource);
    VelocytoMexData rawOutput = rawData;
    stripBarcodeSuffixesInPlace(rawOutput.barcodes);

    vector<string> filteredBarcodes = PfMultiMerge::readLines(filteredBarcodesPath);
    VelocytoMexData filteredOutput = subsetByBarcodeList(rawData, filteredBarcodes, true);

    const string rawOut = outputRoot + "/raw_velocyto_feature_bc_matrix";
    const string filteredOut = outputRoot + "/filtered_velocyto_feature_bc_matrix";
    if (writeVelocytoGzDir(rawOut, rawOutput) != 0) {
        throw std::runtime_error("Velocyto MEX materializer: failed writing " + rawOut);
    }
    if (writeVelocytoGzDir(filteredOut, filteredOutput) != 0) {
        throw std::runtime_error("Velocyto MEX materializer: failed writing " + filteredOut);
    }

    writeManifest(outputRoot + "/velocyto_feature_bc_matrix_manifest.json",
                  runDir,
                  outputRoot,
                  rawSource,
                  filteredSource,
                  filteredBarcodesSource,
                  summarize(rawOutput),
                  summarize(filteredOutput));
    logStream << timeMonthDayTime() << " ..... finished Velocyto MEX materialization\n";
    return 0;
}

int runVelocytoMexMaterialize(Parameters& P) {
    if (!P.pSolo.featureYes[SoloFeatureTypes::Velocyto]) {
        return 0;
    }
    if (P.pSolo.useInlineReplayer) {
        P.inOut->logMain << "WARNING: skipping native Velocyto outs/ MEX materialization "
                         << "because inline replayer mode does not expose the standard Solo.out/Velocyto tree\n";
        return 0;
    }

    const string runDir = normalizeDir(P.outFileNamePrefix);
    const string outputRoot = runDir + "/outs";
    try {
        return materializeRunVelocytoMex(runDir, outputRoot, P.inOut->logMain);
    } catch (const std::exception& ex) {
        P.inOut->logMain << "EXITING because of fatal VELOCYTO MEX MATERIALIZATION error: "
                         << ex.what() << "\n";
        return 1;
    }
}

} // namespace VelocytoMexWriter

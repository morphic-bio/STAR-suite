#include "FlexGdna.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <sys/stat.h>

namespace {

bool regularFile(const std::string& path)
{
    struct stat st {};
    return !path.empty() && stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

std::string directoryName(const std::string& path)
{
    const size_t slash = path.find_last_of("/\\");
    return slash == std::string::npos ? std::string(".") : path.substr(0, slash);
}

std::string lower(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

std::string trim(const std::string& value)
{
    const size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return std::string();
    const size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::vector<std::string> parseCsvRow(const std::string& line)
{
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                field.push_back('"');
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (c == ',' && !quoted) {
            fields.push_back(trim(field));
            field.clear();
        } else {
            field.push_back(c);
        }
    }
    fields.push_back(trim(field));
    return fields;
}

bool includedValue(const std::string& value)
{
    const std::string v = lower(trim(value));
    return v.empty() || v == "true" || v == "1" || v == "yes";
}

bool deprecatedProbeId(const std::string& probeId)
{
    return lower(probeId).find("deprecated") != std::string::npos;
}

FlexGdnaRegion parseRegion(const std::string& value)
{
    const std::string v = lower(trim(value));
    if (v == "spliced")
        return FlexGdnaSpliced;
    if (v == "unspliced")
        return FlexGdnaUnspliced;
    return FlexGdnaUnknown;
}

struct Model {
    double constant = 0.0;
    double slope = 0.0;
    double pivot = 0.0;
    double rss = std::numeric_limits<double>::infinity();
};

Model fitAtPivot(const std::vector<double>& x,
                 const std::vector<double>& y,
                 size_t pivotIndex)
{
    const size_t n = x.size();
    const double pivot = x[pivotIndex];
    double sumZ = 0.0;
    double sumZZ = 0.0;
    double sumY = 0.0;
    double sumZY = 0.0;

    for (size_t i = 0; i < n; ++i) {
        const double z = i < pivotIndex ? 0.0 : x[i] - pivot;
        sumZ += z;
        sumZZ += z * z;
        sumY += y[i];
        sumZY += z * y[i];
    }

    const double a = static_cast<double>(n);
    const double determinant = a * sumZZ - sumZ * sumZ;
    Model out;
    out.pivot = pivot;
    if (determinant == 0.0) {
        out.constant = a == 0.0 ? 0.0 : sumY / a;
        out.slope = 0.0;
    } else {
        out.constant = (sumZZ * sumY - sumZ * sumZY) / determinant;
        out.slope = (-sumZ * sumY + a * sumZY) / determinant;
    }

    out.rss = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double z = i < pivotIndex ? 0.0 : x[i] - pivot;
        const double residual = y[i] - out.constant - out.slope * z;
        out.rss += residual * residual;
    }
    return out;
}

} // namespace

FlexGdnaProbeMetadata& FlexGdnaProbeMetadata::instance()
{
    static FlexGdnaProbeMetadata metadata;
    return metadata;
}

void FlexGdnaProbeMetadata::reset()
{
    ready_ = false;
    probeCsvPath_.clear();
    probeRegionById_.clear();
    geneProbeCounts_.clear();
    totalProbes_ = 0;
    totalUnsplicedProbes_ = 0;
    controlGeneCount_ = 0;
}

bool FlexGdnaProbeMetadata::load(const std::string& probeCsvPath,
                                 const std::string& probeListPath,
                                 std::string* errorOut)
{
    reset();

    std::ifstream probeList(probeListPath.c_str());
    if (!probeList) {
        if (errorOut)
            *errorOut = "cannot open probe list: " + probeListPath;
        return false;
    }

    std::unordered_map<std::string, uint32_t> geneIndex;
    std::string line;
    uint32_t index = 0;
    while (std::getline(probeList, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#')
            continue;
        ++index;
        if (index > 0x7FFFu) {
            if (errorOut)
                *errorOut = "probe list exceeds 15-bit gene index";
            return false;
        }
        geneIndex[line] = index;
    }
    if (geneIndex.empty()) {
        if (errorOut)
            *errorOut = "probe list contains no genes";
        return false;
    }
    geneProbeCounts_.assign(static_cast<size_t>(index) + 1u, FlexGdnaGeneProbeCounts());

    std::ifstream csv(probeCsvPath.c_str());
    if (!csv) {
        if (errorOut)
            *errorOut = "cannot open filtered probe CSV: " + probeCsvPath;
        reset();
        return false;
    }

    int geneCol = -1;
    int probeIdCol = -1;
    int includedCol = -1;
    int regionCol = -1;
    bool sawHeader = false;
    while (std::getline(csv, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        const std::vector<std::string> fields = parseCsvRow(line);
        if (!sawHeader) {
            for (size_t i = 0; i < fields.size(); ++i) {
                const std::string name = lower(fields[i]);
                if (name == "gene_id")
                    geneCol = static_cast<int>(i);
                else if (name == "probe_id")
                    probeIdCol = static_cast<int>(i);
                else if (name == "included")
                    includedCol = static_cast<int>(i);
                else if (name == "region")
                    regionCol = static_cast<int>(i);
            }
            if (geneCol < 0 || probeIdCol < 0 || regionCol < 0) {
                if (errorOut)
                    *errorOut = "probe CSV requires gene_id, probe_id, and region columns";
                reset();
                return false;
            }
            sawHeader = true;
            continue;
        }

        const size_t required = static_cast<size_t>(
            std::max(std::max(geneCol, probeIdCol), std::max(includedCol, regionCol)));
        if (fields.size() <= required)
            continue;

        const std::string& geneId = fields[static_cast<size_t>(geneCol)];
        const std::string& probeId = fields[static_cast<size_t>(probeIdCol)];
        if ((includedCol >= 0 && !includedValue(fields[static_cast<size_t>(includedCol)]))
            || deprecatedProbeId(probeId)) {
            continue;
        }
        const FlexGdnaRegion region = parseRegion(fields[static_cast<size_t>(regionCol)]);
        const auto geneIt = geneIndex.find(geneId);
        if (geneIt == geneIndex.end())
            continue;

        const auto inserted = probeRegionById_.emplace(probeId, region);
        if (!inserted.second && inserted.first->second != region) {
            if (errorOut)
                *errorOut = "probe ID has conflicting region annotations: " + probeId;
            reset();
            return false;
        }
        if (!inserted.second)
            continue;

        ++totalProbes_;
        FlexGdnaGeneProbeCounts& counts = geneProbeCounts_[geneIt->second];
        if (region == FlexGdnaSpliced) {
            ++counts.spliced;
        } else if (region == FlexGdnaUnspliced) {
            ++counts.unspliced;
            ++totalUnsplicedProbes_;
        }
    }

    if (!sawHeader || totalProbes_ == 0u) {
        if (errorOut)
            *errorOut = "probe CSV contains no included probes on the active gene axis";
        reset();
        return false;
    }

    for (size_t i = 1; i < geneProbeCounts_.size(); ++i) {
        if (geneProbeCounts_[i].spliced > 0 && geneProbeCounts_[i].unspliced > 0)
            ++controlGeneCount_;
    }

    probeCsvPath_ = probeCsvPath;
    ready_ = true;
    return true;
}

FlexGdnaRegion FlexGdnaProbeMetadata::regionForProbeId(const std::string& probeId) const
{
    const auto it = probeRegionById_.find(probeId);
    return it == probeRegionById_.end() ? FlexGdnaUnknown : it->second;
}

std::string FlexGdnaProbeMetadata::discoverProbeCsv(
    const std::string& configuredPath,
    const std::string& probeListPath,
    const std::string& genomeDir)
{
    if (!configuredPath.empty() && configuredPath != "-" && lower(configuredPath) != "auto")
        return regularFile(configuredPath) ? configuredPath : std::string();

    std::vector<std::string> candidates;
    if (!probeListPath.empty() && probeListPath != "-") {
        const std::string probeDir = directoryName(probeListPath);
        candidates.push_back(probeDir + "/filtered_probe_set.csv");
        candidates.push_back(probeDir + "/flex_probe_artifacts/filtered_probe_set.csv");
    }
    if (!genomeDir.empty() && genomeDir != "-") {
        candidates.push_back(genomeDir + "/flex_probe_artifacts/filtered_probe_set.csv");
        candidates.push_back(genomeDir + "/filtered_probe_set.csv");
    }
    for (const std::string& candidate : candidates) {
        if (regularFile(candidate))
            return candidate;
    }
    return std::string();
}

FlexGdnaEstimate flexGdnaEstimate(
    const FlexGdnaProbeMetadata& metadata,
    const std::vector<FlexGdnaGeneMoleculeCounts>& molecules,
    uint64_t classifiedMolecules,
    uint64_t unknownMolecules,
    uint64_t conflictingMolecules,
    uint64_t unassignedMolecules)
{
    FlexGdnaEstimate out;
    out.totalFilteredMolecules =
        classifiedMolecules + unknownMolecules + conflictingMolecules;
    out.classifiedMolecules = classifiedMolecules;
    out.unknownMolecules = unknownMolecules;
    out.conflictingMolecules = conflictingMolecules;
    out.unassignedMolecules = unassignedMolecules;

    if (!metadata.ready()) {
        out.status = "missing_probe_metadata";
        return out;
    }

    const std::vector<FlexGdnaGeneProbeCounts>& probes = metadata.geneProbeCounts();
    const size_t n = std::min(probes.size(), molecules.size());
    std::vector<std::pair<double, double>> points;
    points.reserve(metadata.controlGeneCount());
    for (size_t i = 1; i < n; ++i) {
        if (probes[i].spliced == 0 || probes[i].unspliced == 0)
            continue;
        const double spliced = static_cast<double>(molecules[i].spliced)
                             / static_cast<double>(probes[i].spliced);
        const double unspliced = static_cast<double>(molecules[i].unspliced)
                               / static_cast<double>(probes[i].unspliced);
        points.push_back(std::make_pair(std::log1p(spliced), std::log1p(unspliced)));
    }
    out.controlGenes = static_cast<uint32_t>(points.size());
    if (points.size() < 10u) {
        out.status = "insufficient_control_genes";
        return out;
    }
    if (out.totalFilteredMolecules == 0u) {
        out.status = "no_gene_assigned_filtered_molecules";
        return out;
    }

    std::stable_sort(points.begin(), points.end(),
                     [](const std::pair<double, double>& lhs,
                        const std::pair<double, double>& rhs) {
                         return lhs.first < rhs.first;
                     });
    std::vector<double> x(points.size());
    std::vector<double> y(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        x[i] = points[i].first;
        y[i] = points[i].second;
    }

    Model best;
    for (size_t pivot = 1; pivot + 1 < x.size(); ++pivot) {
        const Model fit = fitAtPivot(x, y, pivot);
        if (fit.rss < best.rss)
            best = fit;
    }
    if (!std::isfinite(best.rss)) {
        out.status = "model_fit_failed";
        return out;
    }

    out.modelConstant = best.constant;
    out.modelSlope = best.slope;
    out.modelCriticalPoint = best.pivot;
    out.modelRss = best.rss;
    out.estimatedGdnaPerProbe = std::exp(best.constant) - 1.0;
    if (out.estimatedGdnaPerProbe < 0.0)
        out.estimatedGdnaPerProbe = 0.0;
    out.threshold = std::ceil(out.estimatedGdnaPerProbe);

    const double estimatedMolecules = std::round(
        out.estimatedGdnaPerProbe * static_cast<double>(metadata.totalUnsplicedProbes()));
    const double boundedMolecules = std::min(
        static_cast<double>(out.totalFilteredMolecules),
        std::max(0.0, estimatedMolecules));
    out.estimatedGdnaFraction =
        boundedMolecules / static_cast<double>(out.totalFilteredMolecules);
    out.status = "ok";
    out.valid = true;
    return out;
}

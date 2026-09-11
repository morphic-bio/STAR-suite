// Compile with LOADER_SOURCE pointing to the current or preserved FlexFilter.cpp.
#include LOADER_SOURCE
#include <iostream>
int main(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        std::vector<std::string> barcodes, features;
        SampleMatrixData data;
        const bool ok = loadMEXFiles(argv[i], barcodes, features, data);
        std::cout << "CASE " << i << ' ' << ok << '\n';
        if (!ok) continue;
        std::cout << data.nCells << ' ' << data.nGenes << ' ' << data.countMatStride << '\n';
        for (const auto& s : data.barcodes) std::cout << "B " << s << '\n';
        for (const auto& s : data.features) std::cout << "F " << s << '\n';
        for (size_t c = 0; c < data.nCells; ++c)
            std::cout << "C " << c << ' ' << data.nUMIperCB[c] << ' ' << data.nGenePerCB[c] << '\n';
        for (auto i : data.countCellGeneUMIindex) std::cout << "I " << i << '\n';
        for (auto v : data.countCellGeneUMI) std::cout << "V " << v << '\n';
    }
}

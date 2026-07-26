#include "ProbeListIndex.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>

namespace {

void require(bool condition, const std::string &message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        std::exit(1);
    }
}

} // namespace

int main() {
    char pathTemplate[] = "/tmp/star-probe-list-index.XXXXXX";
    const int fd = mkstemp(pathTemplate);
    require(fd >= 0, "could not create temporary probe list");
    close(fd);

    {
        std::ofstream out(pathTemplate);
        require(out.good(), "could not open temporary probe list");
        out << "# ordered Flex feature axis\n"
            << "  ENSG_A  \n"
            << "ENSG_DEPRECATED_probe\n"
            << "ENSG_B\n";
    }

    ProbeListIndex index;
    uint32_t deprecated = 0;
    require(index.load(pathTemplate, true, &deprecated), "filtered load failed");
    require(deprecated == 1, "filtered load did not count deprecated row");
    require(index.size() == 2, "filtered ordered axis has wrong size");
    require(index.orderedGeneIds()[0] == "ENSG_A", "first filtered feature changed order");
    require(index.orderedGeneIds()[1] == "ENSG_B", "second filtered feature changed order");
    require(index.geneIndex15("ENSG_A") == 1, "first filtered feature has wrong index");
    require(index.geneIndex15("ENSG_B") == 2, "second filtered feature has wrong index");

    deprecated = 99;
    require(index.load(pathTemplate, false, &deprecated), "unfiltered load failed");
    require(deprecated == 0, "unfiltered load reported deprecated removal");
    require(index.size() == 3, "unfiltered ordered axis has wrong size");
    require(index.orderedGeneIds()[1] == "ENSG_DEPRECATED_probe",
            "unfiltered ordered axis dropped the deprecated row");
    require(index.geneIndex15("ENSG_B") == 3, "unfiltered index does not match axis order");

    std::remove(pathTemplate);
    std::cout << "PASS: ProbeListIndex preserves the ordered Flex feature axis\n";
    return 0;
}

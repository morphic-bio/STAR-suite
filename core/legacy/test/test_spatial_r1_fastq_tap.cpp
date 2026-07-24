#include "SpatialR1FastqTap.h"

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <string>
#include <thread>

#include <sys/stat.h>
#include <unistd.h>

namespace tap = spatial_r1_fastq_tap;

int main()
{
    char directoryTemplate[] = "/tmp/star-spatial-r1-tap-test-XXXXXX";
    const char *created = ::mkdtemp(directoryTemplate);
    assert(created != nullptr);
    const std::string directory = created;
    const std::string fifo = directory + "/r1.fastq.fifo";
    assert(::mkfifo(fifo.c_str(), 0600) == 0);

    std::string observed;
    std::thread reader([&]() {
        std::ifstream input(fifo.c_str());
        assert(input.good());
        observed.assign(std::istreambuf_iterator<char>(input),
                        std::istreambuf_iterator<char>());
    });

    tap::Writer writer;
    std::string error;
    assert(writer.open(fifo, true, error));
    assert(writer.write("read0", "ACGT", "IIII", error));
    assert(writer.write("read1", "TGCA", "JJJJ", error));
    assert(writer.recordsWritten() == 2);
    assert(!writer.write("bad name", "ACGT", "IIII", error));
    assert(!writer.write("read2", "ACGT", "III", error));
    assert(writer.close(error));
    reader.join();
    assert(observed == "@read0\nACGT\n+\nIIII\n@read1\nTGCA\n+\nJJJJ\n");

    const std::string regular = directory + "/regular.fastq";
    {
        std::ofstream output(regular.c_str());
        assert(output.good());
    }
    tap::Writer rejectsRegular;
    assert(!rejectsRegular.open(regular, true, error));

    assert(::unlink(fifo.c_str()) == 0);
    assert(::unlink(regular.c_str()) == 0);
    assert(::rmdir(directory.c_str()) == 0);
    return 0;
}

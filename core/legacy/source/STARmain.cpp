// The STAR executable: STAR's run with no host hooks (host/StarHost.h).
#include "host/StarHost.h"

int main(int argc, char *argv[]) {
    return star::host::runMain(argc, argv, nullptr);
}

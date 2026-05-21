#include "SoloMemoryProfile.h"
#include "systemFunctions.h"
#include <cstdlib>

bool soloMemoryProfileEnabled()
{
    static int cached = -1;
    if (cached < 0) {
        const char *env = std::getenv("STAR_SOLO_MEMORY_PROFILE");
        cached = (env != nullptr && env[0] != '\0') ? 1 : 0;
    }
    return cached > 0;
}

void soloMemoryProfileCheckpoint(std::ostream &log, const std::string &label)
{
    if (!soloMemoryProfileEnabled()) {
        return;
    }
    log << "Solo memory: " << label << " | " << linuxProcMemory() << std::flush;
}

void soloMemoryProfileCheckpoint(std::ostream &log,
                                 const std::string &label,
                                 const std::string &extra)
{
    if (!soloMemoryProfileEnabled()) {
        return;
    }
    log << "Solo memory: " << label << " | " << extra << " | " << linuxProcMemory() << std::flush;
}

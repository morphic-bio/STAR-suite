#ifndef SOLO_MEMORY_PROFILE_H
#define SOLO_MEMORY_PROFILE_H

#include <ostream>
#include <string>

/** Enabled when STAR_SOLO_MEMORY_PROFILE is set in the environment. */
bool soloMemoryProfileEnabled();

/** Log VmRSS/VmHWM/VmPeak with a phase label (to logMain or stderr). */
void soloMemoryProfileCheckpoint(std::ostream &log, const std::string &label);

/** Same with extra key=value counters (no heap allocation if extra is literal). */
void soloMemoryProfileCheckpoint(std::ostream &log,
                                 const std::string &label,
                                 const std::string &extra);

#endif

#ifndef STAR_HOST_STARHOST_INTERNAL_H
#define STAR_HOST_STARHOST_INTERNAL_H

// STAR-internal helpers behind the host interface. Not for host programs.

#include <fstream>
#include <string>

#include "StarHost.h"

class Parameters;

namespace star {
namespace host {
namespace detail {

// False if runMain was entered before in this process.
bool enterRunMain();
// Label of the External domain for STAR's logs; stable for the process.
const char* externalLabel(const Hooks* hooks);
// Log.out stream used by logMain(); nullptr detaches it.
void attachLog(std::ofstream* logMain);
RunView makeRunView(const Parameters& P);
// The host's message, or the fallback when the host gave none.
std::string hostErrorText(const std::string& error, const char* fallback);
std::string upperCopy(const std::string& text);

}  // namespace detail
}  // namespace host
}  // namespace star

#endif  // STAR_HOST_STARHOST_INTERNAL_H

#include "star_chromap_orchestration.h"

#include "IncludeDefine.h"
#include "Parameters.h"

#include <cctype>
#include <string>

namespace {

std::string lowerCopy(std::string s) {
  for (size_t i = 0; i < s.size(); ++i) {
    s[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(s[i])));
  }
  return s;
}

std::string trimCopy(const std::string &input) {
  size_t first = input.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    return "";
  }
  size_t last = input.find_last_not_of(" \t\r\n");
  return input.substr(first, last - first + 1);
}

bool isYesToken(const std::string &input) {
  const std::string value = lowerCopy(trimCopy(input));
  return value == "yes" || value == "y" || value == "true" || value == "1";
}

}  // namespace

StarChromapAtacAsyncRun::StarChromapAtacAsyncRun() : impl_(nullptr) {}
StarChromapAtacAsyncRun::~StarChromapAtacAsyncRun() {}

bool preflightStarChromapAtacIfEnabled(Parameters &P, bool /*batchModeActive*/) {
  if (P.chromapAtac.enabled == 0 &&
      !isYesToken(P.multiomeAtacPeakMex.inlineMode)) {
    return true;
  }
  P.inOut->logMain
      << "ERROR: --chromapAtacEnable 1 or --multiomeAtacPeakMexInline yes requires a STAR binary built with "
         "Chromap support. Rebuild with WITH_CHROMAP=1 (see docs/LIBCHROMAP_CONTRACT.md).\n";
  return false;
}

bool startStarChromapAtacIfEnabled(Parameters &P,
                                   bool batchModeActive,
                                   StarChromapAtacAsyncRun & /*run*/) {
  return preflightStarChromapAtacIfEnabled(P, batchModeActive);
}

bool runStarChromapAtacIfEnabled(Parameters &P,
                                 bool batchModeActive,
                                 StarChromapAtacAsyncRun & /*run*/) {
  return preflightStarChromapAtacIfEnabled(P, batchModeActive);
}

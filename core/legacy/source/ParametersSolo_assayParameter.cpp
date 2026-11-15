#include "Parameters.h"
#include "ErrorWarning.h"
#include <sstream>

// Sets one registered parameter on behalf of an assay mode and logs the
// effective value. A value given on the command line is kept; when the mode
// owns the parameter, a different command-line value is an error rather than
// a silently ignored request.
void ParametersSolo::assayParameter(const string& assay, const string& mode, const string& logTag,
                                    const string& name, const string& wanted, bool owned) {
    const auto fail = [&](const string& message) {
        exitWithError("EXITING because of " + assay + " parameter error: " + message + "\n",
                      std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    for (auto* p : pP->parArray) {
        if (p->nameString != name) continue;
        if (p->inputLevel > 0) {
            std::ostringstream out; out << *p;
            std::istringstream values(out.str());
            string token, actual;
            while (values >> token) { if (!actual.empty()) actual += ' '; actual += token; }
            if (owned && actual != wanted)
                fail("--" + name + " conflicts with " + mode + " mode (effective value: " + wanted + ")");
        } else {
            std::istringstream in(wanted);
            p->inputValues(in);
        }
        pP->inOut->logMain << "[" << logTag << " effective] " << name << " = " << *p << '\n';
        return;
    }
    fail("internal parameter registration missing: " + name);
}

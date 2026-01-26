/**
 * @file scrna_types.h
 * @brief Common type definitions for libscrna
 * 
 * This header provides standalone type definitions that don't depend on
 * STAR-specific headers (IncludeDefine.h), making libscrna usable by
 * external tools like process_features.
 * 
 * When compiled within STAR context (IncludeDefine.h included first),
 * these types are already defined as macros, so we skip the typedefs.
 */

#ifndef H_scrna_types
#define H_scrna_types

#include <cstdint>
#include <vector>
#include <string>

// Type aliases matching STAR conventions for compatibility
// Skip if already defined as macros (STAR's IncludeDefine.h)
#ifndef uint32
typedef uint32_t uint32;
#endif
#ifndef int64
typedef int64_t int64;
#endif
#ifndef uint64
typedef uint64_t uint64;
#endif

// Use std:: explicitly for clarity
using std::vector;
using std::string;

#endif // H_scrna_types

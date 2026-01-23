/**
 * @file scrna_types.h
 * @brief Common type definitions for libscrna
 * 
 * This header provides standalone type definitions that don't depend on
 * STAR-specific headers (IncludeDefine.h), making libscrna usable by
 * external tools like process_features.
 */

#ifndef H_scrna_types
#define H_scrna_types

#include <cstdint>
#include <vector>
#include <string>

// Type aliases matching STAR conventions for compatibility
typedef uint32_t uint32;
typedef int64_t int64;
typedef uint64_t uint64;

// Use std:: explicitly for clarity
using std::vector;
using std::string;

#endif // H_scrna_types

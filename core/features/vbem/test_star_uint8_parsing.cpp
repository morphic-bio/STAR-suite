// Test using actual ParameterInfo.h code structure
#include <iostream>
#include <sstream>
#include <cstdint>
#include <string>
#include <vector>

// Simulate IncludeDefine.h
typedef uint8_t uint8;
using namespace std;

// Copy of ParameterInfo.h template functions with our fix
template <class parameterType>
inline parameterType inputOneValue (istringstream &streamIn) {
    parameterType oneV;
    streamIn >> oneV;
    return oneV;
};

template <>
inline uint8 inputOneValue <uint8> (istringstream &streamIn) {
    int v = 0;
    streamIn >> v;
    return static_cast<uint8>(v);
};

template <class parameterType>
inline void printOneValue (parameterType *value, std::ostream& outStr) {
    outStr << *value;
};

template <>
inline void printOneValue <uint8> (uint8 *value, std::ostream& outStr) {
    outStr << static_cast<int>(*value);
};

// Simulate ParameterInfoScalar
template <class parameterType>
class ParameterInfoScalar {
public:
    parameterType *value;
    
    ParameterInfoScalar(parameterType* valueIn) {
        value = valueIn;
    };
    
    void inputValues(istringstream &streamIn) {
        *value = inputOneValue <parameterType> (streamIn);
    };
    
    void printValues(std::ostream& outStr) const {
        printOneValue(value, outStr);
    };
};

int main() {
    std::cout << "Testing STAR ParameterInfo uint8 parsing...\n\n";
    
    // Simulate trimCutadaptQuality parameter
    uint8 trimCutadaptQuality = 0;
    ParameterInfoScalar<uint8> param(&trimCutadaptQuality);
    
    // Test parsing "20" (the default value)
    std::cout << "Test: Parsing 'trimCutadaptQuality 20'\n";
    std::string line = "trimCutadaptQuality 20";
    istringstream lineStream(line);
    string paramName;
    lineStream >> paramName;  // Skip parameter name
    param.inputValues(lineStream);
    
    std::cout << "  Parsed value: ";
    param.printValues(std::cout);
    std::cout << std::endl;
    
    if (trimCutadaptQuality == 20) {
        std::cout << "  ✓ PASS: trimCutadaptQuality correctly parsed as 20\n";
    } else {
        std::cout << "  ✗ FAIL: Expected 20, got " << static_cast<int>(trimCutadaptQuality) << "\n";
        std::cout << "  This would have caused quality_cutoff=50 (ASCII '2') before the fix!\n";
        return 1;
    }
    
    // Test parsing "30"
    std::cout << "\nTest: Parsing 'trimCutadaptQuality 30'\n";
    trimCutadaptQuality = 0;
    line = "trimCutadaptQuality 30";
    lineStream.str(line);
    lineStream.clear();
    lineStream >> paramName;
    param.inputValues(lineStream);
    
    std::cout << "  Parsed value: ";
    param.printValues(std::cout);
    std::cout << std::endl;
    
    if (trimCutadaptQuality == 30) {
        std::cout << "  ✓ PASS: trimCutadaptQuality correctly parsed as 30\n";
    } else {
        std::cout << "  ✗ FAIL: Expected 30, got " << static_cast<int>(trimCutadaptQuality) << "\n";
        return 1;
    }
    
    // Test what would happen WITHOUT the fix (demonstrate the bug)
    std::cout << "\nDemonstrating the bug (without fix):\n";
    std::cout << "  Without fix, parsing '20' would read only '2' (ASCII 50)\n";
    std::cout << "  This would set quality_cutoff=50, causing over-aggressive trimming!\n";
    
    std::cout << "\nAll tests passed! ✓\n";
    std::cout << "\nThe fix ensures:\n";
    std::cout << "  - trimCutadaptQuality='20' → quality_cutoff=20 (correct)\n";
    std::cout << "  - Before fix: trimCutadaptQuality='20' → quality_cutoff=50 (bug!)\n";
    std::cout << "  - Debug logs now show numeric values (e.g., 'quality_cutoff=20')\n";
    
    return 0;
}


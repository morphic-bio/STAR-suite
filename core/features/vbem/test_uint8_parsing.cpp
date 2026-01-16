#include <iostream>
#include <sstream>
#include <cstdint>
#include <cassert>

// Simulate the fix
typedef uint8_t uint8;

template <class parameterType>
inline parameterType inputOneValue (std::istringstream &streamIn) {
    parameterType oneV;
    streamIn >> oneV;
    return oneV;
};

// Template specialization for uint8 - THE FIX
template <>
inline uint8 inputOneValue <uint8> (std::istringstream &streamIn) {
    int v = 0;
    streamIn >> v;
    return static_cast<uint8>(v);
};

template <class parameterType>
inline void printOneValue (parameterType *value, std::ostream& outStr) {
    outStr << *value;
};

// Template specialization for uint8 - THE FIX
template <>
inline void printOneValue <uint8> (uint8 *value, std::ostream& outStr) {
    outStr << static_cast<int>(*value);
};

int main() {
    std::cout << "Testing uint8 parameter parsing fix...\n\n";
    
    // Test 1: Parse "20" as uint8
    std::istringstream stream1("20");
    uint8 value1 = inputOneValue<uint8>(stream1);
    std::cout << "Test 1: Parsing '20' -> ";
    printOneValue(&value1, std::cout);
    std::cout << std::endl;
    
    if (value1 == 20) {
        std::cout << "  ✓ PASS: Correctly parsed as 20\n";
    } else {
        std::cout << "  ✗ FAIL: Expected 20, got " << static_cast<int>(value1) << "\n";
        return 1;
    }
    
    // Test 2: Parse "50" as uint8 (should be 50, not '5')
    std::istringstream stream2("50");
    uint8 value2 = inputOneValue<uint8>(stream2);
    std::cout << "Test 2: Parsing '50' -> ";
    printOneValue(&value2, std::cout);
    std::cout << std::endl;
    
    if (value2 == 50) {
        std::cout << "  ✓ PASS: Correctly parsed as 50\n";
    } else {
        std::cout << "  ✗ FAIL: Expected 50, got " << static_cast<int>(value2) << "\n";
        return 1;
    }
    
    // Test 3: Parse "0" as uint8
    std::istringstream stream3("0");
    uint8 value3 = inputOneValue<uint8>(stream3);
    std::cout << "Test 3: Parsing '0' -> ";
    printOneValue(&value3, std::cout);
    std::cout << std::endl;
    
    if (value3 == 0) {
        std::cout << "  ✓ PASS: Correctly parsed as 0\n";
    } else {
        std::cout << "  ✗ FAIL: Expected 0, got " << static_cast<int>(value3) << "\n";
        return 1;
    }
    
    // Test 4: Parse "255" as uint8 (max value)
    std::istringstream stream4("255");
    uint8 value4 = inputOneValue<uint8>(stream4);
    std::cout << "Test 4: Parsing '255' -> ";
    printOneValue(&value4, std::cout);
    std::cout << std::endl;
    
    if (value4 == 255) {
        std::cout << "  ✓ PASS: Correctly parsed as 255\n";
    } else {
        std::cout << "  ✗ FAIL: Expected 255, got " << static_cast<int>(value4) << "\n";
        return 1;
    }
    
    // Test 5: Verify printOneValue shows numeric values
    uint8 test_val = 20;
    std::ostringstream oss;
    printOneValue(&test_val, oss);
    std::string output = oss.str();
    std::cout << "Test 5: Printing uint8(20) -> '" << output << "'\n";
    
    if (output == "20") {
        std::cout << "  ✓ PASS: Correctly prints as '20'\n";
    } else {
        std::cout << "  ✗ FAIL: Expected '20', got '" << output << "'\n";
        return 1;
    }
    
    std::cout << "\nAll tests passed! ✓\n";
    std::cout << "\nThe fix correctly:\n";
    std::cout << "  1. Parses multi-digit numbers (e.g., '20') as integers, not single characters\n";
    std::cout << "  2. Prints uint8 values as numeric strings (e.g., '20'), not character codes\n";
    
    return 0;
}


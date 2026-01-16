#include "digamma.h"
#include <cmath>

// Stable digamma implementation based on cephes-style algorithms
// Uses recurrence for x < 1, series expansion for [1, 2], asymptotic for x >= 10
// Reference: Cephes Math Library, Netlib SPECFUN, adapted for self-contained implementation

namespace {
    // Euler-Mascheroni constant
    const double EULER = 0.57721566490153286060651209008240243104215933593992;
    
    // Bernoulli numbers for asymptotic expansion
    const double B2 = 1.0 / 6.0;
    const double B4 = -1.0 / 30.0;
    const double B6 = 1.0 / 42.0;
    const double B8 = -1.0 / 30.0;
    const double B10 = 5.0 / 66.0;
    const double B12 = -691.0 / 2730.0;
    
    // Series expansion for digamma(1+x) for x in [0, 1]
    // digamma(1+x) = -EULER + sum_{k>=1} (-1)^{k+1} * zeta(k+1) * x^k
    // Using first few terms of the series
    const double PSI_SERIES[] = {
        -EULER,                    // constant term
        1.64493406684822643647,    // zeta(2) * x
        -1.20205690315959428540,   // -zeta(3) * x^2
        1.08232323371113819152,    // zeta(4) * x^3
        -1.03692775514336992633,   // -zeta(5) * x^4
        1.01734306198444913971,    // zeta(6) * x^5
        -1.00834927738192282684    // -zeta(7) * x^6
    };
}

double digamma(double x) {
    if (x <= 0.0) {
        // Invalid input - return NaN
        return std::nan("");
    }
    
    double result = 0.0;
    
    // Use recurrence relation for x < 1: digamma(x) = digamma(x+1) - 1/x
    while (x < 1.0) {
        result -= 1.0 / x;
        x += 1.0;
    }
    
    // For x >= 1, use appropriate method
    if (x >= 10.0) {
        // Asymptotic expansion for large x: digamma(x) ≈ log(x) - 1/(2x) - sum(B_{2k}/(2k*x^{2k}))
        double log_x = std::log(x);
        double inv_x = 1.0 / x;
        double inv_x2 = inv_x * inv_x;
        double inv_x4 = inv_x2 * inv_x2;
        double inv_x6 = inv_x4 * inv_x2;
        double inv_x8 = inv_x6 * inv_x2;
        double inv_x10 = inv_x8 * inv_x2;
        double inv_x12 = inv_x10 * inv_x2;
        
        result += log_x 
                  - 0.5 * inv_x
                  - B2 * inv_x2 / 2.0
                  - B4 * inv_x4 / 4.0
                  - B6 * inv_x6 / 6.0
                  - B8 * inv_x8 / 8.0
                  - B10 * inv_x10 / 10.0
                  - B12 * inv_x12 / 12.0;
    } else {
        // For 1 <= x < 10, reduce to [1, 2] using recurrence
        // digamma(x+1) = digamma(x) + 1/x
        while (x >= 2.0) {
            result += 1.0 / (x - 1.0);
            x -= 1.0;
        }
        
        // Now x is in [1, 2], use series expansion around x=1
        // digamma(1+w) where w = x - 1, w in [0, 1]
        double w = x - 1.0;
        double w_power = w;
        
        result += PSI_SERIES[0];  // constant term
        
        for (int k = 1; k < 7; ++k) {
            result += PSI_SERIES[k] * w_power;
            w_power *= w;
        }
    }
    
    return result;
}

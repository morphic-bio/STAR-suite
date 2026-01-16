#ifndef DIGAMMA_H
#define DIGAMMA_H

// Compute digamma function (psi function, derivative of log-gamma)
// digamma(x) = d/dx log(gamma(x))
// Stable, self-contained implementation based on cephes-style algorithms
// No external dependencies (no Boost)
double digamma(double x);

#endif // DIGAMMA_H

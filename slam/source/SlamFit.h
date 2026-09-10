#ifndef STAR_SLAM_FIT_H
#define STAR_SLAM_FIT_H
#include "SlamSolver.h"
#include "libem/slam_vb_overdisp.h"
struct SlamFitParameters {
    double error, conversion, phi, alpha, beta;
    bool overdispersed;
    bool operator==(const SlamFitParameters& b) const {
        return error == b.error && conversion == b.conversion && phi == b.phi &&
            alpha == b.alpha && beta == b.beta && overdispersed == b.overdispersed;
    }
};
struct SlamFit {
    double mean = 0, map = 0, sigma = 0, likelihood = 0;
    bool converged = false;
};
inline SlamFit fitSlamHistogram(const MismatchHistogram& histogram, const SlamFitParameters& p) {
    SlamFit fit;
    if (p.overdispersed) {
        const auto result = SlamVbOverdispSolver(p.error, p.conversion, p.phi, p.alpha, p.beta).solve(histogram);
        fit.mean = result.ntr_mean; fit.map = result.ntr_map;
        fit.likelihood = result.log_likelihood; fit.converged = result.converged;
    } else {
        const auto result = SlamSolver(p.error, p.conversion).solve(histogram);
        fit.mean = fit.map = result.ntr; fit.sigma = result.sigma;
        fit.likelihood = result.log_likelihood; fit.converged = result.converged;
    }
    return fit;
}
#endif

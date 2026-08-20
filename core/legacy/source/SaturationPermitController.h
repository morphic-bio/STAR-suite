#ifndef SATURATION_PERMIT_CONTROLLER_H
#define SATURATION_PERMIT_CONTROLLER_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace star {
namespace multiome {

// Deterministic saturation policy for the shared MAP/FEATURE/ATAC permit pool.
// Input transport and work decomposition are deliberately outside this class.
//
// Active domains are probed from the largest supplied work estimate to the
// smallest. If estimates are absent, the historical ATAC -> MAP -> FEATURE
// order is retained. Every unprobed live domain keeps one permit while the
// current domain receives all remaining capacity. Observed sustained occupancy
// becomes a borrowable floor. ETA is consulted only when at least one probe was
// capacity limited, and a decision never raises a domain above a known
// saturation point.
class SaturationPermitController {
 public:
  enum class Domain : uint8_t { MAP = 0, FEATURE = 1, ATAC = 2 };

  enum class Phase : uint8_t {
    PROBE_ATAC = 0,
    PROBE_MAP = 1,
    PROBE_FEATURE = 2,
    STEADY = 3
  };

  enum class Reason : uint8_t {
    NONE = 0,
    ATAC_PROBE_COMPLETE,
    MAP_PROBE_COMPLETE,
    FEATURE_PROBE_COMPLETE,
    MAP_ETA_LATE,
    MAP_ETA_EARLY,
    FEATURE_ETA_LATE,
    ATAC_ETA_LATE,
    MAP_COMPLETE,
    FEATURE_COMPLETE,
    ATAC_COMPLETE
  };

  struct WorkEstimates {
    uint64_t map = 0;
    uint64_t feature = 0;
    uint64_t atac = 0;
  };

  struct Config {
    int configuredPermits = 2;
    int fixedFeatureFloor = 0;
    int probeWindows = 2;
    int completionWindows = 2;
    bool featureActive = false;
    WorkEstimates workEstimates;
  };

  struct Observation {
    double mapOccupancy = 0.0;
    double featureOccupancy = 0.0;
    double atacOccupancy = 0.0;
    uint64_t mapUnitsDelta = 0;
    uint64_t featureUnitsDelta = 0;
    uint64_t atacUnitsDelta = 0;
    int mapInUse = 0;
    int featureInUse = 0;
    int atacInUse = 0;
    int mapWaiters = 0;
    int featureWaiters = 0;
    int atacWaiters = 0;
    double mapEtaSec = std::numeric_limits<double>::infinity();
    double featureEtaSec = std::numeric_limits<double>::infinity();
    double atacEtaSec = std::numeric_limits<double>::infinity();
    bool mapEstimateComplete = false;
    bool featureEstimateComplete = false;
    bool atacEstimateComplete = false;
  };

  struct Decision {
    Phase phase = Phase::PROBE_ATAC;
    Reason reason = Reason::NONE;
    Domain probeDomain = Domain::ATAC;
    bool floorsChanged = false;
    bool capacityLimited = false;
    bool mapSaturationKnown = false;
    bool featureSaturationKnown = false;
    bool atacSaturationKnown = false;
    int mapFloor = 1;
    int featureFloor = 0;
    int atacFloor = 1;
    int mapSaturation = 0;
    int featureSaturation = 0;
    int atacSaturation = 0;
  };

  // Backward-compatible two-domain construction. A configured FEATURE floor
  // remains fixed and is excluded from the MAP/ATAC probe budget.
  explicit SaturationPermitController(int configuredPermits,
                                      int fixedFeatureFloor = 0,
                                      int probeWindows = 2)
      : SaturationPermitController(makeLegacyConfig(
            configuredPermits, fixedFeatureFloor, probeWindows)) {}

  explicit SaturationPermitController(const Config &config)
      : configuredPermits_(std::max(config.featureActive ? 3 : 2,
                                    config.configuredPermits)),
        featureActive_(config.featureActive),
        fixedFeatureFloor_(config.featureActive
                               ? 0
                               : std::max(0, std::min(
                                     config.fixedFeatureFloor,
                                     configuredPermits_ - 2))),
        budget_(configuredPermits_ - fixedFeatureFloor_),
        probeWindows_(std::max(1, config.probeWindows)),
        completionWindows_(std::max(1, config.completionWindows)),
        estimates_(config.workEstimates) {
    floors_.fill(0);
    saturation_.fill(0);
    saturationKnown_.fill(false);
    complete_.fill(false);
    completionStreak_.fill(0);
    floors_[domainIndex(Domain::FEATURE)] = fixedFeatureFloor_;
    buildProbeOrder();
    beginProbe(0);
  }

  Decision initialDecision() const { return makeDecision(Reason::NONE, false); }

  Decision observe(const Observation &observation) {
    Reason completionReason = Reason::NONE;
    bool completionChanged = false;
    for (Domain domain : probeOrder_) {
      const size_t index = domainIndex(domain);
      if (complete_[index]) {
        continue;
      }
      const bool quiescent = estimateComplete(observation, domain) &&
          unitsDelta(observation, domain) == 0 &&
          inUse(observation, domain) == 0 && waiters(observation, domain) == 0;
      completionStreak_[index] = quiescent ? completionStreak_[index] + 1 : 0;
      if (completionStreak_[index] >= completionWindows_) {
        complete_[index] = true;
        floors_[index] = 0;
        completionChanged = true;
        if (completionReason == Reason::NONE) {
          completionReason = domainCompleteReason(domain);
        }
      }
    }

    if (phase_ != Phase::STEADY) {
      advancePastCompletedProbes();
      if (phase_ == Phase::STEADY) {
        recomputeCapacityLimited();
        return makeDecision(completionReason, completionChanged);
      }
      const Domain domain = probeOrder_[probeIndex_];
      const size_t index = domainIndex(domain);
      const bool allocationRealized =
          inUse(observation, domain) >= floors_[index] ||
          waiters(observation, domain) == 0 ||
          noOtherDomainAboveFloor(observation, domain);
      if (allocationRealized && unitsDelta(observation, domain) > 0) {
        occupancySum_ += clampOccupancy(occupancy(observation, domain));
        ++validProbeWindows_;
      }
      if (validProbeWindows_ < probeWindows_) {
        return makeDecision(completionReason, completionChanged);
      }

      saturation_[index] = saturationFromAverage(
          occupancySum_, validProbeWindows_, probeCapacity_);
      saturationKnown_[index] = saturation_[index] < probeCapacity_;
      floors_[index] = saturation_[index];
      const Reason reason = probeCompleteReason(domain);

      if (probeIndex_ + 1 < probeOrder_.size()) {
        beginProbe(probeIndex_ + 1);
      } else {
        phase_ = Phase::STEADY;
        recomputeCapacityLimited();
      }
      return makeDecision(reason, true);
    }

    if (completionChanged) {
      recomputeCapacityLimited();
      return makeDecision(completionReason, true);
    }

    if (!capacityLimited_) {
      return makeDecision(Reason::NONE, false);
    }

    std::vector<Domain> etaDomains;
    for (Domain domain : probeOrder_) {
      if (unitsDelta(observation, domain) > 0 &&
          std::isfinite(eta(observation, domain))) {
        etaDomains.push_back(domain);
      }
    }
    if (etaDomains.size() < 2) {
      return makeDecision(Reason::NONE, false);
    }

    const auto etaLess = [&observation](Domain lhs, Domain rhs) {
      return eta(observation, lhs) < eta(observation, rhs);
    };
    const Domain target = *std::max_element(
        etaDomains.begin(), etaDomains.end(), etaLess);
    std::vector<Domain> donors;
    for (Domain domain : etaDomains) {
      if (domain != target && floors_[domainIndex(domain)] > 1) {
        donors.push_back(domain);
      }
    }
    if (donors.empty()) {
      return makeDecision(Reason::NONE, false);
    }
    const Domain donor = *std::min_element(
        donors.begin(), donors.end(), etaLess);
    const double donorEta = eta(observation, donor);
    const double targetEta = eta(observation, target);

    constexpr double kEtaHysteresis = 0.10;
    if (target == donor || targetEta <= donorEta * (1.0 + kEtaHysteresis)) {
      return makeDecision(Reason::NONE, false);
    }

    const size_t targetIndex = domainIndex(target);
    const size_t donorIndex = domainIndex(donor);
    const int targetCeiling = saturationKnown_[targetIndex]
        ? saturation_[targetIndex]
        : budget_ - static_cast<int>(probeOrder_.size()) + 1;
    if (floors_[targetIndex] >= targetCeiling || floors_[donorIndex] <= 1) {
      return makeDecision(Reason::NONE, false);
    }

    ++floors_[targetIndex];
    --floors_[donorIndex];
    return makeDecision(etaReason(target), true);
  }

  static const char *domainName(Domain domain) {
    switch (domain) {
      case Domain::MAP: return "map";
      case Domain::FEATURE: return "feature";
      case Domain::ATAC:
      default: return "atac";
    }
  }

  static const char *phaseName(Phase phase) {
    switch (phase) {
      case Phase::PROBE_MAP: return "probe-map";
      case Phase::PROBE_FEATURE: return "probe-feature";
      case Phase::STEADY: return "steady";
      case Phase::PROBE_ATAC:
      default: return "probe-atac";
    }
  }

  static const char *reasonName(Reason reason) {
    switch (reason) {
      case Reason::ATAC_PROBE_COMPLETE: return "atac-probe-complete";
      case Reason::MAP_PROBE_COMPLETE: return "map-probe-complete";
      case Reason::FEATURE_PROBE_COMPLETE: return "feature-probe-complete";
      case Reason::MAP_ETA_LATE: return "map-eta-late";
      case Reason::MAP_ETA_EARLY: return "map-eta-early";
      case Reason::FEATURE_ETA_LATE: return "feature-eta-late";
      case Reason::ATAC_ETA_LATE: return "atac-eta-late";
      case Reason::MAP_COMPLETE: return "map-complete";
      case Reason::FEATURE_COMPLETE: return "feature-complete";
      case Reason::ATAC_COMPLETE: return "atac-complete";
      case Reason::NONE:
      default: return "none";
    }
  }

 private:
  static constexpr size_t kDomainCount = 3;

  static Config makeLegacyConfig(int configuredPermits,
                                 int fixedFeatureFloor,
                                 int probeWindows) {
    Config config;
    config.configuredPermits = configuredPermits;
    config.fixedFeatureFloor = fixedFeatureFloor;
    config.probeWindows = probeWindows;
    return config;
  }

  static size_t domainIndex(Domain domain) {
    return static_cast<size_t>(domain);
  }

  static int domainPriority(Domain domain) {
    switch (domain) {
      case Domain::ATAC: return 0;
      case Domain::MAP: return 1;
      case Domain::FEATURE:
      default: return 2;
    }
  }

  uint64_t estimate(Domain domain) const {
    switch (domain) {
      case Domain::MAP: return estimates_.map;
      case Domain::FEATURE: return estimates_.feature;
      case Domain::ATAC:
      default: return estimates_.atac;
    }
  }

  static double occupancy(const Observation &observation, Domain domain) {
    switch (domain) {
      case Domain::MAP: return observation.mapOccupancy;
      case Domain::FEATURE: return observation.featureOccupancy;
      case Domain::ATAC:
      default: return observation.atacOccupancy;
    }
  }

  static uint64_t unitsDelta(const Observation &observation, Domain domain) {
    switch (domain) {
      case Domain::MAP: return observation.mapUnitsDelta;
      case Domain::FEATURE: return observation.featureUnitsDelta;
      case Domain::ATAC:
      default: return observation.atacUnitsDelta;
    }
  }

  static int inUse(const Observation &observation, Domain domain) {
    switch (domain) {
      case Domain::MAP: return observation.mapInUse;
      case Domain::FEATURE: return observation.featureInUse;
      case Domain::ATAC:
      default: return observation.atacInUse;
    }
  }

  static int waiters(const Observation &observation, Domain domain) {
    switch (domain) {
      case Domain::MAP: return observation.mapWaiters;
      case Domain::FEATURE: return observation.featureWaiters;
      case Domain::ATAC:
      default: return observation.atacWaiters;
    }
  }

  static double eta(const Observation &observation, Domain domain) {
    switch (domain) {
      case Domain::MAP: return observation.mapEtaSec;
      case Domain::FEATURE: return observation.featureEtaSec;
      case Domain::ATAC:
      default: return observation.atacEtaSec;
    }
  }

  static bool estimateComplete(const Observation &observation, Domain domain) {
    switch (domain) {
      case Domain::MAP: return observation.mapEstimateComplete;
      case Domain::FEATURE: return observation.featureEstimateComplete;
      case Domain::ATAC:
      default: return observation.atacEstimateComplete;
    }
  }

  void buildProbeOrder() {
    probeOrder_.push_back(Domain::ATAC);
    probeOrder_.push_back(Domain::MAP);
    if (featureActive_) {
      probeOrder_.push_back(Domain::FEATURE);
    }
    bool anyEstimate = false;
    for (Domain domain : probeOrder_) {
      anyEstimate = anyEstimate || estimate(domain) > 0;
    }
    if (anyEstimate) {
      std::stable_sort(probeOrder_.begin(), probeOrder_.end(),
          [this](Domain lhs, Domain rhs) {
            const uint64_t lhsEstimate = estimate(lhs);
            const uint64_t rhsEstimate = estimate(rhs);
            if (lhsEstimate != rhsEstimate) {
              return lhsEstimate > rhsEstimate;
            }
            return domainPriority(lhs) < domainPriority(rhs);
          });
    }
  }

  void beginProbe(size_t index) {
    probeIndex_ = index;
    validProbeWindows_ = 0;
    occupancySum_ = 0.0;
    for (Domain domain : probeOrder_) {
      const size_t domainIndexValue = domainIndex(domain);
      if (!complete_[domainIndexValue] && saturation_[domainIndexValue] == 0) {
        floors_[domainIndexValue] = 1;
      }
    }

    const Domain probeDomain = probeOrder_[probeIndex_];
    int otherFloors = 0;
    for (Domain domain : probeOrder_) {
      if (domain != probeDomain && !complete_[domainIndex(domain)]) {
        otherFloors += floors_[domainIndex(domain)];
      }
    }
    probeCapacity_ = std::max(1, budget_ - otherFloors);
    floors_[domainIndex(probeDomain)] = probeCapacity_;
    phase_ = phaseForDomain(probeDomain);
  }

  void advancePastCompletedProbes() {
    while (phase_ != Phase::STEADY &&
           complete_[domainIndex(probeOrder_[probeIndex_])]) {
      if (probeIndex_ + 1 < probeOrder_.size()) {
        beginProbe(probeIndex_ + 1);
      } else {
        phase_ = Phase::STEADY;
      }
    }
  }

  void recomputeCapacityLimited() {
    capacityLimited_ = false;
    for (Domain domain : probeOrder_) {
      const size_t index = domainIndex(domain);
      if (!complete_[index] && !saturationKnown_[index]) {
        capacityLimited_ = true;
        return;
      }
    }
  }

  bool noOtherDomainAboveFloor(const Observation &observation,
                               Domain probeDomain) const {
    for (Domain domain : probeOrder_) {
      if (domain != probeDomain &&
          inUse(observation, domain) > floors_[domainIndex(domain)]) {
        return false;
      }
    }
    return true;
  }

  static double clampOccupancy(double value) {
    return std::max(0.0, value);
  }

  static int saturationFromAverage(double occupancySum, int windows,
                                   int ceiling) {
    if (windows <= 0) return 1;
    const int observed = static_cast<int>(
        std::ceil(occupancySum / static_cast<double>(windows)));
    return std::max(1, std::min(observed, std::max(1, ceiling)));
  }

  static Phase phaseForDomain(Domain domain) {
    switch (domain) {
      case Domain::MAP: return Phase::PROBE_MAP;
      case Domain::FEATURE: return Phase::PROBE_FEATURE;
      case Domain::ATAC:
      default: return Phase::PROBE_ATAC;
    }
  }

  static Reason probeCompleteReason(Domain domain) {
    switch (domain) {
      case Domain::MAP: return Reason::MAP_PROBE_COMPLETE;
      case Domain::FEATURE: return Reason::FEATURE_PROBE_COMPLETE;
      case Domain::ATAC:
      default: return Reason::ATAC_PROBE_COMPLETE;
    }
  }

  static Reason domainCompleteReason(Domain domain) {
    switch (domain) {
      case Domain::MAP: return Reason::MAP_COMPLETE;
      case Domain::FEATURE: return Reason::FEATURE_COMPLETE;
      case Domain::ATAC:
      default: return Reason::ATAC_COMPLETE;
    }
  }

  Reason etaReason(Domain target) const {
    if (!featureActive_ && target == Domain::ATAC) {
      return Reason::MAP_ETA_EARLY;
    }
    switch (target) {
      case Domain::MAP: return Reason::MAP_ETA_LATE;
      case Domain::FEATURE: return Reason::FEATURE_ETA_LATE;
      case Domain::ATAC:
      default: return Reason::ATAC_ETA_LATE;
    }
  }

  Decision makeDecision(Reason reason, bool changed) const {
    Decision decision;
    decision.phase = phase_;
    decision.reason = reason;
    decision.probeDomain = phase_ == Phase::STEADY
        ? probeOrder_.back()
        : probeOrder_[probeIndex_];
    decision.floorsChanged = changed;
    decision.capacityLimited = capacityLimited_;
    decision.mapSaturationKnown = saturationKnown_[domainIndex(Domain::MAP)];
    decision.featureSaturationKnown =
        saturationKnown_[domainIndex(Domain::FEATURE)];
    decision.atacSaturationKnown = saturationKnown_[domainIndex(Domain::ATAC)];
    decision.mapFloor = floors_[domainIndex(Domain::MAP)];
    decision.featureFloor = floors_[domainIndex(Domain::FEATURE)];
    decision.atacFloor = floors_[domainIndex(Domain::ATAC)];
    decision.mapSaturation = saturation_[domainIndex(Domain::MAP)];
    decision.featureSaturation = saturation_[domainIndex(Domain::FEATURE)];
    decision.atacSaturation = saturation_[domainIndex(Domain::ATAC)];
    return decision;
  }

  int configuredPermits_;
  bool featureActive_;
  int fixedFeatureFloor_;
  int budget_;
  int probeWindows_;
  int completionWindows_;
  WorkEstimates estimates_;
  std::array<int, kDomainCount> floors_{};
  std::array<int, kDomainCount> saturation_{};
  std::array<bool, kDomainCount> saturationKnown_{};
  std::array<bool, kDomainCount> complete_{};
  std::array<int, kDomainCount> completionStreak_{};
  std::vector<Domain> probeOrder_;
  size_t probeIndex_ = 0;
  Phase phase_ = Phase::PROBE_ATAC;
  int probeCapacity_ = 1;
  int validProbeWindows_ = 0;
  double occupancySum_ = 0.0;
  bool capacityLimited_ = false;
};

}  // namespace multiome
}  // namespace star

#endif

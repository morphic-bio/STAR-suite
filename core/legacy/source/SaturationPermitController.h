#ifndef SATURATION_PERMIT_CONTROLLER_H
#define SATURATION_PERMIT_CONTROLLER_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace star {
namespace multiome {

// Small, deterministic policy object for the combined MAP/ATAC permit pool.
// It deliberately learns sustained runnable occupancy rather than a noisy
// marginal-throughput curve:
//   1. give the larger ATAC arm all but one permit and observe it;
//   2. protect that observed demand and let MAP probe the remainder;
//   3. retain the learned values as borrowable floors when they fit;
//   4. only when capacity limited, use completion ETA to move one floor at a
//      time while never exceeding a learned saturation point.
//
// Floors are reservations, not partitions. ThreadControl lends unused floor
// capacity immediately, so an input/decode stall in one arm cannot idle work
// that is ready in the other arm.
class SaturationPermitController {
 public:
  enum class Phase : uint8_t { PROBE_ATAC = 0, PROBE_MAP = 1, STEADY = 2 };

  enum class Reason : uint8_t {
    NONE = 0,
    ATAC_PROBE_COMPLETE,
    MAP_PROBE_COMPLETE,
    MAP_ETA_LATE,
    MAP_ETA_EARLY,
    MAP_COMPLETE,
    ATAC_COMPLETE
  };

  struct Observation {
    double mapOccupancy = 0.0;
    double atacOccupancy = 0.0;
    uint64_t mapUnitsDelta = 0;
    uint64_t atacUnitsDelta = 0;
    int mapInUse = 0;
    int atacInUse = 0;
    int mapWaiters = 0;
    int atacWaiters = 0;
    double mapEtaSec = std::numeric_limits<double>::infinity();
    double atacEtaSec = std::numeric_limits<double>::infinity();
  };

  struct Decision {
    Phase phase = Phase::PROBE_ATAC;
    Reason reason = Reason::NONE;
    bool floorsChanged = false;
    bool capacityLimited = false;
    bool mapSaturationKnown = false;
    bool atacSaturationKnown = false;
    int mapFloor = 1;
    int atacFloor = 1;
    int mapSaturation = 0;
    int atacSaturation = 0;
  };

  explicit SaturationPermitController(int configuredPermits,
                                      int featureFloor = 0,
                                      int probeWindows = 2)
      : budget_(std::max(2, configuredPermits - std::max(0, featureFloor))),
        probeWindows_(std::max(1, probeWindows)),
        mapFloor_(1),
        atacFloor_(std::max(1, budget_ - 1)) {}

  Decision initialDecision() const { return makeDecision(Reason::NONE, false); }

  Decision observe(const Observation &observation) {
    if (phase_ == Phase::PROBE_ATAC) {
      if (observation.atacUnitsDelta > 0) {
        atacOccupancySum_ += clampOccupancy(observation.atacOccupancy);
        ++validProbeWindows_;
      }
      if (validProbeWindows_ >= probeWindows_) {
        atacSaturation_ = saturationFromAverage(
            atacOccupancySum_, validProbeWindows_, budget_ - 1);
        atacSaturationKnown_ = atacSaturation_ < budget_ - 1;
        atacFloor_ = atacSaturation_;
        mapProbeCapacity_ = std::max(1, budget_ - atacFloor_);
        mapFloor_ = mapProbeCapacity_;
        phase_ = Phase::PROBE_MAP;
        validProbeWindows_ = 0;
        mapOccupancySum_ = 0.0;
        return makeDecision(Reason::ATAC_PROBE_COMPLETE, true);
      }
      return makeDecision(Reason::NONE, false);
    }

    if (phase_ == Phase::PROBE_MAP) {
      // Floor changes are deliberately non-preemptive. Do not mistake MAP's
      // forced low occupancy for saturation while ATAC still holds permits
      // from the preceding 31-permit probe. Count a window only after MAP has
      // either reached its requested probe floor, has no waiter (natural
      // saturation), or ATAC has drained to its new learned floor.
      const bool allocationRealized =
          observation.mapInUse >= mapFloor_ || observation.mapWaiters == 0 ||
          observation.atacInUse <= atacFloor_;
      if (allocationRealized && observation.mapUnitsDelta > 0) {
        mapOccupancySum_ += clampOccupancy(observation.mapOccupancy);
        ++validProbeWindows_;
      }
      if (validProbeWindows_ >= probeWindows_) {
        mapSaturation_ = saturationFromAverage(
            mapOccupancySum_, validProbeWindows_, mapProbeCapacity_);
        mapSaturationKnown_ = mapSaturation_ < mapProbeCapacity_;
        mapFloor_ = mapSaturation_;
        // A capacity-limited probe consumed all capacity left after ATAC's
        // learned demand. ETA becomes the simple tie-breaker in steady state.
        capacityLimited_ = !mapSaturationKnown_ || !atacSaturationKnown_ ||
                           mapSaturation_ + atacSaturation_ > budget_;
        phase_ = Phase::STEADY;
        return makeDecision(Reason::MAP_PROBE_COMPLETE, true);
      }
      return makeDecision(Reason::NONE, false);
    }

    // No rate comparison is used here. When both saturation demands fit,
    // their borrowable floors are already the stable answer. Remaining work
    // matters only when the probe was capacity limited.
    // Do not infer completion from a zero-work sampling window: both input
    // providers can have multi-second decode/load gaps. Borrowable floors
    // already disappear operationally when a completed domain has no waiter,
    // so no explicit completion transition is necessary.
    if (!capacityLimited_ || observation.mapUnitsDelta == 0 ||
        observation.atacUnitsDelta == 0 ||
        !std::isfinite(observation.mapEtaSec) ||
        !std::isfinite(observation.atacEtaSec)) {
      return makeDecision(Reason::NONE, false);
    }

    // Ten-percent hysteresis is intentionally coarse. The objective is a
    // stable good-enough allocation, not chasing noisy marginal rates.
    constexpr double kEtaHysteresis = 0.10;
    if (observation.mapEtaSec >
        observation.atacEtaSec * (1.0 + kEtaHysteresis)) {
      const int mapCeiling = mapSaturationKnown_ ? mapSaturation_ : budget_ - 1;
      if (mapFloor_ < mapCeiling && atacFloor_ > 1) {
        ++mapFloor_;
        --atacFloor_;
        return makeDecision(Reason::MAP_ETA_LATE, true);
      }
    } else if (observation.mapEtaSec * (1.0 + kEtaHysteresis) <
               observation.atacEtaSec) {
      const int atacCeiling =
          atacSaturationKnown_ ? atacSaturation_ : budget_ - 1;
      if (atacFloor_ < atacCeiling && mapFloor_ > 1) {
        ++atacFloor_;
        --mapFloor_;
        return makeDecision(Reason::MAP_ETA_EARLY, true);
      }
    }
    return makeDecision(Reason::NONE, false);
  }

  static const char *phaseName(Phase phase) {
    switch (phase) {
      case Phase::PROBE_MAP: return "probe-map";
      case Phase::STEADY: return "steady";
      case Phase::PROBE_ATAC:
      default: return "probe-atac";
    }
  }

  static const char *reasonName(Reason reason) {
    switch (reason) {
      case Reason::ATAC_PROBE_COMPLETE: return "atac-probe-complete";
      case Reason::MAP_PROBE_COMPLETE: return "map-probe-complete";
      case Reason::MAP_ETA_LATE: return "map-eta-late";
      case Reason::MAP_ETA_EARLY: return "map-eta-early";
      case Reason::MAP_COMPLETE: return "map-complete";
      case Reason::ATAC_COMPLETE: return "atac-complete";
      case Reason::NONE:
      default: return "none";
    }
  }

 private:
  static double clampOccupancy(double occupancy) {
    return std::max(0.0, occupancy);
  }

  static int saturationFromAverage(double occupancySum, int windows,
                                   int ceiling) {
    if (windows <= 0) return 1;
    const int observed = static_cast<int>(
        std::ceil(occupancySum / static_cast<double>(windows)));
    return std::max(1, std::min(observed, std::max(1, ceiling)));
  }

  Decision makeDecision(Reason reason, bool changed) const {
    Decision decision;
    decision.phase = phase_;
    decision.reason = reason;
    decision.floorsChanged = changed;
    decision.capacityLimited = capacityLimited_;
    decision.mapSaturationKnown = mapSaturationKnown_;
    decision.atacSaturationKnown = atacSaturationKnown_;
    decision.mapFloor = mapFloor_;
    decision.atacFloor = atacFloor_;
    decision.mapSaturation = mapSaturation_;
    decision.atacSaturation = atacSaturation_;
    return decision;
  }

  int budget_;
  int probeWindows_;
  Phase phase_ = Phase::PROBE_ATAC;
  int validProbeWindows_ = 0;
  int mapProbeCapacity_ = 1;
  int mapFloor_;
  int atacFloor_;
  int mapSaturation_ = 0;
  int atacSaturation_ = 0;
  bool mapSaturationKnown_ = false;
  bool atacSaturationKnown_ = false;
  bool capacityLimited_ = false;
  double mapOccupancySum_ = 0.0;
  double atacOccupancySum_ = 0.0;
};

}  // namespace multiome
}  // namespace star

#endif

#ifndef DataProducts_CFOEventTiming_hh
#define DataProducts_CFOEventTiming_hh
//
// Conventions for the event timing in the CFO event record, shared by simulation
// (EventWindowMarkerProducer, CFOFragmentFromSim) and data (EventHeaderFromCFOFragment).
//
// DR marker TDC: the time from a 40 MHz clock edge to the delivery ring (DR) marker arriving
// at the CFO, in 1.25 ns (800 MHz) ticks (Mu2e-doc-55726, Mu2e-doc-4097). A valid TDC is
// 0 to clockTick/tdcTick - 1 (0 to 19). The proton bunch arrives 0 to one clock tick before
// -timeFromProtonsToDRMarker relative to the event window start, and the TDC measures where
// in that tick it is.
//
// Event duration: the CFO records the event window length in clock ticks minus one
// (CFO firmware EventWindowDataGenerator).
//

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace mu2e {
  namespace CFOEventTiming {

    constexpr double tdcTick = 1.25; // DR marker TDC tick in ns

    // Largest valid DR marker TDC for a clock tick (in ns)
    inline int maxTDC(double clockTick) {
      return static_cast<int>(std::lround(clockTick / tdcTick)) - 1;
    }

    inline bool validTDC(int tdc, double clockTick) {
      return tdc >= 0 && tdc <= maxTDC(clockTick);
    }

    // DR marker TDC for a proton bunch arriving at phaseShift (in [-clockTick, 0] ns)
    // relative to -timeFromProtonsToDRMarker
    inline int tdcFromPhase(double phaseShift, double clockTick) {
      return std::clamp(static_cast<int>(std::floor((phaseShift + clockTick) / tdcTick)), 0, maxTDC(clockTick));
    }

    // Proton bunch time relative to the event window start from a valid DR marker TDC: the TDC bin center
    inline double pbtimeFromTDC(int tdc, double timeFromProtonsToDRMarker, double clockTick) {
      return -timeFromProtonsToDRMarker - clockTick + tdcTick * (tdc + 0.5);
    }
    inline double pbtimeFromTDCError() { return tdcTick / std::sqrt(12.); }

    // Proton bunch time without a TDC: the center of the clock tick the bunch can be in
    inline double nominalPBTime(double timeFromProtonsToDRMarker, double clockTick) {
      return -timeFromProtonsToDRMarker - clockTick / 2.;
    }
    inline double nominalPBTimeError(double clockTick) { return clockTick / std::sqrt(12.); }

    // Event window length in clock ticks <-> event duration in the CFO record
    inline uint16_t recordedDuration(uint16_t windowTicks) { return windowTicks > 0 ? windowTicks - 1 : 0; }
    inline uint16_t windowTicks(uint64_t recordedDuration) { return static_cast<uint16_t>(recordedDuration + 1); }

  } // namespace CFOEventTiming
} // namespace mu2e

#endif

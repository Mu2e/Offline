// Ed Callaghan
// A fixed-length sinusoid
// February 2025

#ifndef TrackerMC_SinusoidalWireSignal_hh
#define TrackerMC_SinusoidalWireSignal_hh

// stl
#include <cmath>

// clhep
#include "CLHEP/Units/GlobalPhysicalConstants.h"

// mu2e
#include "Offline/TrackerMC/inc/AnalogWireSignal.hh"

namespace mu2e{
  class SinusoidalWireSignal: public AnalogWireSignal{
    public:
      SinusoidalWireSignal(double, double, double, double, double);
     ~SinusoidalWireSignal() = default;

      double Evaluate(double) override;

    protected:
      double _amplitude;
      double _frequency;
      double _phase;
      double _time_lo;
      double _time_hi;

    private:
      /**/
  };
} // namespace mu2e

#endif

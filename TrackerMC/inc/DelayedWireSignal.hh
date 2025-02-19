// Ed Callaghan
// A delay added to a preexisting signal
// February 2025

#ifndef TrackerMC_DelayedWireSignal_hh
#define TrackerMC_DelayedWireSignal_hh

// stl
#include <memory>
#include <vector>

// mu2e
#include "Offline/TrackerMC/inc/AnalogWireSignal.hh"

namespace mu2e{
  class DelayedWireSignal: public AnalogWireSignal{
    public:
      DelayedWireSignal(AnalogWireSignalPtr&, double);
     ~DelayedWireSignal() = default;

    protected:
      AnalogWireSignalPtr& _signal;
      double _delay;

      double evaluate_shape(double) override;

    private:
      /**/
  };
} // namespace mu2e

#endif

// Ed Callaghan
// A summable of multiple signals
// February 2025

#ifndef TrackerMC_SummedWireSignal_hh
#define TrackerMC_SummedWireSignal_hh

// stl
#include <memory>
#include <vector>

// mu2e
#include "Offline/TrackerMC/inc/AnalogWireSignal.hh"

namespace mu2e{
  class SummedWireSignal: public AnalogWireSignal{
    public:
      SummedWireSignal();
     ~SummedWireSignal() = default;

      double Evaluate(double) override;

      SummedWireSignal operator+ (const AnalogWireSignalPtr&);

    protected:
      std::vector<AnalogWireSignalPtr> _components;

    private:
      /**/
  };
} // namespace mu2e

#endif

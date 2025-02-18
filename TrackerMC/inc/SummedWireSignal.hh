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

      using SigPtr = std::shared_ptr<AnalogWireSignal>;
      SummedWireSignal operator+ (const SigPtr&);

    protected:
      std::vector<SigPtr> _components;

    private:
      /**/
  };
} // namespace mu2e

#endif

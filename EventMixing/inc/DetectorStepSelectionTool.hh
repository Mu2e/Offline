// Ed Callaghan
// Interface for selecting detector step data products to propagate after mixing
// November 2024

#ifndef EventMixing_DetectorStepSelectionTool_hh
#define EventMixing_DetectorStepSelectionTool_hh

// stl
#include <algorithm>
#include <string>

// mu2e
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"

namespace mu2e{
  class DetectorStepSelectionTool{
    public:
      DetectorStepSelectionTool() = default;
     ~DetectorStepSelectionTool() = default;

      // non-templated interface to accomodate virtual methods,
      // necessary for use as configurable tools
      virtual bool Select(const CaloShowerStep&) = 0;
      virtual bool Select(const CrvStep&)        = 0;
      virtual bool Select(const StrawGasStep&)   = 0;

    protected:
      /**/

    private:
      /**/
  };
} // namespace mu2e

#endif

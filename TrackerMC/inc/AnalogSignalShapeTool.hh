// Ed Callaghan
// Interface to produce analog signals
// February 2025

#ifndef TrackerMC_AnalogSignalShapeTool_hh
#define TrackerMC_AnalogSignalShapeTool_hh

// stl
#include <memory>

// art
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/GeneralUtilities/inc/UnaryFunction.hh"

namespace mu2e{
  using UnaryFunctionPtr = std::shared_ptr<UnaryFunction>;

  class AnalogSignalShapeTool{
    public:
      AnalogSignalShapeTool() = default;
     ~AnalogSignalShapeTool() = default;

      // produce a signal shape, which need not be constant
      virtual UnaryFunctionPtr Sample() = 0;

    protected:
      /**/

    private:
      /**/
  };
} // namespace mu2e

#endif

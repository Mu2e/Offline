// Ed Callaghan
// Interface to produce analog signals
// February 2025

#ifndef TrackerMC_AnalogWireSignalTool_hh
#define TrackerMC_AnalogWireSignalTool_hh

// art
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/TrackerMC/inc/AnalogWireSignal.hh"

namespace mu2e{
  class AnalogWireSignalTool{
    public:
      AnalogWireSignalTool() = default;
     ~AnalogWireSignalTool() = default;

      virtual AnalogWireSignalPtr Sample() = 0;

    protected:
      /**/

    private:
      /**/
  };
} // namespace mu2e

#endif

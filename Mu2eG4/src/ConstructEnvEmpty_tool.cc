// A tool to create an "empty" G4study geometry.
// The calling code in Mu2eStudyWorld already creates a "box with a
// skin" world, and sometimes one does not need any other volumes.
//
// “Perfection is achieved, not when there is nothing more to add, but
// when there is nothing left to take away.”
//   - Antoine de Saint-Exupéry
//
// Andrei Gaponenko, 2024

#include "art/Utilities/ToolMacros.h"

#include "Offline/Mu2eG4/inc/InitEnvToolBase.hh"

namespace mu2e {

  class ConstructEnvEmpty: public InitEnvToolBase {
  public:
    ConstructEnvEmpty(const fhicl::ParameterSet&) {}

    int construct(VolumeInfo const&, SimpleConfig const&) override {
      return 0;
    }
  };

}

DEFINE_ART_CLASS_TOOL(mu2e::ConstructEnvEmpty)

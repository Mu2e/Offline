// Andrei Gaponenko, 2011

#ifndef EXTMONFNALBUILDINGMAKER_HH
#define EXTMONFNALBUILDINGMAKER_HH

#include <memory>
#include <string>

#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuilding.hh"

namespace mu2e {
  class SimpleConfig;
  class ProtonBeamDump;

  class ExtMonFNALBuildingMaker {

    static ExtMonFNALBuilding::CollimatorExtMonFNAL readCollimatorExtMonFNAL(const std::string& name,
                                                                             double zLength, // along dump Z
                                                                             double angleH,
                                                                             double angleV,
                                                                             const SimpleConfig& c);

    static ExtMonFNALBuilding::FilterMagnetExtMonFNAL readFilterMagnetExtMonFNAL(const SimpleConfig& c);

  public:
    static std::auto_ptr<ExtMonFNALBuilding> make(const SimpleConfig& config, const ProtonBeamDump& dump);
  };
}

#endif/*EXTMONFNALBUILDINGMAKER_HH*/

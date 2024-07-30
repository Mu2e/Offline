// Andrei Gaponenko, 2011

#ifndef EXTMONFNALBUILDINGMAKER_HH
#define EXTMONFNALBUILDINGMAKER_HH

#include <memory>
#include <string>

#include "Offline/Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"

namespace mu2e {
  class SimpleConfig;
  class ProtonBeamDump;

  class ExtMonFNALBuildingMaker {

    static ExtMonFNALBuilding::CollimatorExtMonFNAL readCollimatorExtMonFNAL(const std::string& name,
                                                                             double angleH,
                                                                             double angleV,
                                                                             const SimpleConfig& c);

  public:
    static std::unique_ptr<ExtMonFNALBuilding> make(const SimpleConfig& config,
                                                    const Mu2eHall& hall,
                                                    const ProtonBeamDump& dump);
  };
}

#endif/*EXTMONFNALBUILDINGMAKER_HH*/

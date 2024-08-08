// Andrei Gaponenko, 2011

#ifndef EXTMONFNALFILTERMAKER_HH
#define EXTMONFNALFILTERMAKER_HH

#include <string>

#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALFilter.hh"

namespace mu2e {
  class SimpleConfig;
  class ProtonBeamDump;
  class ExtMonFNALBuilding;

  class ExtMonFNALFilterMaker {

    static ExtMonFNALCollimator readExtMonFNALCollimator(const std::string& name,
                                                         const SimpleConfig& c);

    static void positionCollimatorAbsolute(ExtMonFNALFilter *filter,
                                           const std::string& prefix,
                                           const ProtonBeamDump& dump,
                                           const SimpleConfig& c);

    static void positionEntranceCollimatorRelative(ExtMonFNALCollimator *col,
                                                   const ProtonBeamDump& dump,
                                                   const SimpleConfig& c);

    static void computeChannelTail(ExtMonFNALFilter *filter,
                                   const ExtMonFNALBuilding& emfb,
                                   const ProtonBeamDump& dump,
                                   const SimpleConfig& c);

  public:
    static ExtMonFNALFilter read(const SimpleConfig& c,
                                 const ExtMonFNALBuilding& emfb,
                                 const ProtonBeamDump& dump);
  };
}

#endif/*EXTMONFNALFILTERMAKER_HH*/

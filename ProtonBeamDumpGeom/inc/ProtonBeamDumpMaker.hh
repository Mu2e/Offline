// Andrei Gaponenko, 2011

#ifndef PROTONBEAMDUMPMAKER_HH
#define PROTONBEAMDUMPMAKER_HH

#include <memory>
#include <string>

#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
namespace mu2e { class SimpleConfig; }

namespace mu2e {
  class ProtonBeamDumpMaker {
    static ProtonBeamDump::CollimatorExtMonFNAL readCollimatorExtMonFNAL(const std::string& name, double angleH, double angleV, const SimpleConfig& c);
    static ProtonBeamDump::FilterMagnetExtMonFNAL readFilterMagnetExtMonFNAL(const SimpleConfig& c);
  public:
    static std::auto_ptr<ProtonBeamDump> make(const SimpleConfig& config);
  };
}

#endif/*PROTONBEAMDUMPMAKER_HH*/

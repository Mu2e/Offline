// Andrei Gaponenko, 2011

#ifndef PROTONBEAMDUMPMAKER_HH
#define PROTONBEAMDUMPMAKER_HH

#include <memory>
#include <string>

#include "GeometryService/inc/ProtonBeamDump.hh"
namespace mu2e { class SimpleConfig; }

namespace mu2e {
  class ProtonBeamDumpMaker {
    std::auto_ptr<ProtonBeamDump> m_det;

    ProtonBeamDump::CollimatorExtMonFNAL readCollimatorExtMonFNAL(const std::string& name, double angleH, double angleV, const SimpleConfig& c);
    ProtonBeamDump::FilterMagnetExtMonFNAL readFilterMagnetExtMonFNAL(const SimpleConfig& c);
  public:
    explicit ProtonBeamDumpMaker(const SimpleConfig& config);

    // interface to GeometryService
    std::auto_ptr<ProtonBeamDump> getPtr() { return m_det; }
  };
}

#endif/*PROTONBEAMDUMPMAKER_HH*/

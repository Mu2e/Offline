// Andrei Gaponenko, 2011

#ifndef PROTONBEAMDUMPMAKER_HH
#define PROTONBEAMDUMPMAKER_HH

#include <memory>
#include <string>

#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
namespace mu2e { class SimpleConfig; }

namespace mu2e {
  class ProtonBeamDumpMaker {
  public:
    static std::auto_ptr<ProtonBeamDump> make(const SimpleConfig& config,
                                              double frontShieldingYmin,
                                              double frontShieldingYmax);
  };
}

#endif/*PROTONBEAMDUMPMAKER_HH*/

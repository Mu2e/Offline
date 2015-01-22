// Andrei Gaponenko, 2011

#ifndef PROTONBEAMDUMPMAKER_HH
#define PROTONBEAMDUMPMAKER_HH

#include <memory>
#include <string>

#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
namespace mu2e { class SimpleConfig; }

namespace mu2e {
  class ProtonBeamDumpMaker {
  public:
    static std::unique_ptr<ProtonBeamDump> make(const SimpleConfig& config,
                                                const Mu2eHall& hall);
  };
}

#endif/*PROTONBEAMDUMPMAKER_HH*/

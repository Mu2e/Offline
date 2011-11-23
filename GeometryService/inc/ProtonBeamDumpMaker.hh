// Andrei Gaponenko, 2011

#ifndef PROTONBEAMDUMPMAKER_HH
#define PROTONBEAMDUMPMAKER_HH

#include <memory>

namespace mu2e { class SimpleConfig; }
namespace mu2e { class ProtonBeamDump; }

namespace mu2e {
  class ProtonBeamDumpMaker {
    std::auto_ptr<ProtonBeamDump> m_det;
  public:
    explicit ProtonBeamDumpMaker(const SimpleConfig& config);
    
    // interface to GeometryService
    std::auto_ptr<ProtonBeamDump> getPtr() { return m_det; }
  };
}

#endif/*PROTONBEAMDUMPMAKER_HH*/

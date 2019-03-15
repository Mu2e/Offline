//
// structs used to record GenParticle information in TTrees.
// All momenta are in units of MeV/c, time in nsec WRT when the proton bunch pulse peak hits the production target,
// positions are in mm WRT the center of the tracker.
// Andy Edmonds (March 2019)
// 
#ifndef GenInfo_HH
#define GenInfo_HH
#include "DataProducts/inc/XYZVec.hh"
#include "TrkDiag/inc/helixpar.hh"
#include "Rtypes.h"
namespace mu2e
{

// general info about the gen particle which was simulated
  struct GenInfo {
    Int_t _pdg, _gen; // true PDG code, generator code
    Float_t _time;  // time of this step
    XYZVec _mom;   // momentum at the start of this step
    XYZVec _pos;  // particle position at the start of this step
    GenInfo() { reset(); }
    void reset() { _pdg = _gen = -1; _time = -1.0; _mom=XYZVec(); _pos = XYZVec(); }
    static std::string leafnames() { static std::string leaves; 
      leaves = std::string("pdg/I:gen/I:t0/F:")+Geom::XYZnames("mom") + std::string(":") + Geom::XYZnames("pos");
      return leaves;
    }
  };  
}
#endif


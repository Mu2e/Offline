#ifndef MCDataProducts_PrimaryParticle_hh
#define MCDataProducts_PrimaryParticle_hh
//
//  Define the primary particle
// Original author: David Brown (LBNL) Jan 2019
//
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include <vector>
namespace mu2e {
  class PrimaryParticle {
    public:
      typedef std::vector<art::Ptr<SimParticle> > SPPV;
      PrimaryParticle() {}
      PrimaryParticle(SPPV const& simps) : _simps(simps) {}
      SPPV const& primarySimParticles() const { return _simps; }
      ProcessCode primaryProcess() const { return _simps.size() >0 ? _simps.front()->creationCode(): ProcessCode(ProcessCode::NoProcess); }
      SPPV& modifySimParticles() { return _simps; } // needed for compression
    private:
      SPPV _simps; // associated SimParticles (can be >1).  All must have the same creation code
  };

  inline std::ostream& operator<<(std::ostream& ost, const PrimaryParticle& pp ){
    ost << "Primary process " << pp.primaryProcess();
    return ost;
  }


  // don't create a vector of these: there should only ever be one/event!
}
#endif /* MCDataProducts_PrimaryParticle_hh */

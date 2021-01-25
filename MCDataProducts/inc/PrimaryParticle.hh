#ifndef MCDataProducts_PrimaryParticle_hh
#define MCDataProducts_PrimaryParticle_hh
//
//  Define the primary particle
// Original author: David Brown (LBNL) Jan 2019
//
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include <vector>
namespace mu2e {
  class PrimaryParticle {
    public:
      typedef std::vector<art::Ptr<SimParticle> > SPPV;
      PrimaryParticle() {}
      PrimaryParticle(PrimaryParticle const& other ) : _genp(other._genp), _simps(other._simps)
      {}
      PrimaryParticle(GenParticle const& genp,SPPV const& simps) :
	_genp(genp), _simps(simps) {}
      GenParticle const& primary () const { return _genp; }
      SPPV const& primarySimParticles() const { return _simps; }
      SPPV& modifySimParticles() { return _simps; }

    private:
      GenParticle _genp; // primary particle 
      SPPV _simps; // associated SimParticles (can be >1)
// should also record the VD steps associated with the primaries TODO
  };

  inline std::ostream& operator<<(std::ostream& ost, const PrimaryParticle& pp ){
    ost << pp.primary();
    return ost;
  }
 

  // don't create a vector of these: there should only ever be one/event!
}
#endif /* MCDataProducts_PrimaryParticle_hh */

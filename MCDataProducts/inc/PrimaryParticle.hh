#ifndef MCDataProducts_PrimaryParticle_hh
#define MCDataProducts_PrimaryParticle_hh
//
//  Define the primary particle
// Original author: David Brown (LBNL) Jan 2019
//
#include "MCDataProducts/inc/SimParticle.hh"
#include <vector>
namespace mu2e {
  class PrimaryParticle {
    public:
      typedef std::vector<art::Ptr<SimParticle> > SPPV;
      PrimaryParticle() {}
      PrimaryParticle(PrimaryParticle const& other ) : _simps(other._simps)
      { // check these all have the same creation code
	auto isimp = _simps.begin()++;
	while(isimp != _simps.end()){
	  if((*isimp)->creationCode() != _simps.front()->creationCode())
	    throw cet::exception("Simulation")<<"PrimaryParticle: creation codes don't match" << std::endl;
	  ++isimp;
	}
	_pcode = _simps.front()->creationCode().id();
      }
      PrimaryParticle(SPPV const& simps) : _simps(simps) {}
      SPPV const& primarySimParticles() const { return _simps; }
      ProcessCode const& primaryProcess() const { return _pcode; }
      SPPV& modifySimParticles() { return _simps; } // needed for compression
    private:
      SPPV _simps; // associated SimParticles (can be >1).  All must have the same creation code
      ProcessCode _pcode; // process code of this primary
  };

  inline std::ostream& operator<<(std::ostream& ost, const PrimaryParticle& pp ){
    ost << "Primary process " << pp.primaryProcess();
    return ost;
  }
 

  // don't create a vector of these: there should only ever be one/event!
}
#endif /* MCDataProducts_PrimaryParticle_hh */

#ifndef Mu2eUtilities_checkSimParticleCollection_hh
#define Mu2eUtilities_checkSimParticleCollection_hh

//
// Within a SimParticleColleciton, check that all mother/daughter pointers are self-consistent.
//
// Contact person Rob Kutschke
//

#include "Offline/MCDataProducts/inc/SimParticle.hh"

namespace mu2e{

  bool checkSimParticleCollection ( SimParticleCollection const& out, bool doThrow = false );

}

#endif /* Mu2eUtilities_checkSimParticleCollection_hh */

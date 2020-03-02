#ifndef Mu2eUtilities_checkSimParticleCollection_hh
#define Mu2eUtilities_checkSimParticleCollection_hh

//
// Within a SimParticleColleciton, check that all mother/daughter pointers are self-consistent.
//
// $Id: checkSimParticleCollection.hh,v 1.1 2011/12/16 23:13:50 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/12/16 23:13:50 $
//
// Contact person Rob Kutschke
//

#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace mu2e{

  bool checkSimParticleCollection ( SimParticleCollection const& out, bool doThrow = false );

}

#endif /* Mu2eUtilities_checkSimParticleCollection_hh */

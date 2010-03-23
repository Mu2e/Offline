#ifndef ToyDP_SimParticleCollection_hh
#define ToyDP_SimParticleCollection_hh

//
// Define a type for a collection of SimParticle objects.
//
// $Id: SimParticleCollection.hh,v 1.1 2010/03/23 20:34:30 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/03/23 20:34:30 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/SimParticle.hh"

namespace mu2e {
   typedef std::vector<mu2e::SimParticle> SimParticleCollection;
}

#endif

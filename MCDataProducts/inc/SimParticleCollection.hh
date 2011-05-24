#ifndef MCDataProducts_SimParticleCollection_hh
#define MCDataProducts_SimParticleCollection_hh

//
// Define a type for a collection of SimParticle objects.
//
// $Id: SimParticleCollection.hh,v 1.2 2011/05/24 20:03:31 wb Exp $
// $Author: wb $
// $Date: 2011/05/24 20:03:31 $
//
// Original author Rob Kutschke
//

#include "MCDataProducts/inc/SimParticle.hh"
#include "cetlib/map_vector.h"

namespace mu2e {
   typedef cet::map_vector<mu2e::SimParticle> SimParticleCollection;
}

#endif /* MCDataProducts_SimParticleCollection_hh */

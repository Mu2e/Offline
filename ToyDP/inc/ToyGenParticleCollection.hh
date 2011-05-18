#ifndef ToyDP_ToyGenParticleCollection_hh
#define ToyDP_ToyGenParticleCollection_hh

//
// Define a type for a collection of ToyGenParticle.
//
// $Id: ToyGenParticleCollection.hh,v 1.3 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/ToyGenParticle.hh"

namespace mu2e {
   typedef std::vector<mu2e::ToyGenParticle> ToyGenParticleCollection;
}

#endif /* ToyDP_ToyGenParticleCollection_hh */

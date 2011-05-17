#ifndef ToyDP_ToyGenParticleCollection_hh
#define ToyDP_ToyGenParticleCollection_hh

//
// Define a type for a collection of ToyGenParticle.
//
// $Id: ToyGenParticleCollection.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/ToyGenParticle.hh"

namespace mu2e {
   typedef std::vector<mu2e::ToyGenParticle> ToyGenParticleCollection;
}

#endif /* ToyDP_ToyGenParticleCollection_hh */

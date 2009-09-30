#ifndef ToyDP_ToyGenParticleCollection_hh
#define ToyDP_ToyGenParticleCollection_hh

//
// Define a type for a collection of ToyGenParticle.
//
// $Id: ToyGenParticleCollection.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/ToyGenParticle.hh"

namespace mu2e {
   typedef std::vector<mu2e::ToyGenParticle> ToyGenParticleCollection;
}

#endif

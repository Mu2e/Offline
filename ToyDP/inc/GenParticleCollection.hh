#ifndef ToyDP_GenParticleCollection_hh
#define ToyDP_GenParticleCollection_hh

//
// Define a type for a collection of GenParticle.
//
// $Id: GenParticleCollection.hh,v 1.1 2011/05/18 22:44:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/18 22:44:23 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/GenParticle.hh"

namespace mu2e {
   typedef std::vector<mu2e::GenParticle> GenParticleCollection;
}

#endif /* ToyDP_GenParticleCollection_hh */

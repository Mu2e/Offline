#ifndef MCDataProducts_GenParticleCollection_hh
#define MCDataProducts_GenParticleCollection_hh

//
// Define a type for a collection of GenParticle.
//
// $Id: GenParticleCollection.hh,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "MCDataProducts/inc/GenParticle.hh"

namespace mu2e {
   typedef std::vector<mu2e::GenParticle> GenParticleCollection;
}

#endif /* MCDataProducts_GenParticleCollection_hh */

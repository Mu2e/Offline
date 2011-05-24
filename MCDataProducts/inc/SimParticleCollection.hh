#ifndef MCDataProducts_SimParticleCollection_hh
#define MCDataProducts_SimParticleCollection_hh

//
// Define a type for a collection of SimParticle objects.
//
// $Id: SimParticleCollection.hh,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author Rob Kutschke
//

#include "GeneralUtilities/inc/MapVector.hh"
#include "MCDataProducts/inc/SimParticle.hh"

namespace mu2e {
   typedef MapVector<mu2e::SimParticle> SimParticleCollection;
}

#endif /* MCDataProducts_SimParticleCollection_hh */

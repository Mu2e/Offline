#ifndef ToyDP_SimParticleCollection_hh
#define ToyDP_SimParticleCollection_hh

//
// Define a type for a collection of SimParticle objects.
//
// $Id: SimParticleCollection.hh,v 1.2 2010/11/08 23:51:55 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/11/08 23:51:55 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "GeneralUtilities/inc/MapVector.hh"
#include "ToyDP/inc/SimParticle.hh"

namespace mu2e {
   typedef MapVector<mu2e::SimParticle> SimParticleCollection;
}

#endif

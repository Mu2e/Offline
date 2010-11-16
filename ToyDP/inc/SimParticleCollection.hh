#ifndef ToyDP_SimParticleCollection_hh
#define ToyDP_SimParticleCollection_hh

//
// Define a type for a collection of SimParticle objects.
//
// $Id: SimParticleCollection.hh,v 1.3 2010/11/16 21:18:10 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/11/16 21:18:10 $
//
// Original author Rob Kutschke
//

#include "GeneralUtilities/inc/MapVector.hh"
#include "ToyDP/inc/SimParticle.hh"

namespace mu2e {
   typedef MapVector<mu2e::SimParticle> SimParticleCollection;
}

#endif

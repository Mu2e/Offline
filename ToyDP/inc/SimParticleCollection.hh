#ifndef ToyDP_SimParticleCollection_hh
#define ToyDP_SimParticleCollection_hh

//
// Define a type for a collection of SimParticle objects.
//
// $Id: SimParticleCollection.hh,v 1.4 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

#include "GeneralUtilities/inc/MapVector.hh"
#include "ToyDP/inc/SimParticle.hh"

namespace mu2e {
   typedef MapVector<mu2e::SimParticle> SimParticleCollection;
}

#endif /* ToyDP_SimParticleCollection_hh */

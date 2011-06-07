#ifndef Mu2eG4_finalizeSimParticle_hh
#define Mu2eG4_finalizeSimParticle_hh
//
// Set the art::Ptr<GenParticle> and art::Ptr<SimParticle> data members of a SimParticleCollection.
//
// $Id: finalizeSimParticle.hh,v 1.1 2011/06/07 21:56:50 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/07 21:56:50 $
//
// Original author Rob Kutschke
//

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Common/OrphanHandle.h"

#include <vector>

namespace mu2e {

  void finalizeSimParticle ( SimParticleCollection& v,
			     art::Handle<GenParticleCollection>& genshandle,
			     art::OrphanHandle<SimParticleCollection>& simshandle
			     );

}

#endif /* Mu2eG4_finalizeSimParticle_hh */

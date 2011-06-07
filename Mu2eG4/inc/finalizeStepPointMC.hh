#ifndef Mu2eG4_finalizeStepPointMC_hh
#define Mu2eG4_finalizeStepPointMC_hh
//
// Set the Ptr<SimParticle> data members of a StepPointMCCollection.
//
// $Id: finalizeStepPointMC.hh,v 1.1 2011/06/07 21:56:50 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/07 21:56:50 $
//
// Original author Rob Kutschke
//

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include "art/Persistency/Common/OrphanHandle.h"

#include <vector>

namespace mu2e {

  void finalizeStepPointMC ( StepPointMCCollection& v,
			     art::OrphanHandle<SimParticleCollection>& handle
			     );

}

#endif /* Mu2eG4_finalizeStepPointMC_hh */

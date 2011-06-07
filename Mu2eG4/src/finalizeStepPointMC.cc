//
// Set the Ptr<SimParticle> data members of a StepPointMCCollection.
//
// $Id: finalizeStepPointMC.cc,v 1.1 2011/06/07 21:56:50 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/07 21:56:50 $
//
// Original author Rob Kutschke
//

#include "Mu2eG4/inc/finalizeStepPointMC.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include "art/Persistency/Common/OrphanHandle.h"

#include <vector>

using namespace std;

namespace mu2e {

  void finalizeStepPointMC ( StepPointMCCollection& v,
			     art::OrphanHandle<SimParticleCollection>& handle
			     ){
    for ( StepPointMCCollection::iterator i=v.begin();
	  i !=v.end(); ++i ){
      StepPointMC& step = *i;
      step.setSimParticlePtr( handle );
    }
  }
 
}

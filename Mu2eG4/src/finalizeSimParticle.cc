//
// Set the art::Ptr<GenParticle> and art::Ptr<SimParticle> data members of a SimParticleCollection.
//
// $Id: finalizeSimParticle.cc,v 1.1 2011/06/07 21:56:50 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/07 21:56:50 $
//
// Original author Rob Kutschke
//

#include "Mu2eG4/inc/finalizeSimParticle.hh"

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include "art/Persistency/Common/OrphanHandle.h"

#include <vector>

using namespace std;

namespace mu2e {

  void finalizeSimParticle ( SimParticleCollection& v,
			     art::Handle<GenParticleCollection>& gensHandle,
			     art::OrphanHandle<SimParticleCollection>& simsHandle
			     ){

    for ( SimParticleCollection::iterator i=v.begin();
	  i !=v.end(); ++i ){
      SimParticle& sim = i->second;

      sim.setParentPtr     ( art::Ptr<SimParticle>( simsHandle, sim.getParentId() ) );
      sim.setGenParticlePtr( art::Ptr<GenParticle>( gensHandle, sim.getGenIndex() ) );

      // Make a temporary std::vector<art::Ptr<SimParticle> >.
      std::vector<SimParticleCollection::key_type> const & dauIds = sim.daughterIds();
      std::vector<art::Ptr<SimParticle> > dauSim;
      dauSim.reserve( dauIds.size() );
      for ( std::vector<SimParticleCollection::key_type>::const_iterator i=dauIds.begin();
	    i != dauIds.end(); ++i ){
	dauSim.push_back( art::Ptr<SimParticle>(simsHandle, i->asInt() ) );
      }

      sim.setDaughterPtrs( dauSim );
    }
  }
 
}

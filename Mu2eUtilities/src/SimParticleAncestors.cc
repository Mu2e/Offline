//
// Start with a SimParticle and trace its ancestry back to a generated particle.
//  
// $Id: SimParticleAncestors.cc,v 1.1 2010/11/30 23:50:22 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/11/30 23:50:22 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "Mu2eUtilities/inc/SimParticleAncestors.hh"

// Framework includes
#include "FWCore/Utilities/interface/Exception.h"

using namespace std;

namespace mu2e{

  SimParticleAncestors::SimParticleAncestors( key_type trackId,
                                              SimParticleCollection const& sims,
                                              ToyGenParticleCollection const & gens,
                                              int maxDepth):
    _sim(&sims.at(trackId)),
    _depth(0),
    _maxDepth(maxDepth){

    construct ( sims, gens);

  }

  SimParticleAncestors::SimParticleAncestors( SimParticle const& sim,
                                              SimParticleCollection const& sims,
                                              ToyGenParticleCollection const & gens,
                                              int maxDepth):
    _sim(&sim),
    _depth(0),
    _maxDepth(maxDepth){

    construct ( sims, gens);

  }


  void SimParticleAncestors::construct( SimParticleCollection const& sims,
                                        ToyGenParticleCollection const & gens){

    SimParticle const * s(_sim);
    while ( s->hasParent()){
      s = sims.findOrNull( s->parentId() );
      if ( ++_depth > _maxDepth ){
        throw cms::Exception("LIMIT")
          << "SimParticleAncestors did not find the generated particle after "
          << _depth
          << " generations.";
      }
    }
    _sim0 = s;
    _gen0 = &gens.at(s->generatorIndex());
  }
  

} // end namespace mu2e

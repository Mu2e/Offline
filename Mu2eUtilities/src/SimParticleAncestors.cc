//
// Start with a SimParticle and trace its ancestry back to a generated particle.
//
// $Id: SimParticleAncestors.cc,v 1.5 2011/05/24 20:03:31 wb Exp $
// $Author: wb $
// $Date: 2011/05/24 20:03:31 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "Mu2eUtilities/inc/SimParticleAncestors.hh"

// Framework includes
#include "cetlib_except/exception.h"

using namespace std;

namespace mu2e{

  SimParticleAncestors::SimParticleAncestors( key_type trackId,
                                              SimParticleCollection const& sims,
                                              GenParticleCollection const & gens,
                                              int maxDepth):
    _sim(&sims.at(trackId)),
    _depth(0),
    _maxDepth(maxDepth){

    construct ( sims, gens);

  }

  SimParticleAncestors::SimParticleAncestors( SimParticle const& sim,
                                              SimParticleCollection const& sims,
                                              GenParticleCollection const & gens,
                                              int maxDepth):
    _sim(&sim),
    _depth(0),
    _maxDepth(maxDepth){

    construct ( sims, gens);

  }


  void SimParticleAncestors::construct( SimParticleCollection const& sims,
                                        GenParticleCollection const & gens){

    SimParticle const * s(_sim);
    while ( s->hasParent()){
      s = sims.getOrNull( s->parentId() );
      if ( ++_depth > _maxDepth ){
        throw cet::exception("LIMIT")
          << "SimParticleAncestors did not find the generated particle after "
          << _depth
          << " generations.";
      }
    }
    _sim0 = s;
    _gen0 = &gens.at(s->generatorIndex());
  }


} // end namespace mu2e

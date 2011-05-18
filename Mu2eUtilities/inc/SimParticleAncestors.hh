#ifndef Mu2eUtilities_SimParticleAncestors_hh
#define Mu2eUtilities_SimParticleAncestors_hh
//
// Start with a SimParticle and trace its ancestry back to a generated particle.
//
// $Id: SimParticleAncestors.hh,v 1.3 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) This class is not designed to be persisted because it has bare pointers.
//

#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"

namespace mu2e{


  class SimParticleAncestors{

  public:
    typedef SimParticleCollection::key_type key_type;

    SimParticleAncestors( key_type key,
                          SimParticleCollection const& sims,
                          ToyGenParticleCollection const & gens,
                          int maxDepth=100);

    SimParticleAncestors( SimParticle const& sim,
                          SimParticleCollection const& sims,
                          ToyGenParticleCollection const & gens,
                          int maxDepth=100);

    // Compiler generated code is Ok for:
    //  d'tor, copy c'tor assignment operator.

    SimParticle    const& sim()         const {return *_sim;}
    SimParticle    const& originalSim() const {return *_sim0;}
    ToyGenParticle const& originalGen() const {return *_gen0;}
    int                   depth()       const {return _depth;}

  private:

    // Non-owning pointers to the input particle.
    SimParticle    const* _sim;

    // Non-owning pointers the generated particle that is in the ancestry list
    // of the input particle; also the sim particle that comes
    // directly from the generated particle.
    SimParticle    const* _sim0;
    ToyGenParticle const* _gen0;

    // Number of generations from the generated particle to this one.
    int _depth;

    // Limit on number of generations in search for the generated particle.
    int _maxDepth;

    // A helper function to do the work that is common to the several constructors.
    void construct( SimParticleCollection const& sims,
                    ToyGenParticleCollection const & gens);
  };

} // end namespace mu2e


#endif /* Mu2eUtilities_SimParticleAncestors_hh */

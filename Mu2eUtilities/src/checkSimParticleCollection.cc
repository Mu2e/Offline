//
// Within a SimParticleColleciton, check that all mother/daughter pointers are self-consistent.
//
//
// Contact person Rob Kutschke

#include <cstddef>
#include <utility>
#include <vector>

#include "cetlib/map_vector.h"
#include "cetlib_except/exception.h"
// art includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "Mu2eUtilities/inc/checkSimParticleCollection.hh"

using namespace std;

namespace mu2e {

  bool checkSimParticleCollection ( SimParticleCollection const& out, bool doThrow){

    bool ok(true);

    for ( SimParticleCollection::const_iterator i=out.begin(), e=out.end(); i!=e; ++i ){

      // The next particle to look at.
      SimParticle const& sim = i->second;

      // Check that all daughters see this particle as their mother.
      std::vector<cet::map_vector_key> const& dau = sim.daughterIds();
      for ( size_t j=0; j<dau.size(); ++j ){

        SimParticle const* daughter  = out.getOrNull(cet::map_vector_key(dau[j]));
        if ( daughter == 0 ) {
          ok = false;
          mf::LogError("DATA")
            << "CheckSimsConsisenty:: daughter is not present:"
            << " Particle: "    << i->first
            << " Daughter index and id: " << j << " " << dau[j]
            << "\n";
        } else {

          cet::map_vector_key parentId = daughter->parentId();
          if ( parentId != sim.id() ){
            ok = false;
            mf::LogError("DATA")
              << "CheckSimsConsisenty:: daughter does not point back to mother."
              << " Particle: "    << i->first
              << " Daughter index and id: " << j << " " << dau[j]
              << " Parent:      " << parentId
              << "\n";
          }
        }

      } // end loop over daughters

      // Check that this particle is in the list of its parent's daughters.
      if ( sim.hasParent() ){
        cet::map_vector_key parentId = sim.parentId();

        SimParticle const* parent = out.getOrNull(parentId);
        if ( parent == 0 ){
          ok = false;
          mf::LogError("DATA")
            << "CheckSimsConsisenty:: parent is not present:"
            << " Particle: "    << i->first
            << " ParentId: "    << parentId
            << "\n";
        } else {
          std::vector<cet::map_vector_key> const& mdau = parent->daughterIds();

          bool inList(false);
          for ( size_t j=0; j<mdau.size(); ++j ){
            if ( cet::map_vector_key(mdau.at(j)) == sim.id() ){
              inList = true;
              break;
            }
          }
          if ( !inList ){
            ok = false;
            mf::LogError("DATA")
              << "CheckSimsConsisenty:: mother does not point back to daughter "
              << " Particle: "    << i->first
              << " Parent:      " << parentId
              << "\n";
          }
        }

      } // end hasParent

    } // end loop over all SimParticles

    if ( !ok && doThrow ){
      throw cet::exception("DATA") << "CheckSimsConsisenty:: errors discovered; see above.  Aborting now.\n";
    }

    return ok;

  } // end checkSimsConsistency

} // end namespace mu2e

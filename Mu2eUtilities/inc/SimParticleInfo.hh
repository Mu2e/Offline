#ifndef Mu2eUtilities_SimParticleInfo_HH
#define Mu2eUtilities_SimParticleInfo_HH
//
// Information about one SimParticle and all StrawHits that are
// associated with hit.  This is a building block of the
// the class SimParticlesWithHits.
//
// $Id: SimParticleInfo.hh,v 1.2 2010/12/09 15:30:07 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/12/09 15:30:07 $
//
// Original author Rob Kutschke.
//

// C++ includes.
#include <vector>

// Mu2e includes.
#include "Mu2eUtilities/inc/StrawHitMCInfo.hh"
#include "ToyDP/inc/SimParticleCollection.hh"

namespace edm{
  class Event;
}


using namespace std;

namespace mu2e {
  
  class SimParticleInfo{
  public:
    typedef SimParticleCollection::key_type key_type;

    SimParticleInfo():_simId(-1){}

    SimParticleInfo( key_type simId,
                     SimParticle const& simParticle,
                     edm::Event const& event):
      _simId(simId),
      _simParticle(&simParticle),
      _event(&event),
      _hitInfos()
    {
    }

    key_type id() const { return _simId; }
    SimParticle const& simParticle() const { return *_simParticle; }

    vector<StrawHitMCInfo>&      strawHitInfos()       { return _hitInfos; }
    vector<StrawHitMCInfo>const& strawHitInfos() const { return _hitInfos; }

    // Compiler generated code is Ok for:
    //  d'tor, copy c'tor assignment operator.

  private:

    // ID of this particle in the SimParticleCollection.
    key_type _simId;

    // Pointer to the SimParticle
    SimParticle const* _simParticle;

    // The event in which this information is found.
    edm::Event const* _event;

    // Vector of information about the StrawHits to which this track contributed.
    vector<StrawHitMCInfo>  _hitInfos;

  };

} // namespace mu2e

#endif

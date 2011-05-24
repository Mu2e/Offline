#ifndef Mu2eUtilities_SimParticleInfo_hh
#define Mu2eUtilities_SimParticleInfo_hh
//
// Information about one SimParticle and all StrawHits that are
// associated with hit.  This is a building block of the
// the class SimParticlesWithHits.
//
// $Id: SimParticleInfo.hh,v 1.8 2011/05/24 17:19:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:19:03 $
//
// Original author Rob Kutschke.
//

// C++ includes.
#include <vector>

// Mu2e includes.
#include "Mu2eUtilities/inc/StrawHitMCInfo.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace art{
  class Event;
}


namespace mu2e {

  class SimParticleInfo{

    // This class should only ever be created within SimParticlesWithHits.
    friend class SimParticlesWithHits;

  public:
    typedef SimParticleCollection::key_type key_type;

    key_type id() const { return _simId; }
    SimParticle const& simParticle() const { return *_simParticle; }

    size_t nHits() const { return _hitInfos.size(); }

    std::vector<StrawHitMCInfo>const& strawHitInfos() const { return _hitInfos; }

    StepPointMC const& firstStepPointMCinTracker() const;
    StepPointMC const& lastStepPointMCinTracker()  const;

    // Compiler generated code is Ok for:
    //   d'tor, copy c'tor assignment operator.
    // Once this class mature we will make the copy c'tor and assignment operator private.

  private:

    // This class should only ever be created by SimParticlesWithHits.
    // Therefore c'tor and non-const accessors are private.
    SimParticleInfo():_simId(-1){}

    SimParticleInfo( key_type simId,
                     SimParticle const& simParticle,
                     art::Event const& event);

    std::vector<StrawHitMCInfo>& strawHitInfos()  { return _hitInfos; }

    // ID of this particle in the SimParticleCollection.
    key_type _simId;

    // Pointer to the SimParticle
    SimParticle const* _simParticle;

    // The event in which this information is found.
    art::Event const* _event;

    // Vector of information about the StrawHits to which this track contributed.
    std::vector<StrawHitMCInfo>  _hitInfos;

    // First StepPointMC in tracker.  Lazy evaluated, therefore mutable.
    mutable StepPointMC const* _firstInTracker;

    // Last StepPointMC in tracker.  Lazy evaluated, therefore mutable.
    mutable StepPointMC const* _lastInTracker;

  };

} // namespace mu2e

#endif /* Mu2eUtilities_SimParticleInfo_hh */

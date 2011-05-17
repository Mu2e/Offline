#ifndef Mu2eUtilities_SimParticleInfo_HH
#define Mu2eUtilities_SimParticleInfo_HH
//
// Information about one SimParticle and all StrawHits that are
// associated with hit.  This is a building block of the
// the class SimParticlesWithHits.
//
// $Id: SimParticleInfo.hh,v 1.5 2011/05/17 15:36:01 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:01 $
//
// Original author Rob Kutschke.
//

// C++ includes.
#include <vector>

// Mu2e includes.
#include "Mu2eUtilities/inc/StrawHitMCInfo.hh"
#include "ToyDP/inc/SimParticleCollection.hh"

namespace art{
  class Event;
}


using namespace std;

namespace mu2e {
  
  class SimParticleInfo{

    // This class should only ever be created within SimParticlesWithHits.
    friend class SimParticlesWithHits;

  public:
    typedef SimParticleCollection::key_type key_type;

    key_type id() const { return _simId; }
    SimParticle const& simParticle() const { return *_simParticle; }

    size_t nHits() const { return _hitInfos.size(); }

    vector<StrawHitMCInfo>const& strawHitInfos() const { return _hitInfos; }
   
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

    vector<StrawHitMCInfo>& strawHitInfos()  { return _hitInfos; }

    // ID of this particle in the SimParticleCollection.
    key_type _simId;

    // Pointer to the SimParticle
    SimParticle const* _simParticle;

    // The event in which this information is found.
    art::Event const* _event;

    // Vector of information about the StrawHits to which this track contributed.
    vector<StrawHitMCInfo>  _hitInfos;

    // First StepPointMC in tracker.  Lazy evaluated, therefore mutable.
    mutable StepPointMC const* _firstInTracker; 

    // Last StepPointMC in tracker.  Lazy evaluated, therefore mutable.
    mutable StepPointMC const* _lastInTracker;

  };

} // namespace mu2e

#endif

#ifndef Mu2eUtilities_StrawHitMCInfo_hh
#define Mu2eUtilities_StrawHitMCInfo_hh
//
// Integrated access to all information about a StrawHit that was
// created by a SimParticle.   This class is a building block of
// the SimParticlesWithHits class. If a StrawHit contains contributions
// from two SimParticles, then there will usually be one two StrawHitMCInfo
// objects, one attached to each SimParticle.
//
// $Id: StrawHitMCInfo.hh,v 1.8 2011/06/07 21:41:08 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/07 21:41:08 $
//
// Original author Rob Kutschke.
//
// Notes:
// 1) The time returned by this class is defined as follows
//    Loop over all StepPointMCs that contributed to this hit.  Of these
//    select only those that were created by the trackId of this StrawHitMCInfo.
//    Of these, find the one with the earliest time.

// C++ includes.
#include <limits>

// Mu2e includes.
#include "MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

// Forward declarations.
namespace art{
  class Event;
}


namespace mu2e {

  // Forward declarations
  class StrawHit;
  class StrawHitMCTruth;
  class StepPointMC;

  class StrawHitMCInfo{
  public:

    typedef SimParticleCollection::key_type key_type;

    // No default c'tor by design.

    StrawHitMCInfo( art::Event           const& event,
                    key_type                    trackId,
                    size_t                      index,
                    StrawHit             const& strawHit,
                    StrawHitMCTruth      const& strawHitMCTruth,
                    PtrStepPointMCVector const& strawHitMCPtr,
                    int                         nSimParticles ):
      _index(index),
      _hit(&strawHit),
      _truth(&strawHitMCTruth),
      _mcPtr(&strawHitMCPtr),
      _nSimParticles(nSimParticles),
      _time(std::numeric_limits<double>::max()){
      fillStepPointMCs(event,trackId);
    }

    // Compiler generated code is Ok for:
    //  d'tor, copy c'tor assignment operator.

    size_t                      index()      const { return  _index;   }
    StrawHit             const& hit()        const { return *_hit;}
    StrawHitMCTruth      const& truth()      const { return *_truth;}
    PtrStepPointMCVector const& stepsByPrr() const { return *_mcPtr;}
    bool                        isShared()   const { return  _nSimParticles>1;}
    double                      time()       const { return  _time; }

    std::vector<StepPointMC const *> const& steps() const {
      return  _stepPointMCs;
    }

    // Return true if this hit occured before hit rhs.
    // The metric is the time of the earliest StepPointMC that
    // belongs to this hit and was made by the track to which it
    // is attached.
    bool operator<( StrawHitMCInfo const& rhs) const {
      return ( _time < rhs._time );
    }

  private:

    // Index into:
    // StrawHitCollection, StrawHitMCTruthCollection, PtrStepPointMCVectorCollection.
    size_t _index;

    // Non-owning pointers to the StrawHit, StrawHitMCTruth and PtrStepPointMCVector.
    StrawHit             const* _hit;
    StrawHitMCTruth      const* _truth;
    PtrStepPointMCVector const* _mcPtr;

    // Non-owning pointers to all StepPointMCs that contributed to this track.
    std::vector<StepPointMC const *> _stepPointMCs;

    // The number of distinct SimParticles that contribute to this hit.
    int _nSimParticles;

    // The time of this hit.  See note 1.
    double _time;

    // Fill _stepPointMCs and fill the variable time.
    void fillStepPointMCs(art::Event const& event, key_type trackId );
  };

} // namespace mu2e

#endif /* Mu2eUtilities_StrawHitMCInfo_hh */

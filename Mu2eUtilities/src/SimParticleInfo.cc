//
// Integrated access to all information about a StrawHit that was
// created by a SimParticle.   This class is a building block of
// the SimParticlesWithHits class. If a StrawHit contains contributions
// from two SimParticles, then there will usually be one two StrawHitMCInfo
// objects, one attached to each SimParticle.
//
// $Id: SimParticleInfo.cc,v 1.7 2011/10/28 18:47:07 greenc Exp $
// $Author: greenc $
// $Date: 2011/10/28 18:47:07 $
//
// See the notes in the header file for the meaning of the member datum _time.
//


// Framework includes
#include "art/Framework/Principal/Event.h"

// Mu2e includes
#include "Mu2eUtilities/inc/SimParticleInfo.hh"
#include "Mu2eUtilities/inc/StrawHitMCInfo.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// C++ includes.
#include <cfloat>

namespace mu2e {

    SimParticleInfo::SimParticleInfo(
            key_type simId,
            SimParticle const& simParticle,
            const art::Event& evt):
     _simId(simId),
     _simParticle(&simParticle),
     _event(&evt),
     _hitInfos(),
     _firstInTracker(0),
     _lastInTracker(0){}


    // calculate the first StepPointMC in track
    StepPointMC const& SimParticleInfo::firstStepPointMCinTracker() const {

        if ( _firstInTracker )  return *_firstInTracker;
        else if ( _hitInfos.size() == 0 )  {
          throw cet::exception("HITS")
            << "_hitInfos is empty at " << __func__ << "\n";
        }

        double currentFirstTime = DBL_MAX;
        StepPointMC const* currentStepPoint = 0;

        // Straws are already ordered by earliest time of StepMCPoint
        std::vector<StepPointMC const *> const& steps = _hitInfos.front().steps();

         // loop over all StepPointMC objects
         for ( std::vector<StepPointMC const *>::const_iterator stepsit = steps.begin(); stepsit!=steps.end(); ++stepsit){

             if ( (*stepsit)->trackId() == _simId ) {
                    double time = (*stepsit)->time();
                    if ( time < currentFirstTime ) {
                       currentFirstTime = time;
                       currentStepPoint = *stepsit;
                }
             }
         }
         _firstInTracker =  currentStepPoint;
         return *_firstInTracker;
    }

    // calculate the last StepPointMC in track
    StepPointMC const& SimParticleInfo::lastStepPointMCinTracker() const {

        if ( _lastInTracker )  return *_lastInTracker;
        else if ( _hitInfos.size() == 0 ) {
          throw cet::exception("HITS")
            << "_hitInfos is empty at " << __func__ << "\n";
        }

        double currentLastTime = -1;
        StepPointMC const* currentStepPoint = 0;

            // loop over all straws
            for ( std::vector<StrawHitMCInfo>::const_iterator strawsit = _hitInfos.begin(); strawsit!= _hitInfos.end(); ++strawsit ) {
                 std::vector<StepPointMC const *> const& steps = strawsit->steps();

                 // loop over all StepPointMC objects
                 for ( std::vector<StepPointMC const *>::const_iterator stepsit = steps.begin(); stepsit!=steps.end(); ++stepsit){

                     if ( (*stepsit)->trackId() == _simId ) {
                         double time = (*stepsit)->time();
                         if ( time > currentLastTime ) {
                           currentLastTime = time;
                           currentStepPoint = *stepsit;
                         }
                     }
                  }
              }
       _lastInTracker = currentStepPoint;
       return *_lastInTracker;
    }

} // end of namespace

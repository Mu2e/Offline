//
// Integrated access to all information about a StrawHit that was
// created by a SimParticle.   This class is a building block of
// the SimParticlesWithHits class. If a StrawHit contains contributions
// from two SimParticles, then there will usually be one two StrawHitMCInfo
// objects, one attached to each SimParticle.
//
// $Id: StrawHitMCInfo.cc,v 1.7 2011/10/28 18:47:07 greenc Exp $
// $Author: greenc $
// $Date: 2011/10/28 18:47:07 $
//
// Original author Rob Kutschke.
//
// See the notes in the header file for the meaning of the member datum _time.
//

// Framework includes
#include "art/Framework/Principal/Event.h"

// Mu2e includes
#include "Mu2eUtilities/inc/StrawHitMCInfo.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "MCDataProducts/inc/StepPointMC.hh"

namespace mu2e {

  void StrawHitMCInfo::fillStepPointMCs(art::Event const& event, key_type trackId ){
    for ( size_t i=0; i<_mcPtr->size(); ++i){

      StepPointMC const* step = _mcPtr->at(i).get();
      _stepPointMCs.push_back( step );

      if ( step->trackId() == trackId ){
        double t = step->time();
        _time = ( t < _time ) ? t : _time;
      }

    }

  } // StrawHitMCInfo::fillStepPointMCs

} // namespace mu2e

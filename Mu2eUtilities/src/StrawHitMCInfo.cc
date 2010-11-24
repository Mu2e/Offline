//
// Integrated access to all information about a StrawHit that was 
// created by a SimParticle.   This class is a building block of 
// the SimParticlesWithHits class. If a StrawHit contains contributions
// from two SimParticles, then there will usually be one two StrawHitMCInfo
// objects, one attached to each SimParticle. 
//
// $Id: StrawHitMCInfo.cc,v 1.1 2010/11/24 01:04:28 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/11/24 01:04:28 $
//
// Original author Rob Kutschke.
//
// See the notes in the header file for the meaning of the member datum _time.
// 

// Framework includes
#include "FWCore/Framework/interface/Event.h"

// Mu2e includes
#include "Mu2eUtilities/inc/StrawHitMCInfo.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/StrawHitMCPtr.hh"
#include "ToyDP/inc/StepPointMC.hh"

namespace mu2e {

  void StrawHitMCInfo::fillStepPointMCs(edm::Event const& event, key_type trackId ){
    for ( size_t i=0; i<_mcPtr->size(); ++i){

      StepPointMC const* step = 
        resolveDPIndex<StepPointMCCollection>(event, _mcPtr->at(i));

      _stepPointMCs.push_back( step );

      if ( step->trackId() == trackId ){
        double t = step->time();
        _time = ( t < _time ) ? t : _time;
      }

    }

  } // StrawHitMCInfo::fillStepPointMCs

} // namespace mu2e

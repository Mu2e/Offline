///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include "Mu2eUtilities/inc/McUtilsToolBase.hh"

namespace mu2e {
//-----------------------------------------------------------------------------
// ID of the SimParticle corresponding to the straw hit 'Index' in the 
// StrawHitCollection
//-----------------------------------------------------------------------------
  int McUtilsToolBase::strawHitSimId(const art::Event* Event, int Index) {
    return -1;
  }

//-----------------------------------------------------------------------------
//  double McUtilsToolBase::mcDoca(const art::Event* Event, int Index, const Straw* Straw) {
  double McUtilsToolBase::mcDoca(const art::Event* Event, const TrkStrawHit* StrawHit) { 
    return -99.;
  }


// //-----------------------------------------------------------------------------
//   int McUtilsToolBase::nGenHits(const art::Event*         Event         , 
// 				fhicl::ParameterSet*      TimeOffsets   ,
// 				const StrawHitCollection* Shcol         ) {
//     return -1;
//   }

  
// //-----------------------------------------------------------------------------
//   const StrawDigiMCCollection* McUtilsToolBase::getListOfMcStrawHits(const art::Event* Event,
// 								     const art::InputTag& Tag) {
//     return NULL;
//   }

//-----------------------------------------------------------------------------
  const SimParticle* McUtilsToolBase::getSimParticle(const art::Event* Event, int IHit) {
    return NULL;
  }
  
}

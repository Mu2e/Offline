///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include "Mu2eUtilities/inc/McUtilsToolBase.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"

namespace mu2e {

//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------
  double McUtilsToolBase::mcDoca(const art::Event* Event, const char* MCCollName, const Straw*  Straw) {
    return -99.;
  }


//-----------------------------------------------------------------------------
  int McUtilsToolBase::nGenHits(const art::Event*         Event         , 
				fhicl::ParameterSet*      TimeOffsets   ,
				const char*               MCDigiCollName, 
				const StrawHitCollection* Shcol         ) {
    return -1;
  }

  
//-----------------------------------------------------------------------------
  const PtrStepPointMCVectorCollection* McUtilsToolBase::getListOfMcStrawHits(const art::Event* Event,
									      const art::InputTag& Tag) {
    return NULL;
  }

//-----------------------------------------------------------------------------
  const SimParticle* McUtilsToolBase::getSimParticle(const PtrStepPointMCVectorCollection* List, int IHit) {
    return NULL;
  }
  
}

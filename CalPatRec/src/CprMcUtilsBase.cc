///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/CprMcUtilsBase.hh"

namespace mu2e {

//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------
  double CprMcUtilsBase::mcDoca(const art::Event* Event, const char* MCCollName, const Straw*  Straw) {
    return -99.;
  }


//-----------------------------------------------------------------------------
  int CprMcUtilsBase::nGenHits(const art::Event*         Event         , 
			       fhicl::ParameterSet*      TimeOffsets   ,
			       const char*               MCDigiCollName, 
			       const StrawHitCollection* Shcol         ) {
    return -1;
  }

}

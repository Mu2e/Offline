#ifndef __CalPatRec_CprMcUtilsBase_hh__
#define __CalPatRec_CprMcUtilsBase_hh__

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include "RecoDataProducts/inc/StrawHitCollection.hh"

namespace mu2e {

  class Straw;
  
  class CprMcUtilsBase {
  public:

    CprMcUtilsBase() noexcept = default ;
    virtual ~CprMcUtilsBase()  noexcept = default ;
//-----------------------------------------------------------------------------
// functiosn to be overloaded
//-----------------------------------------------------------------------------
    virtual double mcDoca(const art::Event* Event     , 
			  const char*       MCCollName, 
			  const Straw*      Straw     );

    virtual int    nGenHits(const art::Event*         Event         , 
			    fhicl::ParameterSet*      TimeOffsets   ,
			    const char*               MCDigiCollName, 
			    const StrawHitCollection* Shcol         );
  };
}

#endif

#ifndef __CalPatRec_McUtilsToolsBase_hh__
#define __CalPatRec_McUtilsToolsBase_hh__

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include "TrackerGeom/inc/Straw.hh"
#include "RecoDataProducts/inc/StrawHit.hh"

namespace mu2e {

  class Straw;
#ifndef MCDataProducts_StrawDigiMCCollection_hh
  class StrawDigiMCCollection;
#endif
  class SimParticle;
  
  class McUtilsToolBase {
  public:

    McUtilsToolBase() noexcept = default ;
    virtual ~McUtilsToolBase() noexcept = default ;
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

    virtual const StrawDigiMCCollection* getListOfMcStrawHits(const art::Event* Event,
							      const art::InputTag& Tag);

    virtual const SimParticle* getSimParticle(const StrawDigiMCCollection* List, int IHit);

    virtual int   getID      (const SimParticle* Sim) { return -1;  }
    virtual int   getPdgID   (const SimParticle* Sim) { return -1;  }
    virtual float getStartMom(const SimParticle* Sim) { return -1.; }
  };
}

#endif

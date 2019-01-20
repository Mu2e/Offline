///////////////////////////////////////////////////////////////////////////////
// no-empty implementation stores StrawDigiMCCollection
///////////////////////////////////////////////////////////////////////////////
#ifndef __CalPatRec_McUtilsToolsBase_hh__
#define __CalPatRec_McUtilsToolsBase_hh__

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

// #include "TrackerGeom/inc/Straw.hh"
// #include "RecoDataProducts/inc/StrawHit.hh"

namespace mu2e {

  //  class Straw;
  class TrkStrawHit;
// #ifndef MCDataProducts_StrawDigiMCCollection_hh
//   //  class StrawDigiMCCollection;
// #endif
  class SimParticle;
  
  class McUtilsToolBase {
  public:

    McUtilsToolBase()          noexcept = default ;
    virtual ~McUtilsToolBase() noexcept = default ;
//-----------------------------------------------------------------------------
// functions to be overloaded
//-----------------------------------------------------------------------------
    virtual int    strawHitSimId(const art::Event* Event, int Index);
    //    virtual double mcDoca       (const art::Event* Event, int Index, const Straw* Straw);
    virtual double mcDoca       (const art::Event* Event, const TrkStrawHit* StrawHit);

    // virtual int    nGenHits     (const art::Event*         Event         , 
    // 				 fhicl::ParameterSet*      TimeOffsets   ,
    // 				 const StrawHitCollection* Shcol         );

    // virtual const StrawDigiMCCollection* getListOfMcStrawHits(const art::Event*    Event,
    // 							      const art::InputTag& Tag  );

    virtual const SimParticle* getSimParticle(const art::Event* Event, int IHit);

    virtual int   getID      (const SimParticle* Sim) { return -1;  }
    virtual int   getPdgID   (const SimParticle* Sim) { return -1;  }
    virtual float getStartMom(const SimParticle* Sim) { return -1.; }
  };
}

#endif

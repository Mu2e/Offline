///////////////////////////////////////////////////////////////////////////////
// no-empty implementation stores StrawDigiMCCollection
///////////////////////////////////////////////////////////////////////////////
#ifndef __CalPatRec_McUtilsToolsBase_hh__
#define __CalPatRec_McUtilsToolsBase_hh__

#include "fhiclcpp/types/Atom.h"

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "Offline/DataProducts/inc/GenVector.hh"
// #include "TrackerGeom/inc/Straw.hh"
// #include "RecoDataProducts/inc/StrawHit.hh"

namespace mu2e {

  class TrkStrawHit;
  class SimParticle;

  class McUtilsToolBase {
  public:

    struct Config {
      fhicl::Atom<std::string>   tool_type         {fhicl::Name("tool_type"         ), fhicl::Comment("tool type: McUtils")     };
      fhicl::Atom<art::InputTag> comboHitCollTag   {fhicl::Name("comboHitCollTag"   ), fhicl::Comment("ComboHit collection tag")};
      fhicl::Atom<art::InputTag> strawDigiMCCollTag{fhicl::Name("strawDigiMCCollTag"), fhicl::Comment("StrawDigiMC Coll Tag")   };
    };

    McUtilsToolBase()          noexcept = default ;
    McUtilsToolBase(const Config& config) {};

    virtual ~McUtilsToolBase() noexcept = default ;
//-----------------------------------------------------------------------------
// functions to be overloaded
//-----------------------------------------------------------------------------
    virtual int    strawHitSimId(const art::Event* Event, int Index);
    //    virtual double mcDoca       (const art::Event* Event, int Index, const Straw* Straw);
    virtual double mcDoca       (const art::Event* Event, const TrkStrawHit* StrawHit);

    // virtual int    nGenHits     (const art::Event*         Event         ,
    //                                  fhicl::ParameterSet*      TimeOffsets   ,
    //                                  const StrawHitCollection* Shcol         );

    // virtual const StrawDigiMCCollection* getListOfMcStrawHits(const art::Event*    Event,
    //                                                               const art::InputTag& Tag  );

    virtual const SimParticle* getSimParticle(const art::Event* Event, int IHit);

    virtual const XYZVectorF* getMom(const art::Event* Event, int HitIndex) { return NULL; }
    virtual int   getID         (const SimParticle* Sim) { return -1;  }
    virtual int   getMotherID   (const SimParticle* Sim) { return -1;  }
    virtual int   getMotherPdgID(const SimParticle* Sim) { return -1;  }
    virtual int   getPdgID      (const SimParticle* Sim) { return -1;  }
    virtual float getStartMom   (const SimParticle* Sim) { return -1.; }

  };
}

#endif

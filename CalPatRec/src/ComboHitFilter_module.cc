//////////////////////////////////////////////////////////////////////////////
// creates a new collection of combohits containing only hits of one sim particle
// parameter defaults: CalPatRec/fcl/prolog.fcl
//////////////////////////////////////////////////////////////////////////////
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
// conditions
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
// data
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StereoHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"

#include "Offline/CalPatRec/inc/HlPrint.hh"

// diagnostics

#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "art/Utilities/make_tool.h"

#include <algorithm>
#include <cmath>

using namespace std;

namespace mu2e {

  using CLHEP::Hep3Vector;

  class ComboHitFilter: public art::EDProducer {
  public:
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>     chCollTag  {Name("chCollTag"  ), Comment("ComboHit collection Name"    ) };
      fhicl::Atom<art::InputTag>     chfCollTag {Name("chfCollTag" ), Comment("StrawHitFlag collection Name") };
      fhicl::Atom<art::InputTag>     sdmcCollTag{Name("sdmcCollTag"), Comment("StrawDigiMC collection Name" ) };
      fhicl::Atom<int>               simID      {Name("simID"      ), Comment("Selected sim particle ID"    ) };
      fhicl::Atom<int>               debugLevel {Name("debugLevel" ), Comment("debug level"                 ) };
    };
//-----------------------------------------------------------------------------
// talk-to parameters: input collections and algorithm parameters
//-----------------------------------------------------------------------------
    art::InputTag   _chfCollTag;
    art::InputTag   _chCollTag;
    art::InputTag   _sdmcCollTag;
    int             _simID;
    int             _debugLevel;
//-----------------------------------------------------------------------------
// cache event/geometry objects
//-----------------------------------------------------------------------------
    const ComboHitCollection*      _chColl  ;
    const StrawHitFlagCollection*  _chfColl ;  // combo hit flags
    const StrawDigiMCCollection*   _sdmcColl;

    int                            _nComboHits;
    int                            _nStrawHits;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit ComboHitFilter(const art::EDProducer::Table<ComboHitFilter::Config>& config);

  private:

    bool         findData     (const art::Event&  Evt);
//-----------------------------------------------------------------------------
// overloaded methods of the module class
//-----------------------------------------------------------------------------
    void beginJob()                    override;
    void beginRun(art::Run&   ARun   ) override;
    void produce (art::Event& AnEvent) override;
  };

//-----------------------------------------------------------------------------
  ComboHitFilter::ComboHitFilter(const art::EDProducer::Table<ComboHitFilter::Config>& config):
    art::EDProducer{config},
    _chfCollTag            (config().chfCollTag ()),
    _chCollTag             (config().chCollTag  ()),
    _sdmcCollTag           (config().sdmcCollTag()),
    _simID                 (config().simID      ()),
    _debugLevel            (config().debugLevel ())
  {
    consumes<ComboHitCollection>    (_chCollTag  );
    consumes<StrawHitFlagCollection>(_chfCollTag );
    consumes<StrawDigiMCCollection> (_sdmcCollTag);

    produces<ComboHitCollection>();
    produces<StrawHitFlagCollection>();
  }

  //-----------------------------------------------------------------------------
  void ComboHitFilter::beginJob() {
    if (_debugLevel > 0) {
    }
  }

//-----------------------------------------------------------------------------
// create a Z-ordered map of the tracker
//-----------------------------------------------------------------------------
  void ComboHitFilter::beginRun(art::Run& aRun) {
  }


//-----------------------------------------------------------------------------
  bool ComboHitFilter::findData(const art::Event& Evt) {

    auto chcH   = Evt.getValidHandle<ComboHitCollection>(_chCollTag);
    _chColl     = chcH.product();

    auto chfcH  = Evt.getValidHandle<StrawHitFlagCollection>(_chfCollTag);
    _chfColl    = chfcH.product();

    auto sdmccH = Evt.getValidHandle<StrawDigiMCCollection>(_sdmcCollTag);
    _sdmcColl   = sdmccH.product();

    return (_chColl != nullptr) and (_chfColl != nullptr) and (_sdmcColl != nullptr);
  }

//-----------------------------------------------------------------------------
  void ComboHitFilter::produce(art::Event& Event) {

    if (_debugLevel) printf(">>> ComboHitFilter::produce  event number: %10i\n",Event.event());
//-----------------------------------------------------------------------------
// process event
//-----------------------------------------------------------------------------
    if (! findData(Event)) {
      const char* message = "mu2e::ComboHitFilter_module::produce: data missing or incomplete";
      throw cet::exception("RECO")<< message << endl;
    }

    _nComboHits = _chColl->size();

    unique_ptr<ComboHitCollection>     new_chColl (new ComboHitCollection    ());
    unique_ptr<StrawHitFlagCollection> new_chfColl(new StrawHitFlagCollection());

    for (int i=0; i<_nComboHits; i++) {
      const ComboHit* ch  = &(*_chColl)[i];
      int ind = ch->indexArray().at(0);

      const mu2e::StrawDigiMC*  sdmc = &_sdmcColl->at(ind);
      const mu2e::StrawGasStep* step = sdmc->earlyStrawGasStep().get();
      const mu2e::SimParticle*  simp = step->simParticle().get();

      int sim_id = simp->id().asInt();
      if (_simID == sim_id) {
        new_chColl->push_back(*ch);
        new_chfColl->push_back((*_chfColl)[i]);
      }
    }
    if (_debugLevel > 0) {
      printf(" ComboHitFilter::%s : output collection : %lu hits\n",__func__,new_chColl->size());
    }
//-----------------------------------------------------------------------------
// finally, put the output flag collection into the event
//-----------------------------------------------------------------------------
    Event.put(std::move(new_chColl ));
    Event.put(std::move(new_chfColl));
//-----------------------------------------------------------------------------
// debug, if requested
//-----------------------------------------------------------------------------
    if (_debugLevel > 0) {
      HlPrint* hlp = HlPrint::Instance();
      hlp->SetEvent(&Event);
      hlp->printComboHitCollection("ComboHitFilter",
                                   "ComboHitFilter",
                                   _sdmcCollTag.encode().data());
    }
  }
}
//-----------------------------------------------------------------------------
// magic that makes this class a module.
//-----------------------------------------------------------------------------
DEFINE_ART_MODULE(mu2e::ComboHitFilter)
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------

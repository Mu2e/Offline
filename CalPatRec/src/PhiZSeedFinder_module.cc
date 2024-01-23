/////////////////////////////////////////////////////////////////////////////
// P.Murat
//
// flag combo hits as 'delta'  - hits of identified low energy electrons
//                and 'proton' - hits of identified protons/deuterons
//
// always writes out a ComboHitCollection with correct flags,
// to be used by downstream modules
//
// WriteFilteredComboHits = 0: write out all hits
//                        = 1: write out only hits not flagged as 'delta' or 'proton'
//                             (to be used in trigger)
//
// parameter defaults: CalPatRec/fcl/prolog.fcl
//////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"

#include "art/Utilities/make_tool.h"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

// conditions
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"

// data
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StereoHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"

// diagnostics
#include "Offline/CalPatRec/inc/ChannelID.hh"
#include "Offline/CalPatRec/inc/PhiZSeedFinder_types.hh"
#include "Offline/CalPatRec/inc/PhiZSeedFinderAlg.hh"

#include <algorithm>
#include <cmath>

using namespace std;

using CalPatRec::ChannelID;

namespace mu2e {

  using namespace PhiZSeedFinderTypes;

  class PhiZSeedFinder: public art::EDProducer {
  public:

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>   sschCollTag            {Name("sschCollTag"       )    , Comment("SS ComboHit collection tag" ) };
      fhicl::Atom<art::InputTag>   chCollTag              {Name("chCollTag"         )    , Comment("ComboHit collection tag"    ) };
      fhicl::Atom<art::InputTag>   tcCollTag              {Name("tcCollTag"         )    , Comment("time cluster collection tag") };
      fhicl::Atom<int>             debugLevel             {Name("debugLevel"        )    , Comment("debug level"                ) };
      fhicl::Atom<int>             diagLevel              {Name("diagLevel"         )    , Comment("diag level"                 ) };
      fhicl::Atom<int>             printErrors            {Name("printErrors"       )    , Comment("print errors"               ) };
      fhicl::Atom<int>             writeFilteredComboHits {Name("writeFilteredComboHits"), Comment("0: write all CH, 1: write filtered CH") };
      fhicl::Atom<int>             writeStrawHits         {Name("writeStrawHits"    )    , Comment("1: write all SH, new flags" ) };
      fhicl::Atom<int>             testOrder              {Name("testOrder"         )    , Comment("1: test order"              ) };
      fhicl::Atom<bool>            testHitMask            {Name("testHitMask"       )    , Comment("true: test hit mask"        ) };
      fhicl::Sequence<std::string> goodHitMask            {Name("goodHitMask"       )    , Comment("good hit mask"              ) };
      fhicl::Sequence<std::string> bkgHitMask             {Name("bkgHitMask"        )    , Comment("background hit mask"        ) };

      fhicl::Table<PhiZSeedFinderTypes::Config> diagPlugin      {Name("diagPlugin"      ), Comment("Diag plugin"           ) };
      fhicl::Table<PhiZSeedFinderAlg::Config>   finderParameters{Name("finderParameters"), Comment("finder alg parameters" ) };
    };

  protected:
//-----------------------------------------------------------------------------
// talk-to parameters: input collections and algorithm parameters
//-----------------------------------------------------------------------------
    art::InputTag   _sschCollTag;
    art::InputTag   _chCollTag;
    art::InputTag   _tcCollTag;                // time cluster coll tag

    int             _writeFilteredComboHits;   // write filtered combo hits
    // int             _writeStrawHitFlags;       // obsolete
    int             _writeStrawHits;           // write out filtered (?) straw hits

    int             _debugLevel;
    int             _diagLevel;
    int             _printErrors;
    int             _testOrder;

    StrawHitFlag    _bkgHitMask;

    std::unique_ptr<ModuleHistToolBase> _hmanager;
//-----------------------------------------------------------------------------
// cache event/geometry objects
//-----------------------------------------------------------------------------
    const ComboHitCollection*    _sschColl ;

    const Tracker*               _tracker;
    const DiskCalorimeter*       _calorimeter;

    PhiZSeedFinderTypes::Data_t  _data;              // all data used
    int                          _testOrderPrinted;

    PhiZSeedFinderAlg*           _finder;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit PhiZSeedFinder(const art::EDProducer::Table<Config>& config);

  private:

    bool         findData            (const art::Event&  Evt);
//-----------------------------------------------------------------------------
// overloaded methods of the module class
//-----------------------------------------------------------------------------
    void         beginJob() override;
    void         beginRun(art::Run& ARun) override;
    void         endJob  () override;
    void         produce (art::Event& E ) override;
  };

//-----------------------------------------------------------------------------
  PhiZSeedFinder::PhiZSeedFinder(const art::EDProducer::Table<Config>& config):
    art::EDProducer{config},
    _sschCollTag           (config().sschCollTag()       ),
    _chCollTag             (config().chCollTag()         ),
    _tcCollTag             (config().tcCollTag()         ),
    _debugLevel            (config().debugLevel()        ),
    _diagLevel             (config().diagLevel()         ),
    _printErrors           (config().printErrors()       ),
    _testOrder             (config().testOrder()         ),
    _bkgHitMask            (config().bkgHitMask()        )
  {

    consumes<TimeClusterCollection>(_tcCollTag);
    consumes<ComboHitCollection>   (_chCollTag);

    _finder = new PhiZSeedFinderAlg(config().finderParameters,&_data);

    _testOrderPrinted = 0;

    if (_diagLevel != 0) _hmanager = art::make_tool  <ModuleHistToolBase>(config().diagPlugin,"diagPlugin");
    else                 _hmanager = std::make_unique<ModuleHistToolBase>();

    _data.chCollTag      = _chCollTag;
    _data.tcCollTag      = _tcCollTag;
    _data._finder        = _finder;         // for diagnostics
  }

  //-----------------------------------------------------------------------------
  void PhiZSeedFinder::beginJob() {
    if (_diagLevel > 0) {
      art::ServiceHandle<art::TFileService> tfs;
      _hmanager->bookHistograms(tfs);
    }
  }

  //-----------------------------------------------------------------------------
  void PhiZSeedFinder::endJob() {
  }

//-----------------------------------------------------------------------------
// create a Z-ordered representation of the tracker
//-----------------------------------------------------------------------------
  void PhiZSeedFinder::beginRun(art::Run& aRun) {

    _data.InitGeometry();
//-----------------------------------------------------------------------------
// it is enough to print that once
//-----------------------------------------------------------------------------
    if (_testOrder && (_testOrderPrinted == 0)) {
      ChannelID::testOrderID  ();
      ChannelID::testdeOrderID();
      _testOrderPrinted = 1;
    }

    if (_diagLevel != 0) _hmanager->debug(&_data,1);
  }

//-----------------------------------------------------------------------------
  bool PhiZSeedFinder::findData(const art::Event& Evt) {

    auto tccH   = Evt.getValidHandle<mu2e::TimeClusterCollection>(_tcCollTag);
    _data.tccol = tccH.product();

    auto chcH = Evt.getValidHandle<mu2e::ComboHitCollection>(_chCollTag);
    _data.chcol = chcH.product();

    return (_data.tccol != nullptr) and (_data.chcol != nullptr);
  }

//-----------------------------------------------------------------------------
  void PhiZSeedFinder::produce(art::Event& Event) {
    if (_debugLevel) printf("* >>> PhiZSeedFinder::produce  event number: %10i\n",Event.event());
//-----------------------------------------------------------------------------
// clear memory in the beginning of event processing and cache event pointer
//-----------------------------------------------------------------------------
    _data.InitEvent(&Event,_debugLevel);
//-----------------------------------------------------------------------------
// process event
//-----------------------------------------------------------------------------
    if (! findData(Event)) {
      const char* message = "mu2e::PhiZSeedFinder_module::produce: data missing or incomplete";
      throw cet::exception("RECO")<< message << endl;
    }

    _data._nTimeClusters = _data.tccol->size();
//-----------------------------------------------------------------------------
// run delta finder, it also finds proton time clusters
//-----------------------------------------------------------------------------
    for (int i=0; i<_data._nTimeClusters; i++) {
      const TimeCluster* tc = &_data.tccol->at(i);
      _finder->run(tc);
    }
  }
}
//-----------------------------------------------------------------------------
// macro that makes this class a module.
//-----------------------------------------------------------------------------
DEFINE_ART_MODULE(mu2e::PhiZSeedFinder)
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------

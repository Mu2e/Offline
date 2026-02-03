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
#include "Offline/Mu2eUtilities/inc/StopWatch.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
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
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"

// diagnostics

#include "Offline/CalPatRec/inc/DeltaFinder_types.hh"
#include "Offline/CalPatRec/inc/DeltaFinderAlg.hh"

#include <algorithm>
#include <cmath>

using namespace std;

namespace mu2e {

  using namespace DeltaFinderTypes;

  class DeltaFinder: public art::EDProducer {
  public:

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>   sschCollTag            {Name("sschCollTag"       )    , Comment("SS ComboHit collection name") };
      fhicl::Atom<art::InputTag>   chCollTag              {Name("chCollTag"         )    , Comment("ComboHit collection Name"   ) };
      fhicl::Atom<art::InputTag>   ewmTag                 {Name("ewmTag"            )    , Comment("EventWindowMarker tag") };
      fhicl::Atom<art::InputTag>   sdmcCollTag            {Name("sdmcCollTag"       )    , Comment("StrawDigiMC collection Name") };
      fhicl::Atom<int>             debugLevel             {Name("debugLevel"        )    , Comment("debug level"                ) };
      fhicl::Atom<int>             diagLevel              {Name("diagLevel"         )    , Comment("diag level"                 ) };
      fhicl::Atom<int>             printErrors            {Name("printErrors"       )    , Comment("print errors"               ) };
      fhicl::Atom<int>             writeFilteredComboHits {Name("writeFilteredComboHits"), Comment("0: write all CH, 1: write filtered CH") };
      fhicl::Atom<int>             writeStrawHits         {Name("writeStrawHits"    )    , Comment("1: write all SH, new flags" ) };
      fhicl::Atom<int>             testOrder              {Name("testOrder"         )    , Comment("1: test order"              ) };
      fhicl::Atom<bool>            testHitMask            {Name("testHitMask"       )    , Comment("true: test hit mask"        ) };
      fhicl::Atom<int>             doTiming               {Name("doTiming"          )    , Comment("perform timing analysis"    ), 0};
      fhicl::Sequence<std::string> goodHitMask            {Name("goodHitMask"       )    , Comment("good hit mask"              ) };
      fhicl::Sequence<std::string> bkgHitMask             {Name("bkgHitMask"        )    , Comment("background hit mask"        ) };

      fhicl::Table<DeltaFinderTypes::Config> diagPlugin      {Name("diagPlugin"      ), Comment("Diag plugin"           ) };
      fhicl::Table<DeltaFinderAlg::Config>   finderParameters{Name("finderParameters"), Comment("finder alg parameters" ) };
    };

  protected:
//-----------------------------------------------------------------------------
// talk-to parameters: input collections and algorithm parameters
//-----------------------------------------------------------------------------
    art::InputTag   _sschCollTag;
    art::InputTag   _chCollTag;
    art::InputTag   _ewmTag;
    art::InputTag   _sdmcCollTag;

    int             _writeFilteredComboHits;   // write filtered combo hits
    // int             _writeStrawHitFlags;       // obsolete
    int             _writeStrawHits;           // write out filtered (?) straw hits

    int             _debugLevel;
    int             _diagLevel;
    int             _printErrors;
    int             _testOrder;
    int             _doTiming;

    StrawHitFlag    _bkgHitMask;

    std::unique_ptr<ModuleHistToolBase> _hmanager;
    std::shared_ptr<StopWatch>          _watch;
//-----------------------------------------------------------------------------
// cache event/geometry objects
//-----------------------------------------------------------------------------
    const ComboHitCollection*    _sschColl ;

    const Tracker*               _tracker;
    const DiskCalorimeter*       _calorimeter;

    DeltaFinderTypes::Data_t     _data;              // all data used
    int                          _testOrderPrinted;

    DeltaFinderAlg*              _finder;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit DeltaFinder(const art::EDProducer::Table<Config>& config);

  private:

    bool         findData            (const art::Event&  Evt);
    void         initTimeCluster     (DeltaCandidate* Delta, TimeCluster* Tc);
//-----------------------------------------------------------------------------
// overloaded methods of the module class
//-----------------------------------------------------------------------------
    void         beginJob() override;
    void         beginRun(art::Run& ARun) override;
    void         endJob  () override;
    void         produce (art::Event& E ) override;
  };

//-----------------------------------------------------------------------------
  DeltaFinder::DeltaFinder(const art::EDProducer::Table<Config>& config):
    art::EDProducer{config},
    _sschCollTag           (config().sschCollTag()       ),
    _chCollTag             (config().chCollTag()         ),
    _ewmTag                (config().ewmTag()            ),
    _sdmcCollTag           (config().sdmcCollTag()       ),
    _writeFilteredComboHits(config().writeFilteredComboHits()    ),
    // _writeStrawHitFlags    (config().writeStrawHitFlags()),
    _writeStrawHits        (config().writeStrawHits()    ),
    _debugLevel            (config().debugLevel()        ),
    _diagLevel             (config().diagLevel()         ),
    _printErrors           (config().printErrors()       ),
    _testOrder             (config().testOrder()         ),
    _doTiming              (config().doTiming()          ),
    _bkgHitMask            (config().bkgHitMask()        )
  {

    consumesMany<ComboHitCollection>(); // ??? Necessary because fillStrawHitIndices calls getManyByType.

    produces<IntensityInfoTimeCluster>();

    produces<ComboHitCollection>("");
    if (_writeStrawHits         == 1) produces<ComboHitCollection>("StrawHits");

                                        // this is a list of delta-electron candidates (or proton ones ?)
    produces<TimeClusterCollection>();

    if(_doTiming) {
      _watch = std::make_shared<StopWatch>();
      _data.doTiming = _doTiming; // add the watch to the data so the finder can use it
      _data.watch = _watch;
    }

    _finder = new DeltaFinderAlg(config().finderParameters,&_data);
    _data.timeBin = config().finderParameters().timeBin();

    _testOrderPrinted = 0;

    if (_diagLevel != 0) _hmanager = art::make_tool  <ModuleHistToolBase>(config().diagPlugin,"diagPlugin");
    else                 _hmanager = std::make_unique<ModuleHistToolBase>();

    _data.chCollTag      = _chCollTag;
    _data.sdmcCollTag    = _sdmcCollTag;
    _data._finder        = _finder;         // for diagnostics
  }

  //-----------------------------------------------------------------------------
  void DeltaFinder::beginJob() {
    if (_diagLevel > 0) {
      art::ServiceHandle<art::TFileService> tfs;
      _hmanager->bookHistograms(tfs);
    }
  }

  //-----------------------------------------------------------------------------
  void DeltaFinder::endJob() {
    if(_doTiming) {
      std::cout << "[DeltaFinder::" << __func__ << "::" <<
        moduleDescription().moduleLabel() << "] Timing:\n" << *_watch;
    }
  }

//-----------------------------------------------------------------------------
// create a Z-ordered representation of the tracker
//-----------------------------------------------------------------------------
  void DeltaFinder::beginRun(art::Run& aRun) {

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
  bool DeltaFinder::findData(const art::Event& Evt) {
    if(_doTiming > 0) _watch->SetTime(__func__);

    auto chcH   = Evt.getValidHandle<mu2e::ComboHitCollection>(_chCollTag);
    _data.chcol = chcH.product();

    auto sschcH = Evt.getValidHandle<mu2e::ComboHitCollection>(_sschCollTag);
    _sschColl   = sschcH.product();

    art::Handle<mu2e::EventWindowMarker> ewmHandle;
    if(Evt.getByLabel(_ewmTag, ewmHandle)) _data.ewm = ewmHandle.product();

    if(_doTiming > 0) _watch->StopTime(__func__);
    return (_data.chcol != nullptr) and (_sschColl != nullptr);
  }

//-----------------------------------------------------------------------------
// define the time cluster parameters for found DeltaCandidates
//-----------------------------------------------------------------------------
  void DeltaFinder::initTimeCluster(DeltaCandidate* Dc, TimeCluster* Tc) {
  }

//-----------------------------------------------------------------------------
  void DeltaFinder::produce(art::Event& Event) {
    if(_doTiming > 0) _watch->SetTime(__func__);
    if (_debugLevel) printf("* >>> DeltaFinder::produce  event number: %10i\n",Event.event());
//-----------------------------------------------------------------------------
// clear memory in the beginning of event processing and cache event pointer
//-----------------------------------------------------------------------------
    if (! findData(Event)) { // need to retrieve data before clearing data, to know the event type
      const char* message = "mu2e::DeltaFinder_module::produce: data missing or incomplete";
      throw cet::exception("RECO")<< message << endl;
    }
    if(_doTiming > 0) _watch->SetTime("init-event-data");
    _data.InitEvent(&Event,_debugLevel);
    if(_doTiming > 0) _watch->StopTime("init-event-data");
//-----------------------------------------------------------------------------
// process event
//-----------------------------------------------------------------------------

    _data._nComboHits = _data.chcol->size();
    _data._nStrawHits = _sschColl->size();
//-----------------------------------------------------------------------------
// run delta finder, it also finds proton time clusters
//-----------------------------------------------------------------------------
    _finder->run();
//-----------------------------------------------------------------------------
// done with the pattern recognition part
// create a temporary straw hit flag collection
// need to drop previously set 'energy' flag
// flag collections are going away, keep this one as a temp storage
//-----------------------------------------------------------------------------
    if(_doTiming > 0) _watch->SetTime("finalize-results");
    vector<StrawHitFlag> up_chfcol(_data._nComboHits);

    for (int i=0; i<_data._nComboHits; i++) {
//-----------------------------------------------------------------------------
// initialize output flags to the flags of the input combo hits, in parallel
// count the number of hits in the potentially to be written out straw hit collection
// the number of the output combo hits is the same as the number of input ones,
// no filtering here
//-----------------------------------------------------------------------------
      const ComboHit* ch = &(*_data.chcol)[i];
      StrawHitFlag* flag = &up_chfcol[i];
      flag->merge(ch->flag());
//-----------------------------------------------------------------------------
// always flag delta hits here, don't need previously set delta bits
// if flag protons, flag all hits as good
//-----------------------------------------------------------------------------
      flag->clear(StrawHitFlag::bkg);
      if (_finder->flagProtonHits()) flag->merge(StrawHitFlag::energysel);
    }

    const ComboHit* ch0(0);
    if (_data._nComboHits > 0) ch0 = &_data.chcol->at(0);
//-----------------------------------------------------------------------------
// loop over delta candidates and protons and update the flags
// a) set delta ('bkg') flags
//-----------------------------------------------------------------------------
    int ndeltas = _data.nDeltaCandidates();
    for (int i=0; i<ndeltas; i++) {
      DeltaCandidate* dc = _data.deltaCandidate(i);
//-----------------------------------------------------------------------------
// skip merged in delta candidates
// also require a delta candidate to have at least 5 hits (in the mask, set by
// the delta finder)
// do not consider proton stub candidates (those with <EDep> > 0.004)
//-----------------------------------------------------------------------------
      if (dc->Active() == 0)                                          continue;
      if (dc->Mask()   != 0)                                          continue;
      for (int is=dc->fFirstStation; is<=dc->fLastStation; is++) {
        DeltaSeed* ds = dc->Seed(is);
        if (ds != nullptr) {
//-----------------------------------------------------------------------------
// loop over the hits and flag each of them as delta
//-----------------------------------------------------------------------------
          for (int face=0; face<kNFaces; face++) {
            const HitData_t* hd = ds->HitData(face);
            if (hd == nullptr)                                        continue;
            int loc = hd->fHit-ch0;
            StrawHitFlag* flag = &up_chfcol[loc];
            flag->merge(StrawHitFlag::bkg);           // set delta-electron bit
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// set proton flags, 'energy' now means 'proton'
// tcColl represents the proton time clusters
//-----------------------------------------------------------------------------
    unique_ptr<TimeClusterCollection>  tcColl(new TimeClusterCollection);

    int np15(0);
    int npc = _data.nProtonCandidates();
    for (int i=0; i<npc; i++) {
      ProtonCandidate* pc = _data.protonCandidate(i);
      if (pc->nHitsTot() >= 15) np15++;
      if (pc->nStationsWithHits() == 1) {
//-----------------------------------------------------------------------------
// for proton candidates with just one station require eDep > 4 KeV
// this adds a little inefficiency for proton reco,
// but reduces overefficiency of flagging the CE hits
//-----------------------------------------------------------------------------
        if (pc->eDep() < 0.004) continue;
      }
      for (int is=pc->fFirstStation; is<=pc->fLastStation; is++) {
        for (int face=0; face<kNFaces; face++) {
          int nh = pc->nHits(is,face);
          for (int ih=0; ih<nh; ih++) {
            const HitData_t* hd = pc->hitData(is,face,ih);
            int loc = hd->fHit-ch0;
            StrawHitFlag* flag = &up_chfcol[loc];
            flag->clear(StrawHitFlag::energysel);
          }
        }
      }
//-----------------------------------------------------------------------------
// (later) make a time cluster out of each ProtonCandidate
//-----------------------------------------------------------------------------
      // TimeCluster new_tc;
      // initTimeCluster(dc,&new_tc);
      // tcColl->push_back(new_tc);
    }
//-----------------------------------------------------------------------------
// form the output - flag combo hits -
// if flagged combo hits are written out, likely don't need writing out the flags
// 2023-06-02: from now on, always write the output collection with updated flags
// if _writeFilteredComboHits=true (trigger), write out only "good" hits
//
// write out filtered out collection of ComboHits with right flags, use deep copy
// the output ComboHit collection doesn't need a parallel flag collection
//-----------------------------------------------------------------------------
    auto outputChColl = std::make_unique<ComboHitCollection>();
    outputChColl->reserve(_data._nComboHits);
    _data.outputChColl = outputChColl.get();

    outputChColl->setParent(_data.chcol->parent());
    for (int i=0; i<_data._nComboHits; i++) {
      StrawHitFlag* flag = & up_chfcol[i];
      if (_writeFilteredComboHits and flag->hasAnyProperty(_bkgHitMask)) continue;
//-----------------------------------------------------------------------------
// for the moment, assume bkgHitMask to be empty, so write out all hits
// normally, don't need delta and proton hits
//-----------------------------------------------------------------------------
      const ComboHit* ch = &(*_data.chcol)[i];
      outputChColl->push_back(*ch);
      //      outputChColl->back()._flag.merge(*flag);
      outputChColl->back()._flag = *flag;
    }
//-----------------------------------------------------------------------------
// single-straw hits are also combo hits and they have flags
// flags need to be redefined
//-----------------------------------------------------------------------------
    unique_ptr<ComboHitCollection> outputSschColl;

    if (_writeStrawHits == 1) {
                                                          // first, copy over the original hits

      outputSschColl = std::make_unique<ComboHitCollection>(*_sschColl);

      outputSschColl->setParent(_sschColl->parent());

      for (int i=0; i<_data._nStrawHits; i++) {
        StrawHitFlag* flag = (StrawHitFlag*) &outputSschColl->at(i).flag(); // original
        flag->clear(StrawHitFlag::bkg       );             // clear delta selection
        flag->merge(StrawHitFlag::energysel );             // and assume all hits are not from protons
      }
//-----------------------------------------------------------------------------
// after that, loop over [flagged] combo hits and update 'delta'(bkg) and 'proton'(energysel)
// flags for flagged hits
//-----------------------------------------------------------------------------
      for (int ich=0; ich<_data._nComboHits; ich++) {
        const ComboHit* ch   = &(*_data.chcol )[ich];
        StrawHitFlag    flag = up_chfcol[ich];
        for (auto ish : ch->indexArray()) {
          StrawHitFlag* shflag = (StrawHitFlag*) &(*outputSschColl)[ish].flag();
          if (not flag.hasAnyProperty(StrawHitFlag::energysel)) shflag->clear(StrawHitFlag::energysel);
          if (    flag.hasAnyProperty(StrawHitFlag::bkg      )) shflag->merge(StrawHitFlag::bkg);
        }
      }

      //      Event.put(std::move(outputSschColl),"StrawHits");
    }
    if(_doTiming > 0) _watch->StopTime("finalize-results");
//-----------------------------------------------------------------------------
// 'ppii' - proton counting-based proxy to the stopped muon rate
//-----------------------------------------------------------------------------
    std::unique_ptr<IntensityInfoTimeCluster> ppii(new IntensityInfoTimeCluster(np15));
//-----------------------------------------------------------------------------
// in the end of event processing fill diagnostic histograms
//-----------------------------------------------------------------------------
    if (_diagLevel  > 0) _hmanager->fillHistograms(&_data);
    if (_debugLevel > 0) _hmanager->debug(&_data,2);
//-----------------------------------------------------------------------------
// moving in the end, after diagnostics plugin routines have been called - move
// invalidates the original pointer...
//-----------------------------------------------------------------------------
    if(_doTiming > 0) _watch->SetTime("output");
    Event.put(std::move(outputChColl));
    if (_writeStrawHits == 1) Event.put(std::move(outputSschColl),"StrawHits");

    Event.put(std::move(tcColl));
    Event.put(std::move(ppii));
    if(_doTiming > 0) _watch->StopTime("output");
    if(_doTiming > 0) _watch->StopTime(__func__);
  }
}
//-----------------------------------------------------------------------------
// macro that makes this class a module.
//-----------------------------------------------------------------------------
DEFINE_ART_MODULE(mu2e::DeltaFinder)
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------

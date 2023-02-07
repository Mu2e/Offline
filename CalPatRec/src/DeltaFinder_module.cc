/////////////////////////////////////////////////////////////////////////////
// framework
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
      fhicl::Atom<art::InputTag>   shCollTag         {Name("shCollTag"         ), Comment("SComboHit collection Name"   ) };
      fhicl::Atom<art::InputTag>   chCollTag         {Name("chCollTag"         ), Comment("ComboHit collection Name"    ) };
      fhicl::Atom<art::InputTag>   sdmcCollTag       {Name("sdmcCollTag"       ), Comment("StrawDigiMC collection Name" ) };
      fhicl::Atom<int>             debugLevel        {Name("debugLevel"        ), Comment("debug level"                 ) };
      fhicl::Atom<int>             diagLevel         {Name("diagLevel"         ), Comment("diag level"                  ) };
      fhicl::Atom<int>             printErrors       {Name("printErrors"       ), Comment("print errors"                ) };
      fhicl::Atom<int>             writeComboHits    {Name("writeComboHits"    ), Comment("if 1, write combohit coll"   ) };
      fhicl::Atom<int>             writeStrawHitFlags{Name("writeStrawHitFlags"), Comment("if 1, write SH flag coll"    ) };
      fhicl::Atom<int>             testOrder         {Name("testOrder"         ), Comment("if 1, test order"            ) };
      fhicl::Atom<bool>            testHitMask       {Name("testHitMask"       ), Comment("if true, test hit mask"      ) };
      fhicl::Sequence<std::string> goodHitMask       {Name("goodHitMask"       ), Comment("good hit mask"               ) };
      fhicl::Sequence<std::string> bkgHitMask        {Name("bkgHitMask"        ), Comment("background hit mask"         ) };

      fhicl::Table<DeltaFinderTypes::Config> diagPlugin      {Name("diagPlugin"      ), Comment("Diag plugin"           ) };
      fhicl::Table<DeltaFinderAlg::Config>   finderParameters{Name("finderParameters"), Comment("finder alg parameters" ) };
    };

  protected:
//-----------------------------------------------------------------------------
// talk-to parameters: input collections and algorithm parameters
//-----------------------------------------------------------------------------
    art::InputTag   _shCollTag;
    art::InputTag   _chCollTag;
    art::InputTag   _sdmcCollTag;

    int             _writeComboHits;       // write (filtered ?) combo hits
    int             _writeStrawHitFlags;

    int             _debugLevel;
    int             _diagLevel;
    int             _printErrors;
    int             _testOrder;

    StrawHitFlag    _bkgHitMask;

    std::unique_ptr<ModuleHistToolBase> _hmanager;
//-----------------------------------------------------------------------------
// cache event/geometry objects
//-----------------------------------------------------------------------------
    const StrawHitCollection*    _shColl ;

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
    _shCollTag             (config().shCollTag()         ),
    _chCollTag             (config().chCollTag()         ),
    _sdmcCollTag           (config().sdmcCollTag()       ),
    _writeComboHits        (config().writeComboHits()    ),
    _writeStrawHitFlags    (config().writeStrawHitFlags()),
    _debugLevel            (config().debugLevel()        ),
    _diagLevel             (config().diagLevel()         ),
    _printErrors           (config().printErrors()       ),
    _testOrder             (config().testOrder()         ),
    _bkgHitMask            (config().bkgHitMask()        )
  {

    consumesMany<ComboHitCollection>(); // Necessary because fillStrawHitIndices calls getManyByType.

    produces<StrawHitFlagCollection>("ComboHits");
    if (_writeStrawHitFlags == 1) produces<StrawHitFlagCollection>("StrawHits");
    if (_writeComboHits     == 1) produces<ComboHitCollection>    ("");

                                        // this is a list of delta-electron candidates
    produces<TimeClusterCollection>();

    _finder = new DeltaFinderAlg(config().finderParameters,&_data);

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
      _data.testOrderID  ();
      _data.testdeOrderID();
      _testOrderPrinted = 1;
    }

    if (_diagLevel != 0) _hmanager->debug(&_data,1);
  }

//-----------------------------------------------------------------------------
  bool DeltaFinder::findData(const art::Event& Evt) {

    auto chcH   = Evt.getValidHandle<mu2e::ComboHitCollection>(_chCollTag);
    _data.chcol = chcH.product();

    auto shcH = Evt.getValidHandle<mu2e::StrawHitCollection>(_shCollTag);
    _shColl   = shcH.product();

    return (_data.chcol != nullptr) and (_shColl != nullptr);
  }

//-----------------------------------------------------------------------------
// define the time cluster parameters starting from a DeltaCandidate
//-----------------------------------------------------------------------------
  void DeltaFinder::initTimeCluster(DeltaCandidate* Dc, TimeCluster* Tc) {
  }

//-----------------------------------------------------------------------------
  void DeltaFinder::produce(art::Event& Event) {
    if (_debugLevel) printf(">>> DeltaFinder::produce  event number: %10i\n",Event.event());
//-----------------------------------------------------------------------------
// clear memory in the beginning of event processing and cache event pointer
//-----------------------------------------------------------------------------
    _data.InitEvent(&Event,_debugLevel);
//-----------------------------------------------------------------------------
// process event
//-----------------------------------------------------------------------------
    if (! findData(Event)) {
      const char* message = "mu2e::DeltaFinder_module::produce: data missing or incomplete";
      throw cet::exception("RECO")<< message << endl;
    }

    _data._nComboHits = _data.chcol->size();
    _data._nStrawHits = _shColl->size();

    _finder->run();
//-----------------------------------------------------------------------------
// form output - flag combo hits -
// if flagged combo hits are written out, likely don't need writing out the flags
//-----------------------------------------------------------------------------
    unique_ptr<StrawHitFlagCollection> up_chfcol(new StrawHitFlagCollection(_data._nComboHits));
    _data.outputChfColl = up_chfcol.get();

    for (int i=0; i<_data._nComboHits; i++) {
      const ComboHit* ch = &(*_data.chcol)[i];
      (*_data.outputChfColl)[i].merge(ch->flag());
    }

    const ComboHit* ch0(0);
    if (_data._nComboHits > 0) ch0 = &_data.chcol->at(0);

    StrawHitFlag deltamask(StrawHitFlag::bkg);

    unique_ptr<TimeClusterCollection>  tcColl(new TimeClusterCollection);

    int ndeltas = _data.nDeltaCandidates();

    for (int i=0; i<ndeltas; i++) {
      DeltaCandidate* dc = _data.deltaCandidate(i);
//-----------------------------------------------------------------------------
// skip merged in delta candidates
// also require a delta candidate to have at least 5 hits
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
            _data.outputChfColl->at(loc).merge(deltamask);
          }
        }
      }
//-----------------------------------------------------------------------------
// make a time cluster out of each active DeltaCandidate
//-----------------------------------------------------------------------------
      TimeCluster new_tc;
      initTimeCluster(dc,&new_tc);
      tcColl->push_back(new_tc);
    }

    Event.put(std::move(tcColl));
//-----------------------------------------------------------------------------
// in the end of event processing fill diagnostic histograms
//-----------------------------------------------------------------------------
    if (_diagLevel  > 0) _hmanager->fillHistograms(&_data);
    if (_debugLevel > 0) _hmanager->debug(&_data,2);

    if (_writeComboHits) {
//-----------------------------------------------------------------------------
// write out collection of ComboHits with right flags, use deep copy
//-----------------------------------------------------------------------------
      auto outputChColl = std::make_unique<ComboHitCollection>();
      outputChColl->reserve(_data._nComboHits);

      outputChColl->setParent(_data.chcol->parent());
      for (int i=0; i<_data._nComboHits; i++) {
        StrawHitFlag const* flag = &(*_data.outputChfColl)[i];
        if (flag->hasAnyProperty(_bkgHitMask))                        continue;
//-----------------------------------------------------------------------------
// for the moment, assume bkgHitMask to be empty, so write out all hits
//-----------------------------------------------------------------------------
        const ComboHit* ch = &(*_data.chcol)[i];
        outputChColl->push_back(*ch);
        outputChColl->back()._flag.merge(*flag);
      }
      Event.put(std::move(outputChColl));
    }
//-----------------------------------------------------------------------------
// create the collection of StrawHitFlag for the StrawHitCollection
//-----------------------------------------------------------------------------
    if (_writeStrawHitFlags == 1) {
                                        // first, copy over the original flags

      std::unique_ptr<StrawHitFlagCollection> shfcol(new StrawHitFlagCollection(_data._nStrawHits));

      for(int ich=0; ich<_data._nComboHits; ich++) {
        const ComboHit* ch   = &(*_data.chcol )[ich];
        StrawHitFlag    flag =  (*_data.outputChfColl)[ich];
        flag.merge(ch->flag());
        for (auto ish : ch->indexArray()) {
          (*shfcol)[ish] = flag;
        }
      }

      Event.put(std::move(shfcol),"StrawHits");
    }
//-----------------------------------------------------------------------------
// moving in the end, after diagnostics plugin routines have been called - move
// invalidates the original pointer...
//-----------------------------------------------------------------------------
    Event.put(std::move(up_chfcol),"ComboHits");
  }

}
//-----------------------------------------------------------------------------
// magic that makes this class a module.
//-----------------------------------------------------------------------------
DEFINE_ART_MODULE(mu2e::DeltaFinder)
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------

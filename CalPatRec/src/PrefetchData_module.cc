// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"
// conditions
#include "Offline/TrackerGeom/inc/Tracker.hh"
// root
#include "TMath.h"
#include "TH1F.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2F.h"
#include "TVector2.h"
// data
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
// diagnostics
#include <algorithm>
#include <cmath>
#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"
using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class PrefetchData : public art::EDProducer {

  protected:

  public:
    explicit PrefetchData(fhicl::ParameterSet const&);
    virtual ~PrefetchData();
    virtual void beginJob();
    virtual void produce( art::Event& e);

    void   fake_access(const CaloDigi&    Digi);
    void   fake_access(const CaloHit&     Hit);

    void   fake_access(const StrawHit& Hit, /*const StrawHitFlag& Flag,*/ const StrawHitPosition& Pos);
    void   fake_access(const ComboHit&    Hit);
    // void   fake_access(const StereoHit&   Hit);
    void   fake_access(const StrawDigi&   Digi);
    void   fake_access(const StepPointMC* Step);

  private:

    bool findData(const art::Event& e);
                                        // control flags
    int           _debugLevel;
    bool          _mcDiag;
    int           _fetchCaloDigis;
    int           _fetchCaloHits;
    int           _fetchStrawHits;
    int           _fetchComboHits;
    int           _fetchStrawHitFlags;
    int           _fetchStrawHitPositions;
    // int           _fetchStereoHits;
    int           _fetchStrawDigis;
                                        // data tags
    art::InputTag _cdTag;
    art::InputTag _cchTag;

    art::InputTag _shTag;
    art::InputTag _chTag;
    art::InputTag _sthTag;
    art::InputTag _shfTag;
    art::InputTag _shpTag;
    art::InputTag _sdTag;
                                        // cache of event objects
    const CaloDigiCollection*                   _cdcol;
    const CaloHitCollection*                    _cchcol;

    const StrawHitCollection*                   _shcol;
    const ComboHitCollection*                   _chcol;
    // const StereoHitCollection*                  _sthcol;
    const StrawHitFlagCollection*               _shfcol;
    const StrawHitPositionCollection*           _shpcol;
    const StrawDigiCollection*                  _sdcol;
    int                                         _eventNum;
  };

  //-----------------------------------------------------------------------------
  PrefetchData::PrefetchData(fhicl::ParameterSet const& pset):
    art::EDProducer(pset),
    _debugLevel     (pset.get<int>          ("debugLevel"              )),
    _mcDiag         (pset.get<bool>         ("mcDiag"                  )),

    _fetchCaloDigis (pset.get<int>         ("fetchCaloDigis" )),
    _fetchCaloHits  (pset.get<int>         ("fetchCaloHits" )),
    _fetchStrawHits (pset.get<int>         ("fetchStrawHits" )),
    _fetchComboHits (pset.get<int>         ("fetchComboHits" )),
    _fetchStrawHitFlags (pset.get<int>     ("fetchStrawHitFlags" )),
    _fetchStrawHitPositions (pset.get<int> ("fetchStrawHitPositions" )),
    // _fetchStereoHits (pset.get<int>        ("fetchStereoHits" )),
    _fetchStrawDigis(pset.get<int>         ("fetchStrawDigis")),

    _cdTag     (pset.get<string>       ("caloDigiCollectionTag"        )),
    _cchTag     (pset.get<string>       ("caloHitCollectionTag"        )),
    _shTag     (pset.get<string>       ("strawHitCollectionTag"        )),
    _chTag     (pset.get<string>       ("comboHitCollectionTag"        )),
    // _sthTag    (pset.get<string>       ("stereoHitCollectionTag"       )),
    _shfTag    (pset.get<string>       ("strawHitFlagCollectionTag"    )),
    _shpTag    (pset.get<string>       ("strawHitPositionCollectionTag")),
    _sdTag     (pset.get<art::InputTag>("strawDigiCollection"          ))
  {}

  PrefetchData::~PrefetchData() {
  }

//-----------------------------------------------------------------------------
  void PrefetchData::beginJob() {
  }

//-----------------------------------------------------------------------------
  void PrefetchData::fake_access(const CaloDigi& Hit) {
  }

  void PrefetchData::fake_access(const CaloHit& Hit) {
  }

//-----------------------------------------------------------------------------
  void PrefetchData::fake_access(const StrawHit& Hit,
                                 // const StrawHitFlag& Flag,
                                 const StrawHitPosition& Pos) {
  }

//-----------------------------------------------------------------------------
  void PrefetchData::fake_access(const ComboHit& Hit) {
  }

//-----------------------------------------------------------------------------
  // void PrefetchData::fake_access(const StereoHit& Hit) {
  // }

//-----------------------------------------------------------------------------
  void PrefetchData::fake_access(const StepPointMC* Step) {
  }

//-----------------------------------------------------------------------------
  void PrefetchData::fake_access(const StrawDigi& X) {
  }

//-----------------------------------------------------------------------------
  bool PrefetchData::findData(const art::Event& evt){
    _cdcol   = 0;
    _cchcol   = 0;

    _shcol   = 0;
    // _sthcol  = 0;
    _shfcol  = 0;
    _shpcol  = 0;
    _sdcol   = 0;
    _chcol   = 0;

    if (_fetchCaloDigis) {
      auto cdH = evt.getValidHandle<CaloDigiCollection>(_cdTag);
      _cdcol = cdH.product();
    }

    if (_fetchCaloHits) {
      auto cchH = evt.getValidHandle<CaloHitCollection>(_cchTag);
      _cchcol = cchH.product();
    }

    if (_fetchStrawHits) {
      auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
      _shcol = shH.product();
    }

    if (_fetchComboHits) {
      auto chH = evt.getValidHandle<ComboHitCollection>(_chTag);
      _chcol = chH.product();
    }

    if (_fetchStrawHitFlags) {
      auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
      _shfcol = shfH.product();
    }

    if (_fetchStrawHitPositions) {
      auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shpTag);
      _shpcol = shpH.product();
    }
    // if (_fetchStereoHits) {
    //   auto sthH = evt.getValidHandle<StereoHitCollection>(_sthTag);
    //   _sthcol = sthH.product();
    // }


    if (_fetchStrawDigis) {
      auto sdH = evt.getValidHandle<StrawDigiCollection>(_sdTag);
      _sdcol = sdH.product();
    }
//-----------------------------------------------------------------------------
// prefetch data
//-----------------------------------------------------------------------------
    if (_cdcol) {
      int ncd = _cdcol->size();
      for(int i=0;i<ncd;++i){
        const CaloDigi& cd = _cdcol->at(i);
        fake_access(cd);
      }
    }

    if (_cchcol) {
      int ncch = _cchcol->size();
      for(int i=0;i<ncch;++i){
        const CaloHit& cch = _cchcol->at(i);
        fake_access(cch);
      }
    }

    if (_sdcol) {
      int nsd = _sdcol->size();
      for(int i=0;i<nsd;++i){
        const StrawDigi& sdigi = _sdcol->at(i);
        fake_access(sdigi);
      }
    }

    if (_shcol) {
      int nsh = _shcol->size();
      for(int ish=0;ish<nsh;++ish){
        const StrawHit& sh          = _shcol->at (ish);
        const StrawHitPosition& shp = _shpcol->at(ish);
        // const StrawHitFlag& shf     = _shfcol->at(ish);

         fake_access(sh, /*shf,*/ shp);
      }
    }

    if (_chcol) {
      int nch = _chcol->size();
      for(int ish=0;ish<nch;++ish){
        const ComboHit& ch          = _chcol->at (ish);
         fake_access(ch);
      }
    }

    // if (_sthcol) {
    //   int nsth = _sthcol->size();
    //   for(int i=0;i<nsth;++i){
    //         const StereoHit& sth = _sthcol->at (i);
    //         fake_access(sth);
    //   }
    // }

    // if(_mcDiag){
    //   auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_sdTag);
    //   _mcdigis = mcdH.product();
    //   if (_mcdigis) {
    //         int nd = _mcdigis->size();
    //         for(int i=0;i<nd;++i){
    //           const StrawDigiMC& mcdigi = _mcdigis->at(i);
    //           fake_access(mcdigi);
    //         }
    //   }
    // }

    return true;
  }

//-----------------------------------------------------------------------------
  void PrefetchData::produce(art::Event& Event) {

    _eventNum = Event.event();

    if (_debugLevel > 0) printf(">>> PrefetchData::produce event number: %10i\n",_eventNum);

    findData(Event);
  }

//------------------------------------------------------------------------------
// Part of the magic that makes this class a module.
DEFINE_ART_MODULE(PrefetchData)

}



// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
// root 
#include "TMath.h"
#include "TH1F.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2F.h"
#include "TVector2.h"
// data
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawDigi.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
// Utilities
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
// diagnostics
#include <algorithm>
#include <cmath>
#include "CLHEP/Vector/ThreeVector.h"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
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
    int           _fetchStrawHits;
    int           _fetchComboHits;
    int           _fetchStrawHitFlags;
    int           _fetchStrawHitPositions;
    // int           _fetchStereoHits;
    int           _fetchStrawDigis;
					// data tags
    art::InputTag _cdTag;

    art::InputTag _shTag;
    art::InputTag _chTag;
    art::InputTag _sthTag;
    art::InputTag _shfTag;
    art::InputTag _shpTag;
    art::InputTag _sdTag;
					// cache of event objects
    const CaloDigiCollection*                   _cdcol;

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
    _fetchStrawHits (pset.get<int>         ("fetchStrawHits" )),
    _fetchComboHits (pset.get<int>         ("fetchComboHits" )),
    _fetchStrawHitFlags (pset.get<int>     ("fetchStrawHitFlags" )),
    _fetchStrawHitPositions (pset.get<int> ("fetchStrawHitPositions" )),
    // _fetchStereoHits (pset.get<int>        ("fetchStereoHits" )),
    _fetchStrawDigis(pset.get<int>         ("fetchStrawDigis")),

    _cdTag     (pset.get<string>       ("caloDigiCollectionTag"        )),
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
    // 	const StereoHit& sth = _sthcol->at (i);
    // 	fake_access(sth);
    //   }
    // }

    // if(_mcDiag){
    //   auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_sdTag);
    //   _mcdigis = mcdH.product();
    //   if (_mcdigis) {
    // 	int nd = _mcdigis->size();
    // 	for(int i=0;i<nd;++i){
    // 	  const StrawDigiMC& mcdigi = _mcdigis->at(i);
    // 	  fake_access(mcdigi);
    // 	}
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


   

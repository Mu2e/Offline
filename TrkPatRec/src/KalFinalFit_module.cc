//
// Final Kalman fit
//
// Original author D. Brown and G. Tassielli
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "GeometryService/inc/GeometryService.hh"
// utiliites
#include "GeneralUtilities/inc/Angles.hh"
#include "TrkReco/inc/TrkUtilities.hh"
// data
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/KalSeedCollection.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
// Mu2e BaBar
#include "BTrkData/inc/TrkStrawHit.hh"
#include "TrkReco/inc/KalFit.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/KalRepPayloadCollection.hh"
#include "TrkPatRec/inc/PayloadSaver.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
// root 
#include "TH1F.h"
#include "TTree.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <float.h>
#include <vector>
using namespace std; 
using CLHEP::Hep3Vector;
using CLHEP::HepVector;

namespace mu2e 
{
  class KalFinalFit : public art::EDProducer
  {
    public:
      explicit KalFinalFit(fhicl::ParameterSet const&);
      virtual ~KalFinalFit();
      virtual void produce(art::Event& event ); 
    private:
      unsigned _iev;
      // configuration parameters
      int _debug;
      int _printfreq;
      bool _addhits; 
      art::Handle<mu2e::StrawHitCollection> _strawhitsH;
      art::Handle<KalSeedCollection> _trksSeedH;
      // event object tags
      art::InputTag _shTag;
      art::InputTag _shfTag;
      art::InputTag _ksTag; 
      // flags
      StrawHitFlag _addsel;
      StrawHitFlag _addbkg;
      TrkFitFlag _goodseed;
      double _maxdtmiss;
      // outlier cuts
      double _maxadddoca, _maxaddchi;
      TrkParticle _tpart; // particle type being searched for
      TrkFitDirection _fdir;  // fit direction in search
      // event objects
      const StrawHitCollection* _shcol;
      const StrawHitFlagCollection* _shfcol;
      const KalSeedCollection * _kscol;
      // Kalman fitter
      KalFit _kfit;
      PayloadSaver _payloadSaver;
      // helper functions
      bool findData(const art::Event& e);
      void findMissingHits(KalRep* krep, vector<hitIndex>& indices);

      // flow diagnostic
  };

  KalFinalFit::KalFinalFit(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",101)),
    _addhits(pset.get<bool>("addhits",true)),
    _shTag(pset.get<art::InputTag>("StrawHitCollectionTag","makeSH")),
    _shfTag(pset.get<art::InputTag>("StrawHitFlagCollectionTag","FlagBkgHits")),
    _ksTag(pset.get<art::InputTag>("SeedCollectionTag","KalSeedFit")),
    _addsel(pset.get<vector<string> >("AddHitSelectionBits",vector<string>{} )),
    _addbkg(pset.get<vector<string> >("AddHitBackgroundBits",vector<string>{})),
    _goodseed(pset.get<vector<string> >("GoodKalSeedFitBits",vector<string>{})),
    _maxdtmiss(pset.get<double>("DtMaxMiss",40.0)),
    _maxadddoca(pset.get<double>("MaxAddDoca",2.75)),
    _maxaddchi(pset.get<double>("MaxAddChi",4.0)),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit",fhicl::ParameterSet())),
    _payloadSaver(pset)
  {
    produces<KalRepCollection>();
    produces<KalRepPtrCollection>();
    produces<KalRepPayloadCollection>();
    produces<StrawHitFlagCollection>();
    produces<KalSeedCollection>();
  }

  KalFinalFit::~KalFinalFit(){}

  void KalFinalFit::produce(art::Event& event ) {
    // event printout
    _iev=event.id().event();
    if(_debug > 0 && (_iev%_printfreq)==0)cout<<"KalFinalFit: event="<<_iev<<endl;
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::KalFinalFit: data missing or incomplete"<< endl;
    }
    // create output
    unique_ptr<KalRepCollection>    krcol(new KalRepCollection );
    unique_ptr<KalRepPtrCollection> krPtrcol(new KalRepPtrCollection );
    unique_ptr<KalSeedCollection> kscol(new KalSeedCollection());
    // copy in the existing flags
    unique_ptr<StrawHitFlagCollection> shfcol(new StrawHitFlagCollection(*_shfcol));
    // lookup productID for payload saver
    art::ProductID kalRepsID(getProductID<KalRepCollection>(event));
    // loop over the seed fits
    unsigned iseed(0);
    for (auto kseed : *_kscol) {
    // only process fits which meet the requirements
      if(kseed.status().hasAllProperties(_goodseed)) {
	// check the seed has the same basic parameters as this module expects
	if(kseed.particle() != _tpart || kseed.fitDirection() != _fdir ) {
	  throw cet::exception("RECO")<<"mu2e::KalFinalFit: wrong particle or direction"<< endl;
	}
	// seed should have at least 1 segment
	if(kseed.segments().size() < 1){
	  throw cet::exception("RECO")<<"mu2e::KalFinalFit: no segments"<< endl;
	}
	// build a Kalman rep around this seed
	KalRep *krep(0);
	_kfit.makeTrack(_shcol,kseed,krep);
	// if successfull, try to add missing hits
	if(_addhits && krep != 0 && krep->fitStatus().success()){
	    // first, add back the hits on this track
	  _kfit.unweedHits(krep,_maxaddchi);
	  vector<hitIndex> misshits;
	  findMissingHits(krep,misshits);
	  if(misshits.size() > 0){
	    _kfit.addHits(krep,_shcol,misshits,_maxaddchi);
	  }
	}
	// put successful fits into the event
	if(krep != 0 && krep->fitStatus().success()){
	  // warning about 'fit current': this is not an error
	  if(!krep->fitCurrent()){
	    cout << "Fit not current! " << endl;
	  }
	  // save successful kalman fits in the event
	  krcol->push_back(krep);
	  int index = krcol->size()-1;
	  krPtrcol->emplace_back(kalRepsID, index, event.productGetter(kalRepsID));
	  // convert successful fits into 'seeds' for persistence
	  // FIXME!
	  // flag the hits used in this track. 
	  if(iseed<StrawHitFlag::_maxTrkId){
	    for(auto ihit=krep->hitVector().begin();ihit != krep->hitVector().end();++ihit){
	      if((*ihit)->isActive())shfcol->at(static_cast<TrkStrawHit*>(*ihit)->index()).merge(StrawHitFlag::trackBit(iseed));
	    }
	  }
	} else // fit failure
	  delete krep;
      }
      ++iseed;
    }
    // put the output products into the event
    art::ProductID krcolID(getProductID<KalRepPayloadCollection>(event));
    _payloadSaver.put(*krcol, krcolID, event);
    event.put(move(krcol));
    event.put(move(krPtrcol));
    event.put(move(kscol));
    event.put(move(shfcol));
  }

  // find the input data objects 
  bool KalFinalFit::findData(const art::Event& evt){
    _shcol = 0;
    _shfcol = 0;
    _kscol = 0;

    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
    _shfcol = shfH.product();
    auto ksH = evt.getValidHandle<KalSeedCollection>(_ksTag);
    _kscol = ksH.product();

    return _shcol != 0 && _shfcol != 0 && _kscol != 0;
  }

  void KalFinalFit::findMissingHits(KalRep* krep,vector<hitIndex>& misshits) {
    const Tracker& tracker = getTrackerOrThrow();
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint tpos;
    krep->pieceTraj().getInfo(krep->flt0(),tpos,tdir);
    unsigned nstrs = _shcol->size();
    TrkStrawHitVector tshv;
    convert(krep->hitVector(),tshv);
    for(unsigned istr=0; istr<nstrs;++istr){
      if(_shfcol->at(istr).hasAllProperties(_addsel)&& !_shfcol->at(istr).hasAnyProperty(_addbkg)){
	StrawHit const& sh = _shcol->at(istr);
	if(fabs(_shcol->at(istr).time()-krep->t0()._t0) < _maxdtmiss) {
	  // make sure we haven't already used this hit
	  vector<TrkStrawHit*>::iterator ifnd = find_if(tshv.begin(),tshv.end(),FindTrkStrawHit(sh));
	  if(ifnd == tshv.end()){
	    // good in-time hit.  Compute DOCA of the wire to the trajectory
	    Straw const& straw = tracker.getStraw(sh.strawIndex());
	    CLHEP::Hep3Vector hpos = straw.getMidPoint();
	    CLHEP::Hep3Vector hdir = straw.getDirection();
	    // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
	    HepPoint spt(hpos.x(),hpos.y(),hpos.z());
	    TrkLineTraj htraj(spt,hdir,-straw.getHalfLength(),straw.getHalfLength());
	    // estimate flightlength along track.  This assumes a constant BField!!!
	    double fltlen = (hpos.z()-tpos.z())/tdir.z();
	    // estimate hit length
	    HepPoint tp = krep->pieceTraj().position(fltlen);
	    Hep3Vector tpos(tp.x(),tp.y(),tp.z()); // ugly conversion FIXME!
	    double hitlen = hdir.dot(tpos - hpos);
	    TrkPoca hitpoca(krep->pieceTraj(),fltlen,htraj,hitlen);
	    // flag hits with small residuals
	    if(fabs(hitpoca.doca()) < _maxadddoca){
	      misshits.push_back(istr);
	    }
	  }
	}
      }
    }
  }

}// mu2e
using mu2e::KalFinalFit;
DEFINE_ART_MODULE(KalFinalFit);

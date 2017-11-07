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
#include "Mu2eUtilities/inc/MVATools.hh"
// data
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
// Mu2e BaBar
#include "BTrkData/inc/TrkStrawHit.hh"
#include "TrkReco/inc/KalFit.hh"
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
      bool _saveall,_addhits; 
      vector<double> _zsave;
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
      // trkqual calculation
      std::unique_ptr<MVATools> _trkqualmva;
      // event objects
      const StrawHitCollection* _shcol;
      const StrawHitFlagCollection* _shfcol;
      const KalSeedCollection * _kscol;
      // Kalman fitter
      KalFit _kfit;
      // helper functions
      bool findData(const art::Event& e);
      void findMissingHits(KalRep* krep, vector<StrawHitIndex>& indices);
      void fillTrkQual(KalSeed const& kseed, TrkQual& trkqual);

      // flow diagnostic
  };

  KalFinalFit::KalFinalFit(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",101)),
    _saveall(pset.get<bool>("saveall",false)),
    _addhits(pset.get<bool>("addhits",true)),
    _zsave(pset.get<vector<double> >("ZSavePositions",vector<double>{-1522.0,0.0,1522.0})), // front, middle and back of the tracker
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
    _kfit(pset.get<fhicl::ParameterSet>("KalFit",fhicl::ParameterSet()))
  {
    fhicl::ParameterSet mvapset = pset.get<fhicl::ParameterSet>("TrkQualMVA",fhicl::ParameterSet());
    mvapset.put<string>("MVAWeights",pset.get<string>("TrkQualWeights","TrkDiag/test/TrkQual.weights.xml"));
    _trkqualmva.reset(new MVATools(mvapset));
    _trkqualmva->initMVA();
    if(_debug>0)_trkqualmva->showMVA();

    produces<KalRepCollection>();
    produces<KalRepPtrCollection>();
    produces<StrawHitFlagCollection>();
    produces<KalSeedCollection>();
    produces<TrkQualCollection>();
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
    unique_ptr<TrkQualCollection> tqcol(new TrkQualCollection());
    
    // copy in the existing flags
    unique_ptr<StrawHitFlagCollection> shfcol(new StrawHitFlagCollection(*_shfcol));
    // lookup productID for payload saver
    art::ProductID kalRepsID(getProductID<KalRepCollection>(event));
    // loop over the seed fits.  I need an index loop here to build the Ptr
    for(size_t ikseed=0; ikseed < _kscol->size(); ++ikseed) {
      const auto& kseed = _kscol->at(ikseed);
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
	if(_debug > 1){
	  if(krep == 0)
	    cout << "No Final fit produced " << endl;
	  else{
	    cout << "Seed Fit HelixTraj parameters " << krep->seedTrajectory()->parameters()->parameter()
	      << " covariance " << krep->seedTrajectory()->parameters()->covariance()
	      << " NDOF = " << krep->nDof()
	      << " Final Fit status " << krep->fitStatus()  << endl;
	  }
	}
	// if successfull, try to add missing hits
	if(_addhits && krep != 0 && krep->fitStatus().success()){
	    // first, add back the hits on this track
	  _kfit.unweedHits(krep,_maxaddchi);
	  vector<StrawHitIndex> misshits;
	  findMissingHits(krep,misshits);
	  if(misshits.size() > 0){
	    _kfit.addHits(krep,_shcol,misshits,_maxaddchi);
	  }
	  if(_debug > 1)
	    cout << "AddHits Fit result " << krep->fitStatus()
	    << " NDOF = " << krep->nDof() << endl;
	}
	// put successful fits into the event
	if(krep != 0 && (krep->fitStatus().success() || _saveall)){
	  // warning about 'fit current': this is not an error
	  if(!krep->fitCurrent()){
	    cout << "Fit not current! " << endl;
	  }
	  // flg all hits as belonging to this track
	  // // this is probably obsolete, FIXME!
	  if(ikseed<StrawHitFlag::_maxTrkId){
	    for(auto ihit=krep->hitVector().begin();ihit != krep->hitVector().end();++ihit){
	      if((*ihit)->isActive())shfcol->at(static_cast<TrkStrawHit*>(*ihit)->index()).merge(StrawHitFlag::trackBit(ikseed));
	    }
	  }
	  // save successful kalman fits in the event
	  krcol->push_back(krep);
	  int index = krcol->size()-1;
	  krPtrcol->emplace_back(kalRepsID, index, event.productGetter(kalRepsID));
	  // convert successful fits into 'seeds' for persistence
	  KalSeed fseed(_tpart,_fdir,krep->t0(),krep->flt0(),kseed.status());
	  // reference the seed fit in this fit
	  auto ksH = event.getValidHandle<KalSeedCollection>(_ksTag);
	  fseed._kal = art::Ptr<KalSeed>(ksH,ikseed);
	  // redundant but possibly useful
	  fseed._helix = kseed.helix();
	  // fill with new information
	  fseed._t0 = krep->t0();
	  fseed._flt0 = krep->flt0();
	  fseed._status.merge(TrkFitFlag::kalmanOK);
	  // global fit information
	  fseed._chisq = krep->chisq();
	  // compute the fit consistency.  Note our fit has effectively 6 parameters as t0 is allowed to float and its error is propagated to the chisquared
	  fseed._fitcon =  ChisqConsistency(krep->chisq(),krep->nDof()-1).significanceLevel();
	  if(krep->fitStatus().success()==1) fseed._status.merge(TrkFitFlag::kalmanConverged);
	  TrkUtilities::fillHitSeeds(krep,fseed._hits);
	  TrkUtilities::fillStraws(krep,fseed._straws);
	  // sample the fit at the requested z positions.  Need options here to define a set of
	  // standard points, or to sample each unique segment on the fit FIXME!
	  for(auto zpos : _zsave) {
	    // compute the flightlength for this z
	    double fltlen = krep->pieceTraj().zFlight(zpos);
	    // sample the momentum at this flight.  This belongs in a separate utility FIXME
	    BbrVectorErr momerr = krep->momentumErr(fltlen);
	    // sample the helix
	    double locflt(0.0);
	    const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(krep->localTrajectory(fltlen,locflt));
	    // fill the segment
	    KalSegment kseg;
	    TrkUtilities::fillSegment(*htraj,momerr,kseg);
	    fseed._segments.push_back(kseg);
	  }
	  // save KalSeed for this track
	  kscol->push_back(fseed);
	  // compute TrkQual for this track and save it
	  TrkQual trkqual;
	  fillTrkQual(fseed,trkqual);
	  tqcol->push_back(trkqual);
	} else // fit failure
	  delete krep;
      }
    }
    // put the output products into the event
    event.put(move(krcol));
    event.put(move(krPtrcol));
    event.put(move(kscol));
    event.put(move(shfcol));
    event.put(move(tqcol));
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

  void KalFinalFit::findMissingHits(KalRep* krep,vector<StrawHitIndex>& misshits) {
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

  void KalFinalFit::fillTrkQual(KalSeed const& kseed, TrkQual& trkqual) {
    static StrawHitFlag active(StrawHitFlag::active);
    static TrkFitFlag goodfit(TrkFitFlag::kalmanOK);
    if(kseed.status().hasAllProperties(goodfit)){
      std::vector<TrkStrawHitSeed> const& hits = kseed.hits();
      unsigned nactive(0), ndouble(0), nnull(0);
      for(auto ihit = hits.begin(); ihit != hits.end(); ++ihit){
	if(ihit->flag().hasAllProperties(active)){
	  ++nactive;
	  if(ihit->ambig()==0)++nnull;
	  // look at the adjacent hits; if they are in the same panel, this is a double
	  auto jhit = ihit; ++jhit;
	  auto hhit = ihit; --hhit;
	  if( (jhit != hits.end() && 
		jhit->flag().hasAllProperties(active) &&
		jhit->strawId().getPlane() == ihit->strawId().getPlane() &&
		jhit->strawId().getPanel() == ihit->strawId().getPanel() ) ||
	      (hhit >= hits.begin() && 
	       hhit->flag().hasAllProperties(active) &&
	       hhit->strawId().getPlane() == ihit->strawId().getPlane() &&
	       hhit->strawId().getPanel() == ihit->strawId().getPanel() )
	    ) {
	    ++ndouble;
	  }
	}
      }

      trkqual[TrkQual::nactive] = nactive;
      trkqual[TrkQual::factive] = (float)nactive/(float)hits.size();  // Fraction of active hits
      trkqual[TrkQual::log10fitcon] = kseed.fitConsistency() > FLT_MIN ? log10(kseed.fitConsistency()) : -50.0; // fit chisquared consistency
      trkqual[TrkQual::t0err] = kseed.t0().t0Err();  // estimated t0 error
      trkqual[TrkQual::fdouble] = (float)ndouble/(float)nactive;  // fraction of double hits (2 or more in 1 panel)
      trkqual[TrkQual::fnullambig] = (float)nnull/(float)nactive;  // fraction of hits with null ambiguity
      trkqual[TrkQual::fstraws] = (float)kseed.straws().size()/(float)nactive;  // fraction of straws to hits

      // find the fit segment that best matches the location for testing the quality
      std::vector<KalSegment> const& ksegs = kseed.segments();
      auto bestkseg = ksegs.begin();
      for(auto ikseg = ksegs.begin(); ikseg != ksegs.end(); ++ikseg){ 
	HelixVal const& hel = ikseg->helix();
	// check for a segment whose range includes z=0.  There should be a better way of doing this, FIXME	
	double sind = hel.tanDip()/sqrt(1.0+hel.tanDip()*hel.tanDip());
	if(hel.z0()+sind*ikseg->fmin() < 0.0 && hel.z0()+sind*ikseg->fmax() > 0.0){
	  bestkseg = ikseg;
	  break;
	}
      }
      if(bestkseg != ksegs.end()){
	trkqual[TrkQual::momerr] = bestkseg->momerr(); // estimated momentum error
	trkqual[TrkQual::d0] = bestkseg->helix().d0(); // d0 value
	trkqual[TrkQual::rmax] = bestkseg->helix().d0()+2.0/bestkseg->helix().omega(); // maximum radius of fit
	// calculate the MVA
	trkqual.setMVAValue(_trkqualmva->evalMVA(trkqual.values()));
	trkqual.setMVAStatus(TrkQual::calculated);
      } else {
	trkqual.setMVAStatus(TrkQual::filled);
      } 
    } else {
      trkqual.setMVAStatus(TrkQual::failed);
    }
  }
}// mu2e
using mu2e::KalFinalFit;
DEFINE_ART_MODULE(KalFinalFit);

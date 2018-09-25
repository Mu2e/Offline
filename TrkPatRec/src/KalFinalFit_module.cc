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
#include "art/Utilities/make_tool.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerConditions/inc/StrawResponse.hh"

#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "GeometryService/inc/GeometryService.hh"
// utiliites
#include "GeneralUtilities/inc/Angles.hh"
#include "TrkReco/inc/TrkUtilities.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
// data
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/AlgorithmIDCollection.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "TrkReco/inc/KalFitData.hh"
#include "TrkPatRec/inc/KalFinalFit_types.hh"
#include "TrkReco/inc/DoubletAmbigResolver.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
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
  using namespace KalFinalFitTypes;

  class KalFinalFit : public art::EDProducer
  {
  public:
    explicit KalFinalFit(fhicl::ParameterSet const&);
    virtual ~KalFinalFit();
    void beginRun(art::Run& aRun);
  private:
    void produce(art::Event& event) override;

    unsigned _iev;
    // configuration parameters
    int _debug;
    int _diag;
    int _printfreq;
    int _cprmode;
    bool _saveall,_addhits;
    vector<double> _zsave;
    // event object tokens
    art::ProductToken<ComboHitCollection> const _shToken;
    art::InputTag const _shfTag;
    art::ProductToken<StrawHitFlagCollection> const _shfToken;
    art::ProductToken<KalSeedCollection> const _ksToken;
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
    const ComboHitCollection* _chcol;
    const StrawHitFlagCollection* _shfcol;
    const KalSeedCollection * _kscol;
    // Kalman fitter
    KalFit _kfit;
    KalFitData _result;

    // diagnostic
    Data_t                                _data;
    std::unique_ptr<ModuleHistToolBase>   _hmanager;

    // helper functions
    bool findData(const art::Event& e);
    void findMissingHits(KalFitData&kalData);
    void findMissingHits_cpr(KalFitData&kalData);
    void fillTrkQual(KalSeed const& kseed, TrkQual& trkqual);

    // flow diagnostic
  };

  KalFinalFit::KalFinalFit(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel", 0)),
    _diag(pset.get<int>("diagLevel",0)),
    _printfreq(pset.get<int>("printFrequency", 101)),
    _cprmode(pset.get<int>("cprmode",0)),
    _saveall(pset.get<bool>("saveall", false)),
    _addhits(pset.get<bool>("addhits", true)),
    _zsave(pset.get<vector<double>>("ZSavePositions", vector<double>{-1522.0,0.0,1522.0})), // front, middle and back of the tracker
    _shToken{consumes<ComboHitCollection>(pset.get<art::InputTag>("ComboHitCollection"))},
    _shfTag{pset.get<art::InputTag>("StrawHitFlagCollection", "none")},
    _shfToken{consumes<StrawHitFlagCollection>(_shfTag)},
    _ksToken{consumes<KalSeedCollection>(pset.get<art::InputTag>("SeedCollection"))},
    _addsel(pset.get<vector<string>>("AddHitSelectionBits", vector<string>{})),
    _addbkg(pset.get<vector<string>>("AddHitBackgroundBits", vector<string>{})),
    _goodseed(pset.get<vector<string>>("GoodKalSeedFitBits", vector<string>{})),
    _maxdtmiss(pset.get<double>("DtMaxMiss",40.0)),
    _maxadddoca(pset.get<double>("MaxAddDoca",2.75)),
    _maxaddchi(pset.get<double>("MaxAddChi",4.0)),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle", TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection", TrkFitDirection::downstream))),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit", {})),
    _result()
  {
    auto mvapset = pset.get<fhicl::ParameterSet>("TrkQualMVA", {});
    mvapset.put<string>("MVAWeights",pset.get<string>("TrkQualWeights", "TrkDiag/test/TrkQual.weights.xml"));
    _trkqualmva.reset(new MVATools(mvapset));
    _trkqualmva->initMVA();
    if(_debug>0)_trkqualmva->showMVA();

    produces<KalRepCollection>();
    produces<KalRepPtrCollection>();
    produces<AlgorithmIDCollection>();
    produces<StrawHitFlagCollection>();
    produces<KalSeedCollection>();
    produces<TrkQualCollection>();
//-----------------------------------------------------------------------------
// provide for interactive disanostics
//-----------------------------------------------------------------------------
    _data.result    = &_result;
    
    if (_diag != 0) {
      _hmanager = art::make_tool<ModuleHistToolBase>(pset.get<fhicl::ParameterSet>("diagPlugin"));
      fhicl::ParameterSet ps1 = pset.get<fhicl::ParameterSet>("Fitter.DoubletAmbigResolver");
      _data.dar               = new DoubletAmbigResolver(ps1,0,0,0);
      _data.listOfDoublets    = new std::vector<Doublet>;
    }
    else {
      _hmanager = std::make_unique<ModuleHistToolBase>();
      _data.dar            = NULL;
      _data.listOfDoublets = NULL;
    }
 }

  KalFinalFit::~KalFinalFit(){
    if (_data.dar) {
      delete _data.listOfDoublets;
      delete _data.dar;
    }
  }
//-----------------------------------------------------------------------------
  void KalFinalFit::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::TTracker> th;
    _data.tracker     = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _data.calorimeter = ch.get();
    
    _kfit.setCalorimeter (_data.calorimeter);
    _kfit.setTracker     (_data.tracker);
  }


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
    unique_ptr<AlgorithmIDCollection>  algs     (new AlgorithmIDCollection   );
    unique_ptr<KalSeedCollection> kscol(new KalSeedCollection());
    unique_ptr<TrkQualCollection> tqcol(new TrkQualCollection());
    unique_ptr<StrawHitFlagCollection> shfcol(new StrawHitFlagCollection());
    // lookup productID for payload saver
    art::ProductID kalRepsID(getProductID<KalRepCollection>());
    // copy and merge hit flags
    size_t index(0);
    for(auto const& ch : *_chcol) {
      StrawHitFlag flag(ch.flag());
      if(_shfcol != 0) flag.merge(_shfcol->at(index++));
      shfcol->push_back(flag);
    }
 
    if (_diag){
      _data.event  = &event;
      _data.eventNumber = event.event();
      _data.result = &_result;
      _data.tracks = krcol.get();
    }

    _result.fitType     = 0;
    _result.event       = &event ;
    _result.chcol       = _chcol ;
    _result.shfcol      = _shfcol ;
    _result.tpart       = _tpart ;
    _result.fdir        = _fdir  ;

    // loop over the seed fits.  I need an index loop here to build the Ptr
    for(size_t ikseed=0; ikseed < _kscol->size(); ++ikseed) {
      KalSeed const& kseed(_kscol->at(ikseed));
      _result.kalSeed = & kseed;

      if (kseed.caloCluster()) _result.caloCluster = kseed.caloCluster().get();

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
	//fill the KalFitData variable
	// _result.kalSeed = &kseed;

	// _kfit.makeTrack(_shcol,kseed,krep);
	_result.init();
	_kfit.makeTrack(_result);

	// KalRep *krep = _result.stealTrack();

	if(_debug > 1){
	  if(_result.krep == 0)
	    cout << "No Final fit produced " << endl;
	  else{
	    cout << "Seed Fit HelixTraj parameters " << _result.krep->seedTrajectory()->parameters()->parameter()
	      << " covariance " << _result.krep->seedTrajectory()->parameters()->covariance()
	      << " NDOF = " << _result.krep->nDof()
	      << " Final Fit status " << _result.krep->fitStatus()  << endl;
	  }
	}
	// if successfull, try to add missing hits
	if(_addhits && _result.krep != 0 && _result.krep->fitStatus().success()){
	    // first, add back the hits on this track
	  _result.nunweediter = 0;
	  _kfit.unweedHits(_result,_maxaddchi);
	  if (_debug > 0) _kfit.printHits(_result,"CalTrkFit::produce after unweedHits");

	  if (_cprmode){
	    findMissingHits_cpr(_result);
	  }else {
	    findMissingHits(_result);
	  }

	  if(_result.missingHits.size() > 0){
	    _kfit.addHits(_result,_maxaddchi);
	  }else if (_cprmode){
	    int last_iteration  = -1;
	    _kfit.fitIteration(_result,last_iteration);
	  }
	  if(_debug > 1)
	    cout << "AddHits Fit result " << _result.krep->fitStatus()
	    << " NDOF = " << _result.krep->nDof() << endl;
	  
//-----------------------------------------------------------------------------
// and weed hits again to insure that addHits doesn't add junk
//-----------------------------------------------------------------------------
	  int last_iteration  = -1;
	  if (_cprmode) _kfit.weedHits(_result,last_iteration);
	}
	// put successful fits into the event
	if(_result.krep != 0 && (_result.krep->fitStatus().success() || _saveall)){
//-----------------------------------------------------------------------------
// now evaluate the T0 and its error using the straw hits
//-----------------------------------------------------------------------------
	  if (_cprmode)	_kfit.updateT0(_result);

	  // warning about 'fit current': this is not an error
	  if(!_result.krep->fitCurrent()){
	    cout << "Fit not current! " << endl;
	  }
	  // flg all hits as belonging to a track
	  if(ikseed<StrawHitFlag::_maxTrkId){
	    for(auto ihit=_result.krep->hitVector().begin();ihit != _result.krep->hitVector().end();++ihit){
	      if((*ihit)->isActive())shfcol->at(static_cast<TrkStrawHit*>(*ihit)->index()).merge(StrawHitFlag::track);
	    }
	  }
	  
	  
	  // save successful kalman fits in the event
	  KalRep *krep = _result.stealTrack();
	  krcol->push_back(krep);

	  // save the alorithm bit
	  int best(1),mask(1);
	  if (_cprmode==0) {
	    best = AlgorithmID::TrkPatRecBit;
	    mask = 1 << AlgorithmID::TrkPatRecBit;
	  } else if (_cprmode==1) {
	    best = AlgorithmID::CalPatRecBit;
	    mask = 1 << AlgorithmID::CalPatRecBit;
	  }
	  algs->push_back(AlgorithmID(best,mask));


	  int index = krcol->size()-1;
	  krPtrcol->emplace_back(kalRepsID, index, event.productGetter(kalRepsID));
	  // convert successful fits into 'seeds' for persistence
	  KalSeed fseed(_tpart,_fdir,krep->t0(),krep->flt0(),kseed.status());
	  // reference the seed fit in this fit
	  auto ksH = event.getValidHandle<KalSeedCollection>(_ksToken);
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
	} else {// fit failure
	  _result.deleteTrack();
	  //	  delete krep;
	}
      }
    }

    if (_diag > 0) _hmanager->fillHistograms(&_data);

    // put the output products into the event
    event.put(move(krcol));
    event.put(move(krPtrcol));
    event.put(move(kscol));
    event.put(move(tqcol));
    event.put(move(algs));
    event.put(move(shfcol));
  }
  
  // find the input data objects
  bool KalFinalFit::findData(const art::Event& evt){
    _chcol = 0;
    _kscol = 0;

    auto shH = evt.getValidHandle(_shToken);
    _chcol = shH.product();
    auto ksH = evt.getValidHandle(_ksToken);
    _kscol = ksH.product();
    if(_shfTag.label() != "none"){
      auto shfH = evt.getValidHandle(_shfToken);
      _shfcol = shfH.product();
      if(_shfcol->size() != _chcol->size())
        throw cet::exception("RECO")<<"mu2e::KalFinalFit: inconsistent input collections"<< endl;
    } else {
      _shfcol = 0;
    }

    return _chcol != 0 && _kscol != 0;
  }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void KalFinalFit::findMissingHits_cpr(KalFitData& KRes) {

    const char* oname = "KalFinalFit::findMissingHits_cpr";

    Hep3Vector tdir;
    HepPoint   tpos;
    double     dt;

    KalRep* krep =  KRes.krep;

    krep->pieceTraj().getInfo(0.0,tpos,tdir);

    const TrkDifPieceTraj* reftraj = krep->referenceTraj();

    if (_debug > 0) printf("[%s]      shId    sec     panel       doca        drift     dr added\n",oname);

    KRes.missingHits.clear();
    //    KRes.doca.clear();

    MissingHit_t mh;

    int nstrs = KRes.chcol->size();
    for (int istr=0; istr<nstrs; ++istr) {
      mh.index = istr;
//----------------------------------------------------------------------
// 2015-02-11 change the selection bit for searching for missed hits
//----------------------------------------------------------------------
      ComboHit const& sh    = _chcol->at(istr);
//-----------------------------------------------------------------------------
// I think, we want to check the radial bit: if it is set, than at least one of
// the two measured times is wrong...
//-----------------------------------------------------------------------------
      //      int radius_ok = _shfcol->at(istr).hasAllProperties(StrawHitFlag::radsel);
      dt        = _chcol->at(istr).time()-KRes.krep->t0()._t0;

      //      if (radius_ok && (fabs(dt) < _maxdtmiss)) {
      if (fabs(dt) < _maxdtmiss) {
					// make sure we haven't already used this hit

	TrkStrawHit  *tsh, *closest(NULL);
	bool found = false;

	Straw const&      straw = _data.tracker->getStraw(sh.strawId());
	CLHEP::Hep3Vector hpos  = straw.getMidPoint();

	double            dz_max(1.e12) ; // closest_z(1.e12);
	double            zhit = hpos.z();

	for (std::vector<TrkHit*>::iterator it=krep->hitVector().begin(); it!=krep->hitVector().end(); it++) {
	  tsh = static_cast<TrkStrawHit*> (*it);
	  int tsh_index = tsh->index();
	  if (tsh_index == istr) {
	    found = true;
	    break;
	  }
					// check proximity in Z
          Straw const&  trk_straw = _data.tracker->getStraw(tsh->comboHit().strawId());
          double        ztrk      = trk_straw.getMidPoint().z();

	  double dz  = ztrk-zhit;
	  if (fabs(dz) < fabs(dz_max)) {
	    closest   = tsh;
	    dz_max    = dz;
	  }
	}

        if (! found) {
					// estimate trajectory length to hit 
	  double hflt = 0;
	  TrkHelixUtils::findZFltlen(*reftraj,zhit,hflt);

          // good in-time hit.  Compute DOCA of the wire to the trajectory
	  // also estimate the drift time if this hit were on the track to get hte hit radius

	  TrkT0 hitt0 = closest->hitT0();

	  double s    = closest->fltLen();
	  double mom  = krep->momentum(s).mag();
	  double beta = krep->particleType().beta(mom);

	  double tflt = (hflt-s)/(beta*CLHEP::c_light);
	  hitt0._t0  += tflt;

          CLHEP::Hep3Vector hdir  = straw.getDirection();
          // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
          HepPoint          spt(hpos.x(),hpos.y(),hpos.z());
          TrkLineTraj       htraj(spt,hdir,-20,20);
          // estimate flightlength along track.  This assumes a constant BField!!!

          double      fltlen = (hpos.z()-tpos.z())/tdir.z();

          TrkPoca     hitpoca(krep->pieceTraj(),fltlen,htraj,0.0);

	  double      rdrift;//, hit_error(0.2);

	  TrkStrawHit hit(sh,straw,istr,hitt0,hflt,1.,1.);//hit_error,1.,_maxadddoca,1.);
	  
	  ConditionsHandle<StrawResponse> srep = ConditionsHandle<StrawResponse>("ignored");
	  //	  ConditionsHandle<TrackerCalibrations> tcal("ignored");

	  double tdrift=hit.time()-hit.hitT0()._t0;

	  //	  tcal->TimeToDistance(straw.index(),tdrift,tdir,t2d);

// find the track direction at this hit
	  Hep3Vector tdir  = krep->traj().direction(fltlen);
	  Hep3Vector tperp = tdir - tdir.dot(straw.getDirection())*straw.getDirection();
	  double     phi   = tperp.theta();

	  rdrift = srep->driftTimeToDistance(straw.id(),tdrift,phi);

	  mh.doca   = hitpoca.doca();
	  if (mh.doca > 0) mh.dr = mh.doca-rdrift;
	  else             mh.dr = mh.doca+rdrift;
//-----------------------------------------------------------------------------
// flag hits with small residuals
//-----------------------------------------------------------------------------
	  int added (0);
          if (fabs(mh.dr) < _maxadddoca) {
            KRes.missingHits.push_back(mh);
	    //            KRes.doca.push_back(doca);
	    added = 1;
          }

          if (_debug > 0) {
            printf("[CalTrkFit::findMissingHits] %8i  %6i  %8i  %10.3f %10.3f %10.3f   %3i\n",
                   straw.id().asUint16(),
                   straw.id().getPlane(),
                   straw.id().getPanel(),
                   mh.doca, rdrift, mh.dr,added);
          }
        }
      }
      else {
	if (_debug > 0) {
	  printf("[%s] rejected hit: i, index, flag, dt: %5i %5i %s %10.3f\n",
		 oname,istr,sh.strawId().asUint16(),
		 KRes.shfcol->at(istr).hex().data(),-1.);
	}
      }
    }
    //-----------------------------------------------------------------------------
    // sort hits in accending |dr|, such that the first hit has the smallest |dr|
    //-----------------------------------------------------------------------------
    int nmiss = KRes.missingHits.size();
    for (int i1=0; i1<nmiss-1; i1++) {
      MissingHit_t* h1 = &KRes.missingHits[i1];
      for (int i2=i1+1; i2<nmiss; i2++) {
	MissingHit_t* h2 = &KRes.missingHits[i2];
	if (fabs(h1->dr) > fabs(h2->dr)) {
					// swap hits
	  MissingHit_t tmp = *h1;
	  *h1 = *h2;
	  *h2 = tmp;
	}
      }
    }
  }

  void KalFinalFit::findMissingHits(KalFitData&kalData) {
    KalRep* krep = kalData.krep;

    //clear the array
    kalData.missingHits.clear();

    const Tracker& tracker = getTrackerOrThrow();
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint tpos;
    krep->pieceTraj().getInfo(krep->flt0(),tpos,tdir);
    unsigned nstrs = _chcol->size();
    TrkStrawHitVector tshv;
    convert(krep->hitVector(),tshv);
    for(unsigned istr=0; istr<nstrs;++istr){
      if(_shfcol->at(istr).hasAllProperties(_addsel)&& !_shfcol->at(istr).hasAnyProperty(_addbkg)){
	ComboHit const& sh = _chcol->at(istr);
	if(fabs(_chcol->at(istr).time()-krep->t0()._t0) < _maxdtmiss) {
	  // make sure we haven't already used this hit
	  vector<TrkStrawHit*>::iterator ifnd = find_if(tshv.begin(),tshv.end(),FindTrkStrawHit(sh));
	  if(ifnd == tshv.end()){
	    // good in-time hit.  Compute DOCA of the wire to the trajectory
	    Straw const& straw = tracker.getStraw(sh.strawId());
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
	      MissingHit_t m;
	      m.index = istr;
	      m.doca  = hitpoca.doca();
	      // m.dr = ??;
	      kalData.missingHits.push_back(m);
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

DEFINE_ART_MODULE(mu2e::KalFinalFit);

//
// Starting with a Helix fit, perform a Least-squares fit (using the BTrk
// Kalman fit, appropriately configured) to produce an initial estimate of the parmeters and
// covariance for the final Kalman fit.  This fit uses wire positions only,
// not drift.
//
// Original author Dave Brown (LBNL) 31 Aug 2016
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/DetectorSystem.hh"
// utiliites
#include "GeneralUtilities/inc/Angles.hh"
#include "TrkReco/inc/TrkUtilities.hh"
// data
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/HelixSeedCollection.hh"
#include "RecoDataProducts/inc/KalSeedCollection.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
// BaBar
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
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
  class KalSeedFit : public art::EDProducer
  {
    public:
      explicit KalSeedFit(fhicl::ParameterSet const&);
      virtual ~KalSeedFit();
      virtual void beginRun(art::Run&);
      virtual void produce(art::Event& event ); 
    private:
      unsigned _iev;
      // configuration parameters
      int _debug;
      int _printfreq;
      bool _saveall;
      // event object tags
      art::InputTag _shTag;
      art::InputTag _shfTag;
      art::InputTag _hsTag;
      TrkFitFlag _seedflag; // helix fit flag
  unsigned _minnhits; // minimum # of hits
      double _maxdoca;      // outlier cut
      bool _foutliers; // filter hits far from the helix
      bool _fhoutliers; // filter hits found flagged as outliers in the helix fit
      TrkParticle _tpart; // particle type being searched for
      TrkFitDirection _fdir;  // fit direction in search
      vector<double> _perr; // diagonal parameter errors to use in the fit
      Helicity _helicity; // cached value of helicity expected for this fit
      double _amsign; // cached sign of angular momentum WRT the z axis 
      HepSymMatrix _hcovar; // cache of parameter error covariance matrix
      // cache of event objects
      const StrawHitCollection* _shcol;
      const StrawHitFlagCollection* _shfcol;
      const HelixSeedCollection * _hscol;
      // ouptut collections
      // Kalman fitter.  This will be configured for a least-squares fit (no material or BField corrections).
      KalFit _seedfit;

      const TTracker* _tracker;     // straw tracker geometry

      // helper functions
      bool findData(const art::Event& e);
      void filterOutliers(TrkDef& trkdef);
  };

  KalSeedFit::KalSeedFit(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",101)),
    _saveall(pset.get<bool>("saveall",false)),
    _shTag(pset.get<art::InputTag>("StrawHitCollectionTag","makeSH")),
    _shfTag(pset.get<art::InputTag>("StrawHitFlagCollectionTag","FlagBkgHits")),
    _hsTag(pset.get<art::InputTag>("SeedCollectionTag","RobustHelixFinder")),
    _seedflag(pset.get<vector<string> >("HelixFitFlag",vector<string>{"HelixOK"})),
    _minnhits(pset.get<unsigned>("MinNHits",10)),
    _maxdoca(pset.get<double>("MaxDoca",40.0)),
    _foutliers(pset.get<bool>("FilterOutliers",true)),
    _fhoutliers(pset.get<bool>("FilterHelixOutliers",false)),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _perr(pset.get<vector<double> >("ParameterErrors")),
    _seedfit(pset.get<fhicl::ParameterSet>("SeedFit",fhicl::ParameterSet()))
  {
    produces<KalSeedCollection>();
    // check dimensions
    if(_perr.size() != HelixTraj::NHLXPRM)
      throw cet::exception("RECO")<<"mu2e::KalSeedFit: parameter error vector has wrong size"<< endl;
    // mock covariance matrix, all diagonal
    _hcovar = HepSymMatrix(HelixTraj::NHLXPRM,0);
    for(size_t ipar = 0; ipar < HelixTraj::NHLXPRM; ++ipar){
      _hcovar(ipar+1,ipar+1) = _perr[ipar]*_perr[ipar]; // clhep indexing starts a 1
    }
  }

  KalSeedFit::~KalSeedFit(){}

  void KalSeedFit::beginRun(art::Run& ){
 // calculate the helicity
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    GeomHandle<mu2e::TTracker> th;
    _tracker = th.get();

    // change coordinates to mu2e
    CLHEP::Hep3Vector vpoint(0.0,0.0,0.0);
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(vpoint);
    CLHEP::Hep3Vector field = bfmgr->getBField(vpoint_mu2e);
    // helicity is a purely geometric quantity, however it's easiest
    // to determine it from the kinematics (angular momentum and Z momentum)
    _amsign = copysign(1.0,-_tpart.charge()*field.z());
    _helicity = Helicity(static_cast<float>(_fdir.dzdt()*_amsign));
  }

  void KalSeedFit::produce(art::Event& event ) {
    // create output collection
    unique_ptr<KalSeedCollection> kscol(new KalSeedCollection());
    // event printout
    _iev=event.id().event();
    if(_debug > 0 && (_iev%_printfreq)==0)cout<<"KalSeedFit: event="<<_iev<<endl;
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::KalSeedFit: data missing or incomplete"<< endl;
    }

    // loop over the Helices
    for (size_t iseed=0; iseed<_hscol->size(); ++iseed) {
    // convert the HelixSeed to a TrkDef
      HelixSeed const& hseed(_hscol->at(iseed));
      HepVector hpvec(HelixTraj::NHLXPRM);
      // verify the fit meets requirements and can be translated
      // to a fit trajectory.  This accounts for the physical particle direction
      // helicity.  This could be wrong due to FP effects, so don't treat it as an exception
      if(hseed.status().hasAllProperties(_seedflag) &&
	  _helicity == hseed.helix().helicity() &&
	  TrkUtilities::RobustHelix2Traj(hseed._helix,hpvec,_amsign)){
	HelixTraj hstraj(hpvec,_hcovar);
      // update the covariance matrix
	if(_debug > 1)
//	  hstraj.printAll(cout);
	  cout << "Seed Fit HelixTraj parameters " << hstraj.parameters()->parameter()
	  << "and covariance " << hstraj.parameters()->covariance() <<  endl;
// build a time cluster: exclude the outlier hits
	TimeCluster tclust;
	tclust._t0 = hseed._t0;
	for(auto hhit : hseed._hhits){
	  if((!_fhoutliers) || (!hhit._flag.hasAnyProperty(StrawHitFlag::outlier)))	
	    tclust._strawHitIdxs.push_back(hhit._shidx);
	}
// create a TrkDef; it should be possible to build a fit from the helix seed directly FIXME!
	TrkDef seeddef(tclust,hstraj,_tpart,_fdir);
// filter outliers; this doesn't use drift information, just straw positions
	if(_foutliers)filterOutliers(seeddef);

	const HelixTraj* htraj = &seeddef.helix();
	TrkFitFlag       seedok(TrkFitFlag::seedOK);//FIX ME! is there a bettere flag?
	double           flt0  = htraj->zFlight(0.0);
	double           mom   = TrkMomCalculator::vecMom(*htraj, _seedfit.bField(), flt0).mag();
	double           vflt  = seeddef.particle().beta(mom)*CLHEP::c_light;
	double           helt0 = hseed.t0().t0();
	
	KalSeed kf(_tpart,_fdir, hseed.t0(), flt0, seedok);
	auto hsH = event.getValidHandle<HelixSeedCollection>(_hsTag);
	kf._helix = art::Ptr<HelixSeed>(hsH,iseed);
	// extract the hits from the rep and put the hitseeds into the KalSeed
	int nsh = seeddef.strawHitIndices().size();//tclust._strawHitIdxs.size();
	for (int i=0; i< nsh; ++i){
	  size_t          istraw   = seeddef.strawHitIndices().at(i);
	  const StrawHit& strawhit(_shcol->at(istraw));
	  const Straw&    straw    = _tracker->getStraw(strawhit.strawIndex());	  
	  double          fltlen   = htraj->zFlight(straw.getMidPoint().z());
	  double          propTime = (fltlen-flt0)/vflt;

	  //fill the TrkStrwaHitSeed info
	  TrkStrawHitSeed tshs;
	  tshs._index  = istraw;
	  tshs._t0     = TrkT0(helt0 + propTime, hseed.t0().t0Err());
	  tshs._trklen = fltlen; 
	  kf._hits.push_back(tshs);
	}
	
	if(kf._hits.size() >= _minnhits) kf._status.merge(TrkFitFlag::hitsOK);
	// extract the helix trajectory from the fit (there is just 1)
	// use this to create segment.  This will be the only segment in this track
	if(htraj != 0){
	  KalSegment kseg;
	  // sample the momentum at this point
	  BbrVectorErr momerr;// = seedrep->momentumErr(seedrep->flt0());
	  TrkUtilities::fillSegment(*htraj,momerr,kseg);
	  kf._segments.push_back(kseg);
	} else {
	  throw cet::exception("RECO")<<"mu2e::KalSeedFit: Can't extract helix traj from seed fit" << endl;
	}

	// now, fit the seed helix from the filtered hits
	KalRep *seedrep(0);
	_seedfit.makeTrack(_shcol, kf, seedrep);
	
	if(_debug > 1){
	  if(seedrep == 0)
	    cout << "No Seed fit produced " << endl;
	  else
	    cout << "Seed Fit result " << seedrep->fitStatus()  << endl;
	}
	if(seedrep != 0 && (seedrep->fitStatus().success() || _saveall)){
	// convert the status into a FitFlag
	  TrkFitFlag seedok(TrkFitFlag::seedOK);
	  // create a KalSeed object from this fit, recording the particle and fit direction
	  KalSeed kseed(_tpart,_fdir,seedrep->t0(),seedrep->flt0(),seedok);
	  // fill ptr to the helix seed
	  auto hsH = event.getValidHandle<HelixSeedCollection>(_hsTag);
	  kseed._helix = art::Ptr<HelixSeed>(hsH,iseed);
	  // extract the hits from the rep and put the hitseeds into the KalSeed
	  TrkUtilities::fillHitSeeds(seedrep,kseed._hits);
	  if(kseed._hits.size() >= _minnhits)kseed._status.merge(TrkFitFlag::hitsOK);
	  // extract the helix trajectory from the fit (there is just 1)
	  double locflt;
	  const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(seedrep->localTrajectory(seedrep->flt0(),locflt));
	  // use this to create segment.  This will be the only segment in this track
	  if(htraj != 0){
	    KalSegment kseg;
	    // sample the momentum at this point
	    BbrVectorErr momerr = seedrep->momentumErr(seedrep->flt0());
	    TrkUtilities::fillSegment(*htraj,momerr,kseg);
	    kseed._segments.push_back(kseg);
	    // push this seed into the collection
	    kscol->push_back(kseed);
	    if(_debug > 1){
	      cout << "Seed fit segment parameters " << endl;
		for(size_t ipar=0;ipar<5;++ipar) cout << kseg.helix()._pars[ipar] << " ";
	      cout << " covariance " << endl;
	      for(size_t ipar=0;ipar<15;++ipar)
		cout << kseg.covar()._cov[ipar] << " ";
	      cout << endl;
	    }
	  } else {
	    throw cet::exception("RECO")<<"mu2e::KalSeedFit: Can't extract helix traj from seed fit" << endl;
	  }
	}
      // cleanup the seed fit KalRep.  Optimally the seedrep should be a data member of this module
      // and get reused to avoid thrashing memory, but the BTrk code doesn't support that, FIXME!
	delete seedrep;
      }
    }
    // put the tracks into the event
    event.put(move(kscol));
  }


  // find the input data objects 
  bool KalSeedFit::findData(const art::Event& evt){
    _shcol = 0;
    _shfcol = 0;
    _hscol = 0;

    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
    _shfcol = shfH.product();
    auto hsH = evt.getValidHandle<HelixSeedCollection>(_hsTag);
    _hscol = hsH.product();

    return _shcol != 0 && _shfcol != 0 && _hscol != 0;
  }

  void KalSeedFit::filterOutliers(TrkDef& mydef){
  // for now filter on DOCA.  In future this shoudl be an MVA using time and position FIXME!
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint tposp;
    double flt0 = mydef.helix().zFlight(0.0);
    mydef.helix().getInfo(flt0,tposp,tdir);
    // tracker and conditions
    const Tracker& tracker = getTrackerOrThrow();
    const vector<StrawHitIndex>& indices = mydef.strawHitIndices();
    vector<StrawHitIndex> goodhits;
    for(unsigned ihit=0;ihit<indices.size();++ihit){
      StrawHit const& sh = _shcol->at(indices[ihit]);
      Straw const& straw = tracker.getStraw(sh.strawIndex());
      CLHEP::Hep3Vector hpos = straw.getMidPoint();
      CLHEP::Hep3Vector hdir = straw.getDirection();
      // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
      HepPoint spt(hpos.x(),hpos.y(),hpos.z());
      TrkLineTraj htraj(spt,hdir,-straw.getHalfLength(),straw.getHalfLength());
      // estimate flightlength along track.  This assumes a constant BField!!!
      double fltlen = (hpos.z()-tposp.z())/tdir.z();
      HepPoint tp = mydef.helix().position(fltlen);
      Hep3Vector tpos(tp.x(),tp.y(),tp.z()); // ugly conversion FIXME!
      double hitlen = hdir.dot(tpos - hpos);
      TrkPoca hitpoca(mydef.helix(),fltlen,htraj,hitlen);
      // keep hits with small residuals
      if(fabs(hitpoca.doca()) < _maxdoca){
        goodhits.push_back(indices[ihit]);
      }
    }
    // update track
    mydef.strawHitIndices() = goodhits;
  }

}// mu2e
using mu2e::KalSeedFit;
DEFINE_ART_MODULE(KalSeedFit);

///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven track finding
// P.Murat, G.Pezzullo
// try to order routines alphabetically
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/AlgorithmIDCollection.hh"

// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Utilities/make_tool.h"

// data
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/KalSeed.hh"

// conditions
#include "BFieldGeom/inc/BFieldManager.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "DataProducts/inc/Helicity.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"

#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "TrkReco/inc/TrkUtilities.hh"

//BaBar
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/BaBar/BaBar.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkPoca.hh"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/algorithm/string.hpp>

#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "CalPatRec/inc/KalFitHackNew.hh"
#include "CalPatRec/inc/KalFitResultNew.hh"
#include "CalPatRec/inc/CalSeedFit_types.hh"

//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "TMath.h"
#include "TFile.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <float.h>
#include <vector>
#include <set>
#include <map>

using namespace std;
using namespace boost::accumulators;
using CLHEP::HepVector;
using CLHEP::HepSymMatrix;
using CLHEP::Hep3Vector;

namespace mu2e {

  using namespace CalSeedFitTypes;


  class CalSeedFit : public art::EDFilter {
  protected:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
					// configuration parameters
    int                 _diagLevel; 
    int                 _debugLevel;
    int                 _printfreq;
    int                 _useAsFilter;   // 0: producer, 1: filter
    int                 _rescueHits;
//-----------------------------------------------------------------------------
// event object labels
//-----------------------------------------------------------------------------
    std::string         _shLabel ;      // MakeStrawHit label (makeSH)
    std::string         _shDigiLabel;
    std::string         _shpLabel;
    std::string         _shfLabel;
    std::string         _helixSeedLabel;

    double              _maxdtmiss;
					// outlier cuts
    double              _maxAddDoca;
    double              _maxAddChi;
    TrkParticle         _tpart;	        // particle type being searched for
    TrkFitDirection     _fdir;		// fit direction in search
    std::vector<double> _perr;          // diagonal parameter errors to use in the fit

    int                 _nhits_from_gen;//
//-----------------------------------------------------------------------------
// cache of event objects
//-----------------------------------------------------------------------------
    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const StrawHitPositionCollection*     _shpcol;

    const HelixSeedCollection*            _helixSeeds;

    art::Handle<HelixSeedCollection>      _helixSeedsHandle;

    KalFitHackNew                         _fitter;  // Kalman filter config for the Seed fit ( fit using hit wires)

    KalFitResultNew                       _result; // seed fit result

    const TTracker*                       _tracker;     // straw tracker
    const TrackerCalibrations*            _trackerCalib;
    const Calorimeter*                    _calorimeter; // cached pointer to the calorimeter
//-----------------------------------------------------------------------------
// diagnostics 
//-----------------------------------------------------------------------------
    Data_t                                _data;
    std::unique_ptr<ModuleHistToolBase>   _hmanager;

    double                                _amsign;   // cached sign of angular momentum WRT the z axis 
    CLHEP::HepSymMatrix                   _hcovar;   // cache of parameter error covariance matrix
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    enum fitType { helixFit=0, seedFit=1, kalFit=2 };
    explicit CalSeedFit(const fhicl::ParameterSet& PSet);
    virtual ~CalSeedFit();
    
    virtual void beginJob();
    virtual bool beginRun(art::Run&   run  );
    virtual bool filter  (art::Event& event); 
    virtual void endJob  ();
//-----------------------------------------------------------------------------
// helper functions
//-----------------------------------------------------------------------------
    void   findData      (const art::Event& event);
//----------------------------------------------------------------------
// 2015 - 02 - 16 Gianipez
//----------------------------------------------------------------------
    void findMissingHits(KalFitResultNew& KRes);
  };



//-----------------------------------------------------------------------------
// comparison functor for sorting by Z(wire)
//-----------------------------------------------------------------------------
  struct straw_zcomp : public binary_function<StrawHitIndex,StrawHitIndex,bool> {
    bool operator()(StrawHitIndex const& h1, StrawHitIndex const& h2) {

      mu2e::GeomHandle<mu2e::TTracker> handle;
      const TTracker* t = handle.get();
      const Straw* s1 = &t->getStraw(StrawIndex(h1));
      const Straw* s2 = &t->getStraw(StrawIndex(h2));

      return s1->getMidPoint().z() < s2->getMidPoint().z();
    }
  }; // a semicolumn here is required

//-----------------------------------------------------------------------------
// module constructor, parameter defaults are defiend in CalPatRec/fcl/prolog.fcl
//-----------------------------------------------------------------------------
  CalSeedFit::CalSeedFit(fhicl::ParameterSet const& pset) :
    _diagLevel       (pset.get<int>                ("diagLevel"                      )),
    _debugLevel      (pset.get<int>                ("debugLevel"                     )),
    _printfreq       (pset.get<int>                ("printFrequency"                 )),
    _useAsFilter     (pset.get<int>                ("useAsFilter"                    )),    
    _rescueHits      (pset.get<int>                ("rescueHits"                     )),    
    _shLabel         (pset.get<string>             ("StrawHitCollectionLabel"        )),
    _shDigiLabel     (pset.get<string>             ("StrawDigiCollectionLabel"       )),
    _shpLabel        (pset.get<string>             ("StrawHitPositionCollectionLabel")),
    _shfLabel        (pset.get<string>             ("StrawHitFlagCollectionLabel"    )),
    _helixSeedLabel  (pset.get<string>             ("HelixSeedModuleLabel"           )),
    _maxdtmiss       (pset.get<double>             ("MaxDtMiss"                      )),
    _maxAddDoca      (pset.get<double>             ("MaxAddDoca"                     )),
    _maxAddChi       (pset.get<double>             ("MaxAddChi"                      )),
    _tpart           ((TrkParticle::type)(pset.get<int>("fitparticle"                ))),
    _fdir            ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"   ))),
    _perr            (pset.get<vector<double> >    ("ParameterErrors")),
    _fitter          (pset.get<fhicl::ParameterSet>("Fitter",fhicl::ParameterSet())),
    _result        ()
  {
    produces<KalSeedCollection>     ();
    produces<AlgorithmIDCollection> ();

    _hcovar = HepSymMatrix(HelixTraj::NHLXPRM,0);
    for(size_t ipar = 0; ipar < HelixTraj::NHLXPRM; ++ipar){
      _hcovar(ipar+1,ipar+1) = _perr[ipar]*_perr[ipar]; // clhep indexing starts a 1
    }
//-----------------------------------------------------------------------------
// provide for interactive disanostics
//-----------------------------------------------------------------------------
    _data.result    = &_result;

    if (_debugLevel != 0) _printfreq = 1;

    if (_diagLevel != 0) _hmanager = art::make_tool<ModuleHistToolBase>(pset.get<fhicl::ParameterSet>("diagPlugin"));
    else                 _hmanager = std::make_unique<ModuleHistToolBase>();
  }

//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  CalSeedFit::~CalSeedFit() {
  }

//-----------------------------------------------------------------------------
  void CalSeedFit::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    _hmanager->bookHistograms(tfs);
  }

//-----------------------------------------------------------------------------
  bool CalSeedFit::beginRun(art::Run& ) {
    // calculate the helicity
    GeomHandle<BFieldManager>  bfmgr;
    GeomHandle<DetectorSystem> det;
    // change coordinates to mu2e
    CLHEP::Hep3Vector vpoint(0.0,0.0,0.0);
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(vpoint);
    CLHEP::Hep3Vector field       = bfmgr->getBField(vpoint_mu2e);
    // helicity is a purely geometric quantity, however it's easiest
    // to determine it from the kinematics (angular momentum and Z momentum)
    _amsign   = copysign(1.0,-_tpart.charge()*field.z());
    
    mu2e::GeomHandle<mu2e::TTracker> th;
    _tracker = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
    // calibrations

    mu2e::ConditionsHandle<TrackerCalibrations> tcal("ignored");
    _trackerCalib = tcal.operator ->();

    _fitter.setTracker(_tracker);
    _fitter.setTrackerCalib(_trackerCalib);
    _fitter.setCalorimeter(_calorimeter);
    
    // ConditionsHandle<AcceleratorParams> accPar("ignored");
    // _data.mbtime = accPar->deBuncherPeriod;

    return true;
  }

//-----------------------------------------------------------------------------
// find the input data objects
//-----------------------------------------------------------------------------
  void CalSeedFit::findData(const art::Event& evt) {

    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if (evt.getByLabel(_shLabel,strawhitsH)) {
      _shcol = strawhitsH.product();
    }
    else {
      _shcol  = 0;
      printf(" >>> ERROR in CalSeedFit::findData: StrawHitCollection with label=%s not found. BAIL OUT\n",
             _shLabel.data());
    }

    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if (evt.getByLabel(_shpLabel,shposH)) {
      _shpcol = shposH.product();
    }
    else {
      _shpcol = 0;
      printf(" >>> ERROR in CalSeedFit::findData: StrawHitPositionCollection with label=%s not found. BAIL OUT\n",
             _shpLabel.data());
    }

    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if (evt.getByLabel(_shfLabel,shflagH)) {
      _shfcol = shflagH.product();
    }
    else {
      _shfcol = 0;
      printf(" >>> ERROR in CalSeedFit::findData: StrawHitFlagCollection with label=%s not found. BAIL OUT\n",
             _shfLabel.data());
    }

    if (evt.getByLabel(_helixSeedLabel, _helixSeedsHandle)) {
      _helixSeeds = _helixSeedsHandle.product();
    }
    else {
      printf(" >>> ERROR in CalSeedFit::findData: HelixSeedCollection with label=%s not found. BAIL OUT\n",
             _shpLabel.data());
      _helixSeeds = 0;
    }
//-----------------------------------------------------------------------------
// theow an exception if data not found
//-----------------------------------------------------------------------------
    if ((_shcol == 0) || (_shfcol == 0) || (_shpcol == 0) || (_helixSeeds == 0)) {
      throw cet::exception("RECO")<<"mu2e::KalSeedFit: data missing or incomplete"<< endl;
    }
    
  }

//-----------------------------------------------------------------------------
// event entry point
//-----------------------------------------------------------------------------
  bool CalSeedFit::filter(art::Event& event ) {
    const char*       oname = "CalSeedFit::filter";
    char              message[200];

    const HelixSeed*  helixSeed(0);
    const HelixHit*   hhit(0);

    int iev = event.event();
    if ((iev%_printfreq) == 0) printf("[%s] : START event number %8i\n",oname,iev);

    unique_ptr<KalSeedCollection>      tracks   (new KalSeedCollection     );
    unique_ptr<AlgorithmIDCollection>  algs     (new AlgorithmIDCollection);

    findData(event);
					// initialize the diagnostic structure
    _data.event  = &event;
    _data.result = &_result;
    _data.nrescued.clear();
    _data.mom.clear();

    _data.tracks        = tracks.get();

    _result.fitType     = 0;
    _result.event       = &event ;
    _result.shcol       = _shcol ;
    _result.shpos       = _shpcol;
    _result.shfcol      = _shfcol;
    _result.tpart       = _tpart ;
    _result.fdir        = _fdir  ;
    _result.shDigiLabel = _shDigiLabel;
//-----------------------------------------------------------------------------
// loop over found "time peaks" -  calorimeter clusters above the energy threshold
//-----------------------------------------------------------------------------
    int npeaks = _helixSeeds->size();

    for (int ipeak=0; ipeak<npeaks; ipeak++) {
      helixSeed = &_helixSeeds->at(ipeak);
      int n_helix_hits = helixSeed->hits().size();

      if (n_helix_hits < _fitter.minNStraws())  continue;


      _result.helixSeed   = helixSeed;
      _result.caloCluster = helixSeed->caloCluster().get();
      _result.t0          = helixSeed->t0();

      HepVector hpvec(HelixTraj::NHLXPRM);
      TrkUtilities::RobustHelix2Traj(helixSeed->helix(), hpvec, _amsign);
      if (_result.helixTraj == NULL) _result.helixTraj = new HelixTraj(hpvec,_hcovar); 
      else                          *_result.helixTraj =     HelixTraj(hpvec,_hcovar); 
//-----------------------------------------------------------------------------
// start from hits of the helix candidate
//-----------------------------------------------------------------------------
      _result.hitIndices->clear();

      for (int i=0; i<n_helix_hits; ++i){
	hhit = &helixSeed->hits().at(i);
	if (!hhit->_flag.hasAnyProperty(StrawHitFlag::outlier)) {
	  _result.hitIndices->push_back(hhit->index());
	}
      }
//-----------------------------------------------------------------------------
// seed fit - fit through the wires of found hits, not using the drift times
// fitting through the wires is achieved by using TrkReco/src/FixedAmbigResolver
// which always sets the driftAmbiguity to zero
//-----------------------------------------------------------------------------
      _fitter.makeTrack(_result);
//--------------------------------------------------------------------------------
// 2014-11-24 gianipez added the following diagnostic
//--------------------------------------------------------------------------------
      if (_debugLevel > 0) {
	sprintf(message,
		"CalSeedFit::produce after seedfit::makeTrack: fit_success = %i\n",
		_result.fit.success());
	_fitter.printHits(_result,message);
      }

      if (_result.fit.success()) {
//-----------------------------------------------------------------------------
// track is successfully found, try to pick up hits missed by the helix finder
// at this step, ignore the calorimeter cluster when refitting the track
//-----------------------------------------------------------------------------
	int nrescued = 0;
	if (_rescueHits) { 
	  findMissingHits(_result);
	  nrescued = _result.missingHits.size();
	  if (nrescued > 0) {
	    _result.caloCluster = NULL;
	    _fitter.addHits(_result,_maxAddChi);
	  }
	}
//-----------------------------------------------------------------------------
// final printout
//-----------------------------------------------------------------------------
	if (_debugLevel > 0) {
	  sprintf(message,
		  "CalSeedFit::produce after seedfit::addHits: fit_success = %i\n",
		  _result.fit.success());
	  _fitter.printHits(_result,message);
	}
//-----------------------------------------------------------------------------
// form the output
//-----------------------------------------------------------------------------
	KalRep* krep = _result.stealTrack();

					// convert the status into a FitFlag
	TrkFitFlag seedok(TrkFitFlag::seedOK);
	KalSeed    kseed (_tpart, _fdir, krep->t0(), krep->flt0(), seedok);
	// fill ptr to the helix seed
	kseed._helix = art::Ptr<HelixSeed>(_helixSeedsHandle, ipeak);
	// extract the hits from the rep and put the hitseeds into the KalSeed
	TrkUtilities::fillHitSeeds(krep, kseed._hits);
	kseed._status.merge(TrkFitFlag::hitsOK);

					// extract the helix trajectory from the fit (there is just 1)
	double locflt;
	const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(krep->localTrajectory(krep->flt0(),locflt));
//-----------------------------------------------------------------------------
// use this to create segment.  This will be the only segment in this track
//-----------------------------------------------------------------------------
	if (htraj != 0) {
	  KalSegment kseg;
	  // sample the momentum at this point
	  BbrVectorErr momerr = krep->momentumErr(krep->flt0());
	  TrkUtilities::fillSegment(*htraj, momerr, kseg);
	  kseed._segments.push_back(kseg);

	  kseed._chisq = krep->chisq();
	  // use the default consistency calculation, as t0 is not fit here
	  kseed._fitcon = krep->chisqConsistency().significanceLevel();

	  // push this seed into the collection
	  tracks->push_back(kseed);
	  
	  if (_debugLevel > 1) {
	    printf("Seed fit segment parameters \n");
	    for(size_t ipar=0;ipar<5;++ipar) cout << kseg.helix()._pars[ipar] << " ";
	    printf(" covariance \n");
	    for (size_t ipar=0;ipar<15;++ipar) {
	      cout << kseg.covar()._cov[ipar] << " ";
	    }
	    cout << endl;
	  }

	  if (_diagLevel > 0) {
//-----------------------------------------------------------------------------
// store some info for convenient diagnostics
//-----------------------------------------------------------------------------
	    double  h1_fltlen      = krep->firstHit()->kalHit()->hit()->fltLen();
	    double  hn_fltlen      = krep->lastHit ()->kalHit()->hit()->fltLen();
	    double  entlen         = std::min(h1_fltlen, hn_fltlen);

	    CLHEP::Hep3Vector fitmom = krep->momentum(entlen);
	    _data.mom.push_back(fitmom.mag());
 	    _data.nrescued.push_back(nrescued);
	  }
	}
	else {
	  throw cet::exception("RECO")<<"mu2e::KalSeedFit: Can't extract helix traj from seed fit" << endl;
	}
	
	int best = AlgorithmID::CalPatRecBit;
	int mask = 1 << AlgorithmID::CalPatRecBit;
	
	algs->push_back(AlgorithmID(best,mask));

	delete krep;
      }
      else {
	_result.deleteTrack();
      }
    }

    if (_diagLevel > 0) _hmanager->fillHistograms(&_data);
//-----------------------------------------------------------------------------
// put reconstructed tracks into the event record and do filtering
//-----------------------------------------------------------------------------
    event.put(std::move(tracks));
    event.put(std::move(algs  ));

    if (_useAsFilter == 0) return true;
    else                   return (tracks->size() > 0);
    
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalSeedFit::endJob(){
    // does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }

//-----------------------------------------------------------------------------
// moved from CalSeedFit - for the moment, comment out the source 
// look for hits which were not a part of the helix hit list around the 
// trajectory found by the seed fit
// look at all hits included into the corresponding time cluster
// first reactivate already associated hits
//-----------------------------------------------------------------------------
  void CalSeedFit::findMissingHits(KalFitResultNew& KRes) {

    const char* oname = "CalSeedFit::findMissingHits";

    mu2e::TrkStrawHit*       hit;
    int                      hit_index;
    const StrawHit*          sh;
    const Straw*             straw;

    Hep3Vector               tdir;
    HepPoint                 tpos;
    double                   doca, /*rdrift, */fltlen;

    if (_debugLevel > 0) printf("[%s]: BEGIN\n",oname);

    const KalRep* krep = KRes.krep;

    KRes.missingHits.clear();
    //    KRes.doca.clear();

    const TrkDifTraj& trajectory = krep->traj();
    const vector<TrkHit*>&  trackHits  = krep->hitVector();
//-----------------------------------------------------------------------------
// get track position and direction at S=0
//-----------------------------------------------------------------------------
    trajectory.getInfo(0.0,tpos,tdir);
//-----------------------------------------------------------------------------
// look for so far unused hits around the trajectory
//-----------------------------------------------------------------------------
    const HelixSeed*   hseed = KRes.helixSeed;
    const  std::vector<StrawHitIndex>& tchits = hseed->timeCluster()->hits();
 
    int n = tchits.size();
    for (int i=0; i<n; ++i) {
      hit_index = tchits.at(i);
      sh        = &KRes.shcol->at(hit_index);
      straw     = &_tracker->getStraw(sh->strawId());

      const CLHEP::Hep3Vector& wpos = straw->getMidPoint();
      const CLHEP::Hep3Vector& wdir = straw->getDirection();
	    
      HepPoint      wpt  (wpos.x(),wpos.y(),wpos.z());
      TrkLineTraj   wire (wpt,wdir,-20,20);
//-----------------------------------------------------------------------------
// estimate flightlength along the track for z-coordinate corresponding to the 
// wire position. This assumes a constant BField!!!
// in principle, this should work well enough, however, may want to check
// then determine the distance from the wire to the trajectory
//-----------------------------------------------------------------------------
      fltlen = (wpos.z()-tpos.z())/tdir.z();
      TrkPoca wpoca(trajectory,fltlen,wire,0.0);
      doca   = wpoca.doca();

      int found(-1);

      if (std::fabs(doca) < _maxAddDoca) {
	found = 0;
	for (auto it=trackHits.begin(); it<trackHits.end(); it++) {
	  hit    = static_cast<mu2e::TrkStrawHit*> (*it);
	  int shIndex = int(hit->index());
	  if (hit_index == shIndex) {
	    found = 1;
	    break;
	  }
	}
//-----------------------------------------------------------------------------
// CalSeedFit doesn't look at the hit residuals, only wires
//-----------------------------------------------------------------------------
	if (found == 0) {
	  MissingHit_t mh;
	  mh.index = hit_index;
	  mh.doca  = doca;
	  mh.dr    = doca;
	  KRes.missingHits.push_back(mh);
	  //	  KRes.doca.push_back(doca);
	}
      }

      if (_debugLevel > 0) printf("[%s] %5i %8.3f %2i \n",oname,hit_index,doca,found);

    }
  }

}

using mu2e::CalSeedFit;
DEFINE_ART_MODULE(CalSeedFit);

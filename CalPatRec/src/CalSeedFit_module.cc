///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven track finding
// P.Murat, G.Pezzullo
// try to order routines alphabetically
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/CalSeedFit_module.hh"
#include "CalPatRec/inc/AlgorithmIDCollection.hh"

// framework
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Utilities/make_tool.h"

// conditions
#include "BFieldGeom/inc/BFieldManager.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

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

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/algorithm/string.hpp>

#include "TVector2.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"

using namespace std;
using namespace boost::accumulators;
using CLHEP::HepVector;
using CLHEP::HepSymMatrix;
using CLHEP::Hep3Vector;

namespace mu2e {
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

    if (_diagLevel != 0) _hmanager = art::make_tool<CprModuleHistBase>(pset.get<fhicl::ParameterSet>("histograms"));
    else                 _hmanager = std::make_unique<CprModuleHistBase>();
  }

//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  CalSeedFit::~CalSeedFit() {
  }

//-----------------------------------------------------------------------------
  void CalSeedFit::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    _hmanager->bookHistograms(tfs,&_hist);
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
  bool CalSeedFit::findData(const art::Event& evt) {

    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if (evt.getByLabel(_shLabel,strawhitsH)) {
      _shcol = strawhitsH.product();
      _nhits = _shcol->size();
    }
    else {
      _shcol  = 0;
      _nhits  = 0;
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
//-----------------------------------------------------------------------------
// find list of MC hits - for debugging only
//-----------------------------------------------------------------------------
    art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandle;
    evt.getByLabel(_shDigiLabel,mcptrHandle);
    if (mcptrHandle.isValid()) {
      _listOfMCStrawHits = (mu2e::PtrStepPointMCVectorCollection*) mcptrHandle.product();
    }
    else {
      _listOfMCStrawHits = NULL;
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
// done
//-----------------------------------------------------------------------------
    return (_shcol != 0) && (_shfcol != 0) && (_shpcol != 0) && (_helixSeeds != 0);
  }

//-----------------------------------------------------------------------------
// event entry point
//-----------------------------------------------------------------------------
  bool CalSeedFit::filter(art::Event& event ) {
    const char*    oname = "CalSeedFit::filter";
    char           message[200];
    //    bool           findseed(false);
    int            npeaks;

    ::KalRep*      krep;

    const HelixSeed*   helixSeed(0);
    const HelixHit*    hhit(0);

    _data.event   = &event;
    _data.result  = &_result;
    _data.ntracks = 0;
    _data.nrescued.clear();

    _eventid = event.event();
    _iev     = event.id().event();

    if ((_iev%_printfreq) == 0) printf("[%s] : START event number %8i\n", oname,_iev);

    unique_ptr<KalSeedCollection>      tracks   (new KalSeedCollection     );
    unique_ptr<AlgorithmIDCollection>  algs     (new AlgorithmIDCollection);

    if (!findData(event))  goto END;

    _fitter.setStepPointMCVectorCollection(_listOfMCStrawHits);

    _result._fitType     = 0;
    _result._event       = &event ;
    _result._shcol       = _shcol ;
    _result._shpos       = _shpcol;
    _result._shfcol      = _shfcol;
    _result._tpart       = _tpart ;
    _result._fdir        = _fdir  ;
    _result._shDigiLabel = _shDigiLabel;
//-----------------------------------------------------------------------------
// loop over found "time peaks" -  calorimeter clusters above the energy threshold
//-----------------------------------------------------------------------------
    npeaks = _helixSeeds->size();
    
    for (int ipeak=0; ipeak<npeaks; ipeak++) {
      helixSeed = &_helixSeeds->at(ipeak);
      int n_helix_hits = helixSeed->hits().size();

      if (n_helix_hits < _fitter.minNStraws())  continue;


      _result._helixSeed   = helixSeed;
      _result._caloCluster = helixSeed->caloCluster().get();
      _result._t0          = helixSeed->t0();
      HepVector hpvec(HelixTraj::NHLXPRM);
      TrkUtilities::RobustHelix2Traj(helixSeed->helix(), hpvec, _amsign);
      if (_result._helixTraj == NULL) _result._helixTraj = new HelixTraj(hpvec,_hcovar); 
      else                           *_result._helixTraj =     HelixTraj(hpvec,_hcovar); 
//-----------------------------------------------------------------------------
// start from hits of the helix candidate
//-----------------------------------------------------------------------------
      _result._hitIndices->clear();

      for (int i=0; i<n_helix_hits; ++i){
	hhit = &helixSeed->hits().at(i);
	if (!hhit->_flag.hasAnyProperty(StrawHitFlag::outlier)) {
	  _result._hitIndices->push_back(hhit->index());
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
		_result._fit.success());
	_fitter.printHits(_result,message);
      }

      if (_result._fit.success()) {
	_data.ntracks += 1;
//-----------------------------------------------------------------------------
// track is successfully found, try to pick up additional hits missed 
// by the helix finder
// at this step, ignore presence of the calorimeter cluster when refitting the track
//-----------------------------------------------------------------------------
	int nrescued = 0;
	if (_rescueHits) { 
	  findMissingHits(_result);
	  nrescued = _result._missingHits.size();
	  if (nrescued > 0) {
	    _result._caloCluster = NULL;
	    _fitter.addHits(_result,_maxAddChi);
	  }
	}
	_data.nrescued.push_back(nrescued);
//-----------------------------------------------------------------------------
// final printout
//-----------------------------------------------------------------------------
	if (_debugLevel > 0) {
	  sprintf(message,
		  "CalSeedFit::produce after seedfit::addHits: fit_success = %i\n",
		  _result._fit.success());
	  _fitter.printHits(_result,message);
	}
//-----------------------------------------------------------------------------
// track fit completed, per-track histogramming
//-----------------------------------------------------------------------------
	if (_diagLevel > 0) {
	  _hmanager->fillHistograms(1,&_data,&_hist);
	}
//-----------------------------------------------------------------------------
// form the output
//-----------------------------------------------------------------------------
	krep      = _result.stealTrack();
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

    if (_diagLevel > 0) {
      _data.ntracks = tracks->size();
      _hmanager->fillHistograms(0,&_data,&_hist);
    }
//-----------------------------------------------------------------------------
// put reconstructed tracks into the event record and do filtering
//-----------------------------------------------------------------------------
  END:;

    event.put(std::move(tracks));
    event.put(std::move(algs  ));

    if (_useAsFilter == 0) return true;
    else                   return (_data.ntracks > 0);

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
    for (unsigned istr=0; istr<nstrs;++istr) {
//----------------------------------------------------------------------
// 2015-02-11 gianipez and P. Murat changed the selection bit
//            for searching for missed hits
//----------------------------------------------------------------------
      StrawHit const& sh = _shcol->at(istr);
//-----------------------------------------------------------------------------
// I think, we want to check the radial bit: if it is set, than at least one of
// the two measured times is wrong...
//-----------------------------------------------------------------------------
      radius_ok = _shfcol->at(istr).hasAllProperties(StrawHitFlag::radsel);
      dt        = _shcol->at(istr).time()-kalfit._krep->t0()._t0;

      if (radius_ok && (fabs(dt) < _maxdtmiss)) {
        // make sure we haven't already used this hit
	TrkStrawHitVector tshv;
	convert(kalfit._hits, tshv);
        TrkStrawHitVector::iterator ifnd = find_if(tshv.begin(), tshv.end(),FindTrkStrawHit(sh));
        if(ifnd == tshv.end()){
          // good in-time hit.  Compute DOCA of the wire to the trajectory
          Straw const& straw = _tracker->getStraw(sh.strawIndex());
          CLHEP::Hep3Vector hpos = straw.getMidPoint();
          CLHEP::Hep3Vector hdir = straw.getDirection();
          // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
          HepPoint spt(hpos.x(),hpos.y(),hpos.z());
          TrkLineTraj htraj(spt,hdir,-20,20);
          // estimate flightlength along track.  This assumes a constant BField!!!
          double fltlen = (hpos.z()-tpos.z())/tdir.z();
          TrkPoca hitpoca(kalfit._krep->pieceTraj(),fltlen,htraj,0.0);

          if (_debugLevel > 0) {
            printf("[CalSeedFit::findMissingHits] %8i  %6i  %8i  %10.3f \n",
                   straw.index().asInt(),
                   straw.id().getPlane(),
                   straw.id().getPanel(),
                   hitpoca.doca());
          }
//-----------------------------------------------------------------------------
// flag hits with small residuals
//-----------------------------------------------------------------------------
          if (fabs(hitpoca.doca()) < _maxadddoca) {
            misshits.push_back(istr);
          }
        }
      }
    }

    mu2e::TrkStrawHit*       hit;
    int                      hit_index;
    const StrawHit*          sh;
    const Straw*             straw;

    Hep3Vector               tdir;
    HepPoint                 tpos;
    double                   doca, /*rdrift, */fltlen;

    if (_debugLevel > 0) printf("[%s]: BEGIN\n",oname);

    const KalRep* krep = KRes._krep;

    KRes._missingHits.clear();
    KRes._doca.clear();

    const TrkDifTraj& trajectory = krep->traj();
    const vector<TrkHit*>&  trackHits  = krep->hitVector();
//-----------------------------------------------------------------------------
// get track position and direction at S=0
//-----------------------------------------------------------------------------
    trajectory.getInfo(0.0,tpos,tdir);
//-----------------------------------------------------------------------------
// look for so far unused hits around the trajectory
//-----------------------------------------------------------------------------
    const HelixSeed*   hseed = KRes._helixSeed;
    const  std::vector<StrawHitIndex>& tchits = hseed->timeCluster()->hits();
 
    int n = tchits.size();
    for (int i=0; i<n; ++i) {
      hit_index = tchits.at(i);
      sh        = &KRes._shcol->at(hit_index);
      straw     = &_tracker->getStraw(sh->strawIndex());

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
	  KRes._missingHits.push_back(hit_index);
	  KRes._doca.push_back(doca);
	}
      }

      if (_debugLevel > 0) printf("[%s] %5i %8.3f %2i \n",oname,hit_index,doca,found);

    }
  }

}

using mu2e::CalSeedFit;
DEFINE_ART_MODULE(CalSeedFit);

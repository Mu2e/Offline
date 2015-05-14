///////////////////////////////////////////////////////////////////////////////
// $Id: KalFitHack.cc,v 1.4 2014/04/08 04:25:46 murat Exp $
// $Author: murat $ 
// $Date: 2014/04/08 04:25:46 $
///////////////////////////////////////////////////////////////////////////////
// framework
#include "fhiclcpp/ParameterSet.h"
// the following has to come before other BaBar includes
#include "BaBar/BaBar.hh"
#include "CalPatRec/inc/KalFitHack.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/PanelAmbigResolver.hh"
#include "KalmanTests/inc/PocaAmbigResolver.hh"
#include "KalmanTests/inc/HitAmbigResolver.hh"
#include "KalmanTests/inc/FixedAmbigResolver.hh"
#include "KalmanTests/inc/DoubletAmbigResolver.hh"
#include "KalmanTests/inc/BaBarMu2eField.hh"
//geometry
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "KalmanTrack/KalHit.hh"
#include "TrkBase/HelixTraj.hh"
//09 - 26 - 2013 gianipez added the following include file
#include "TrkBase/HelixParams.hh"
//----------------------------------------
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/TrkHotListFull.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "TrkBase/TrkPoca.hh"
#include "BaBar/ErrLog.hh"
#include "BField/BFieldFixed.hh"
#include "DetectorModel/DetIntersection.hh"
#include "DetectorModel/DetMaterial.hh"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>

using namespace std; 

namespace mu2e 
{
// comparison functor for ordering hits
  struct fltlencomp : public binary_function<TrkStrawHit*, TrkStrawHit*, bool> {
    fltlencomp(TrkFitDirection::FitDirection fdir=TrkFitDirection::downstream) : _fdir(fdir) {}
    bool operator()(TrkStrawHit* x, TrkStrawHit* y) { 
      return _fdir == TrkFitDirection::downstream ? x->fltLen() < y->fltLen() : y->fltLen() < x->fltLen() ;
    }
    TrkFitDirection::FitDirection _fdir;
  };

  struct timecomp : public binary_function<TrkStrawHit*, TrkStrawHit*, bool> {
    timecomp() {}
    bool operator()(TrkStrawHit* x, TrkStrawHit* y) { 
      return x->hitT0()._t0 < y->hitT0()._t0;
    }
    TrkFitDirection::FitDirection _fdir;
  };
// construct from a parameter set  
  KalFitHack::KalFitHack(fhicl::ParameterSet const& pset) :
// KalFitHack parameters
    _debug(pset.get<int>("debugLevel",0)),
    _weedhits(pset.get<bool>("weedhits")),
    _maxhitchi(pset.get<double>("maxhitchi",4.0)),
    _maxweed(pset.get<unsigned>("maxweed",10)),
    _hiterr(pset.get< vector<double> >("hiterr")),
    _maxdriftpull(pset.get<double>("maxDriftPull",10)),
    // t0 parameters
    _initt0(pset.get<bool>("initT0",true)),
    _updateT0(pset.get<bool>("updateT0",true)),
    _updateT0Mode(pset.get<int>("T0UpdateMode")),
    fMinHitDrift(pset.get<double>("HitMinDrift")),
    fRdriftMinusDocaTol(pset.get<double>("RdriftMinusDocaTol")),
    _daveMode(pset.get<int>("daveMode" ,0)),
    _t0tol(pset.get< vector<double> >("t0Tolerance")),
    _t0errfac(pset.get<double>("t0ErrorFactor",1.2)),
    _mint0doca(pset.get<double>("minT0DOCA",-0.2)),
    _t0nsig(pset.get<double>("t0window",2.5)),
    _dtoffset(pset.get<double>("dtOffset")),
    fScaleErrDoublet(pset.get<double>("scaleErrDoublet")),
    //    fUseDoublets(0),
    fMinDriftDoublet  (pset.get<double>("minDriftDoublet")),
    fDeltaDriftDoublet(pset.get<double>("deltaDriftDoublet")),
    _maxDoubletChi2   (pset.get<double>("maxDoubletChi2",9.)),
    fSigmaSlope       (pset.get<double>("sigmaSlope")),
    fMakeStrawHitModuleLabel(pset.get<std::string>("makeStrawHitModuleLabel")),
    //
    _removefailed(pset.get<bool>("RemoveFailedFits",true)),
    _minnstraws(pset.get<unsigned>("minnstraws",15)),
    _ambigstrategy(pset.get< vector<int> >("ambiguityStrategy")),
    _bfield(0),
    _nIter(0)
  {
    //    fStopwatch = new TStopwatch();
    _darPset = new fhicl::ParameterSet(pset.get<fhicl::ParameterSet>("DoubletAmbigResolver",fhicl::ParameterSet()));

// set KalContext parameters
    _disttol = pset.get<double>("IterationTolerance",0.1);
    _intertol = pset.get<double>("IntersectionTolerance",100.0);
    _maxiter = pset.get<long>("MaxIterations",10);
    _maxinter = pset.get<long>("MaxIntersections",0);
    _matcorr = pset.get<bool>("materialCorrection",true);
    _fieldcorr = pset.get<bool>("fieldCorrection",false);
    _smearfactor = pset.get<double>("SeedSmear",1.0e6);
    _sitethresh = pset.get<double>("SiteMomThreshold",0.2);
    _momthresh = pset.get<double>("MomThreshold",10.0);
    _mingap = pset.get<double>("mingap",0.1);
    _minfltlen = pset.get<double>("MinFltLen",0.1);
    _minmom = pset.get<double>("MinMom",10.0);
    _fltepsilon = pset.get<double>("FltEpsilon",0.001);
    _divergeflt = pset.get<double>("DivergeFlt");
    _mindot = pset.get<double>("MinDot",0.0);
    _maxmomdiff = pset.get<double>("MaxMomDiff",0.5);
    _momfac = pset.get<double>("MomFactor",0.0);
    _maxpardif[0] = _maxpardif[1] = pset.get<double>("MaxParameterDifference",1.0);
    // DOF counting subdivision is illogical, FIXME!!!!
    _mindof[0] = _mindof[2] = pset.get<double>("MinNDOF",10);
    _mindof[1] = 0;
//----------------------------------------------------------------------
// 2015-01-09 G.Pezzullo and P.Murat
// Noticed that with respect to KalmanTest/src/KalFit.cc we were using
// different ranges and divisions for the magnetic field. *fixed*
//----------------------------------------------------------------------
    _bintconfig._maxRange       = pset.get<double>("BFieldIntMaxRange",  1.0e5); // 100 m
    _bintconfig._intTolerance   = pset.get<double>("BFieldIntTol"     ,  0.01 ); // 10 KeV
    _bintconfig._intPathMin     = pset.get<double>("BFieldIntMin"     , 20.0  ); // 20 mm
    _bintconfig._divTolerance   = pset.get<double>("BFieldIntDivTol"  ,  0.05 ); // 50 KeV
    _bintconfig._divPathMin     = pset.get<double>("BFieldIntDivMin"  , 50.0  ); // 50 mm
    _bintconfig._divStepCeiling = pset.get<double>("BFieldIntDivMax"  ,500.0  ); // 500 mm
//-----------------------------------------------------------------------------
// initialize sequence of drift signs
//-----------------------------------------------------------------------------
    double s[4][2] = { 1, 1, 1, -1, -1, -1, -1, 1} ;
    for (int i=0; i<4; i++) {
      for (int j=0; j<2; j++) {
	fSign[i][j] = s[i][j];
      }
    }
//-----------------------------------------------------------------------------
// make sure we have at least one entry for additional errors
//-----------------------------------------------------------------------------
    if (_hiterr.size() <= 0) {
      throw cet::exception("RECO") 
	<< "mu2e::KalFitHack: no hit errors specified" << endl;
    }

    if(_hiterr.size() != _ambigstrategy.size()) {
      throw cet::exception("RECO") 
	<< "mu2e::KalFitHack: inconsistent ambiguity resolution" << endl;
    }

    if(_hiterr.size() != _t0tol.size()) {
      throw cet::exception("RECO") 
	<<"mu2e::KalFitHack: inconsistent ambiguity resolution" << endl;
    }
					// construct the ambiguity resolvers
    AmbigResolver* ar;
    double         err;
    int n = _ambigstrategy.size();
    for(int i=0; i<n; ++i) {
      err = _hiterr[i];
      switch (_ambigstrategy[i]) {
      case kFixedAmbig: default:
	ar = new FixedAmbigResolver(pset,err,i);
	break;
      case kHitAmbig:
	ar = new HitAmbigResolver(pset,err,i);
	break;
      case kPanelAmbig:
	ar = new PanelAmbigResolver(pset,err,i);
	break;
      case kPocaAmbig:
	ar = new PocaAmbigResolver(pset,err,i);
	break;
      case kDoubletAmbig: // 4
	ar = new DoubletAmbigResolver(*_darPset,err,i);
	break;
      }
      _ambigresolver.push_back(ar);
    }
  }

//-----------------------------------------------------------------------------
  KalFitHack::~KalFitHack(){
    for(size_t i=0; i<_ambigresolver.size(); ++i){
      delete _ambigresolver[i];
    }

    delete _bfield;
    delete _darPset;

    //    delete fStopwatch;
  }

//------------------------------------------------------------------------------------------
// called once per event from CalPatRec_module::produce
//-----------------------------------------------------------------------------
  void KalFitHack::makeTrack(KalFitResult& kres, CalTimePeak* TPeak) {

    kres._fit = TrkErrCode(TrkErrCode::fail);

					// test if fitable
    if (! fitable(*kres._tdef)) return;
//-----------------------------------------------------------------------------
// first, find t0
//-----------------------------------------------------------------------------
    TrkT0 t0;
    bool caloInitCond(false);
    if (TPeak->Cluster() != NULL) {
      caloInitCond = true;
    }

    if (_initt0) {
      if (!caloInitCond) {
	initT0(*kres._tdef, t0);
      } 
      else {
	initCaloT0(TPeak, *kres._tdef, t0);
      }
    }
    else {
      t0 = kres._tdef->t0();
    }
//-----------------------------------------------------------------------------
// knowing t0, create the hits
//-----------------------------------------------------------------------------
    makeHits(kres, t0);
//-----------------------------------------------------------------------------
// Create the BaBar hit list, and fill it with these hits.  The BaBar list takes ownership
// This will go away when we cleanup the BaBar hit storage, FIXME!!!
//-----------------------------------------------------------------------------
    TrkHotListFull* hotlist = new TrkHotListFull();
    for(std::vector<TrkStrawHit*>::iterator ihit=kres._hits.begin();ihit!=kres._hits.end();ihit++){
      TrkStrawHit* trkhit = *ihit;
      hotlist->append(trkhit);
    }
//-----------------------------------------------------------------------------
// Find the wall and gas material description objects for these hits
//-----------------------------------------------------------------------------
    if (_matcorr) makeMaterials(kres);
// create Kalman rep
    kres._krep = new KalRep(kres._tdef->helix(), hotlist, kres._detinter, *this, kres._tdef->particle());
    assert(kres._krep != 0);
// initialize krep t0; eventually, this should be in the constructor, FIXME!!!
    double flt0 = kres._tdef->helix().zFlight(0.0);
    kres._krep->setT0(t0,flt0);
    
    if (_debug>0){
      printHits(kres,"makeTrack_001");
    }
//-----------------------------------------------------------------------------
// now fit
// 10-07-2013 giani added the following line. It updates the hit times
//            following changes in the t0 value 
//-----------------------------------------------------------------------------
    if(caloInitCond) updateHitTimes(kres);
//-----------------------------------------------------------------------------
// 09 - 26 - 2013 giani 
// include the calorimeter information when it is avaiable
//-----------------------------------------------------------------------------
    if ((_daveMode == 0) && caloInitCond) {
      fitTrack(kres,TPeak);
    } 
    else {
      fitTrack(kres,NULL);
    }
    
    if (_removefailed) kres.removeFailed();
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void KalFitHack::addHits(KalFitResult&              kres   , 
			   const StrawHitCollection*  straws , 
			   std::vector<hitIndex>      indices, 
			   double                     maxchi ,
			   CalTimePeak*               TPeak   ) {

					// there must be a valid Kalman fit to add hits to
    int  activity(0);
    int  final(1);

    if(kres._krep != 0 && kres._fit.success()){
      ConditionsHandle<TrackerCalibrations> tcal("ignored");
      const Tracker& tracker = getTrackerOrThrow();
      std::vector<TrkStrawHit*>::iterator ihigh;
      std::vector<TrkStrawHit*>::reverse_iterator ilow;
//-----------------------------------------------------------------------------
// use the reference trajectory, as that's what all the existing hits do
//-----------------------------------------------------------------------------
      const TrkDifPieceTraj* reftraj = kres._krep->referenceTraj();

      if (_debug>0){
	printf("[KalFitHack::addHits]  shId   A      sec      panel      res        hitRMS       drift      chi2  \n");
      }

      double             hflt, mom, beta, tflt;
      const TrkStrawHit* nearhit;
      TrkStrawHit*       trkhit;
	

      for (unsigned iind=0;iind<indices.size(); ++iind) {
	size_t istraw = indices[iind]._index;
	const StrawHit& strawhit(straws->at(istraw));
	const Straw& straw = tracker.getStraw(strawhit.strawIndex());
					// estimate  initial flightlength
	hflt = 0;
	TrkHelixUtils::findZFltlen(*reftraj,straw.getMidPoint().z(),hflt);
// find the bounding sites near this hit, and extrapolate to get the hit t0
	std::sort(kres._hits.begin(),kres._hits.end(),fltlencomp(kres._tdef->fitdir().fitDirection()));
	findBoundingHits(kres._hits,hflt,ilow,ihigh);

	if(ihigh != kres._hits.end()) nearhit = *ihigh;
	else                          nearhit = *ilow;

	TrkT0 hitt0 = nearhit->hitT0();

	mom  = kres._krep->momentum(nearhit->fltLen()).mag();
	beta = kres._tdef->particle().beta(mom);
	tflt = (hflt-nearhit->fltLen())/(beta*CLHEP::c_light);

					// update the time in the TrkT0 object and create a new 
					// hit object.  Assume we're at the last iteration over added error
	hitt0._t0 += tflt;
	trkhit = new TrkStrawHit(strawhit,straw,istraw,hitt0,hflt,_hiterr.back(),_maxdriftpull);
	assert(trkhit != 0);
					// allow the hit to update its own ambiguity for now: 
					// eventually we should get the resolver to do this, FIXME!!!
	//	trkhit->setAmbigUpdate(true);
					// fixed!
	trkhit->setAmbigUpdate(false);
					// must be initialy active for KalRep to process correctly
	trkhit->setActivity(true);
					// flag the added hit
	trkhit->setUsability(3);
					// add the hit to the track and the fit
	kres._krep->addHot(trkhit);
	kres._hits.push_back(trkhit);
					// create intersections for the material of this hit 
					// and add those to the track
	DetIntersection wallinter;
	if (trkhit->wallElem().reIntersect(reftraj,wallinter)) {
	  kres._krep->addInter(wallinter);
	}
	DetIntersection gasinter;
	if(trkhit->gasElem().reIntersect(reftraj,gasinter))
	  kres._krep->addInter(gasinter);

// check the raw residual: This call works because the HOT isn't yet processed as part of the fit.
//        double chi = fabs(trkhit->residual()/trkhit->hitRms());
//-----------------------------------------------------------------------------
// hit is added with _iamb=0, which means that the position uncertainty 
// includes the drift distance
//-----------------------------------------------------------------------------
	double chi, dr, sig_dr;

	dr     = trkhit->residual();
	sig_dr = sqrt(trkhit->hitRms()*trkhit->hitRms()+trkhit->driftRadius()*trkhit->driftRadius());
	chi    = dr/sig_dr;

	activity = 1;
// if it's outside limits, deactivate the HOT
	if (chi > maxchi || ! trkhit->physicalDrift(maxchi)) {
	  trkhit->setActivity(false);
	  activity = 0;
	}
	if (_debug>0){
	  printf("[KalFitHack::addHits] %5i %3i  %6i  %6i  %10.3f  %10.3f %10.3f  %10.3f \n",
		 straw.index().asInt(),
		 activity,
		 straw.id().getDevice(),
		 straw.id().getSector(),
		 trkhit->residual(),
		 trkhit->hitRms(),
		 trkhit->driftRadius(),
		 chi);
	}
// now that we've got the residual, turn off auto-ambiguity resolution
	trkhit->setAmbigUpdate(false);
      }
 // sort hits by flightlength (or in Z, which is the same)
      std::sort(kres._hits.begin(),kres._hits.end(),fltlencomp(kres._tdef->fitdir().fitDirection()));
//---------------------------------------------------------------------------
// refit the track one more time with minimal external errors
//---------------------------------------------------------------------------
      if ((_daveMode == 0) && TPeak) {
//------------------------------------------------------------------------------------------
// 2015 - 03 - 09 Gainipez added the following line for forcing the fiITeration procedure
// to use findAndUseDoublets
// 2015 - 02 - 27 Gianipez added the loop for including the external errors 
// 2015-04-03 P.Murat: not sure I understand why there are different number of iterations 
//                     in different cases - Giani?
// 2015-04-10        : perform one iteration, -1 means 'use the smallest external error defined'
//------------------------------------------------------------------------------------------
	fitIteration(kres, -1, TPeak,final);
      }
      else {
	fitIteration(kres,-1,NULL,final);
      }
      kres._krep->addHistory(kres._fit,"AddHits");
    }
  }


//-----------------------------------------------------------------------------
// loop over external hit errors, ambiguity assignment, t0 tolerance
// 10-03-2013 giani changed this loop. now it loops on all the stations
// and store the last fit that converges
// assume this is not a final call
//-----------------------------------------------------------------------------
  void KalFitHack::fitTrack(KalFitResult& KRes, CalTimePeak* TPeak) {
    int not_final(0);

    int n = _hiterr.size();
    for (int i=0; i<n; ++i) {
      fitIteration(KRes,i,TPeak,not_final);

      if (! KRes._fit.success()) break; //commented by gianipez
     }

    if(KRes._krep != 0) KRes._krep->addHistory(KRes._fit,"KalFitHack");
  }

//-----------------------------------------------------------------------------
// one step of the track fit
// update external hit errors.  This isn't strictly necessary on the 1st iteration.
//-----------------------------------------------------------------------------
  void KalFitHack::fitIteration(KalFitResult& KRes, int Iteration, CalTimePeak* TPeak, int Final) {

    double       oldt0 = KRes._krep->t0()._t0;
    unsigned     niter(0);
    bool         changed(true);
    bool         fit_success;
    //    double       extError;
    char         msg[100];

    if (_debug >0) {
      printf("------------------------------------------------------------------------------------------\n");
      printf("[KalFitHack::fitIteration] BEGIN Iteration:%i \n", Iteration);
      printf("------------------------------------------------------------------------------------------\n");
    }
//-----------------------------------------------------------------------------------
// 2015 -02 -17 G. Pezzullo: loop over the hits and assign a smaller external error 
// for the doublets
// 2015-04-03: IHErr = -1: special value, invoke last iteration
//-----------------------------------------------------------------------------------
    if (Iteration == -1) Iteration = _hiterr.size()-1;

    _annealingStep     = Iteration;
//--------------------------------------------------------------------------------
// 2015-02-19 G.Pezzu: re-search multiplets using updated fit results
// 2015-03-25 P.Murat: *TODO* I don't think this call is needed, the one in the 
//                     loop should be sufficient - check !
//-----------------------------------------------------------------------------
    KRes._nt0iter = 0;
    KRes._fit     = TrkErrCode::succeed;

    while (KRes._fit.success() && changed && niter < maxIterations()) {
      changed = false;
//-----------------------------------------------------------------------------
// set external errors and start from the standard ambiguity resolution
// if doublets are used, the doulet resolution overrides the results
// reduce external errors for hits from doublets during first iterations
// create list of doublets 
// resolve drift signs for hits in doublets. Choose the combination of drift 
// signs for which the 2-hit segment slope is the closest to that of the track 
// ** all this becomes hiden inside the ambiguity resolver** FIXME
//-----------------------------------------------------------------------------
      _ambigresolver[_annealingStep]->resolveTrk(KRes,Final);
//-----------------------------------------------------------------------------
// perform the track fit
//-----------------------------------------------------------------------------
      KRes._krep->resetFit();
      KRes.fit();

      fit_success = KRes._fit.success();

      if (_debug > 0) {
	sprintf(msg,"KalFitHack::fitIteration::002 IHErr = %2i niter = %2i success = %i",
		_annealingStep,niter,fit_success);
	printHits(KRes,msg);
      }

      if (! fit_success) break;
//-----------------------------------------------------------------------------
// if the fit succeeded, update the track T0, and recalculate the hit T0's 
//-----------------------------------------------------------------------------
      if      (_updateT0Mode == 0) {
//-----------------------------------------------------------------------------
// update T0 mode = 0: when iterating, use the cluster T0 if available
//-----------------------------------------------------------------------------
	if (TPeak != NULL)  updateCalT0(KRes,TPeak);
	else                updateT0   (KRes);
      }
      else if (_updateT0Mode == 1) {
//-----------------------------------------------------------------------------
// mode = 1: when iterating, don't look back at the cluster T0, 
//           in this mode the cluster T0 is used only to seed the process
//-----------------------------------------------------------------------------
	updateT0(KRes);
      }

      changed |= fabs(KRes._krep->t0()._t0-oldt0) > _t0tol[_annealingStep];
      oldt0    = KRes._krep->t0()._t0;
//-----------------------------------------------------------------------------
// drop outliers. weedHits() calls KalRep::fit(), so the fit success code may change
//-----------------------------------------------------------------------------
      if(_weedhits){
	KRes._nweediter = 0;
	changed        |= weedHits(KRes,Iteration,Final);
	fit_success     = KRes._fit.success();
      }
      niter++;
    }
//-----------------------------------------------------------------------------
// done iterating, define drift signs with respect to the final trajectory
// 2015-02-17 G. Pezzu: ::resolveTrk() updates drift signs of ALL hits, 
// so fix the ambiguity of the doublets after that
//-----------------------------------------------------------------------------
    _nIter += niter;

    KRes._ninter = KRes._krep->intersections();
  }


//-----------------------------------------------------------------------------
  bool KalFitHack::fitable(TrkDef const& tdef){
    return tdef.strawHitIndices().size() >= _minnstraws;
  }
  
//-----------------------------------------------------------------------------
  void KalFitHack::makeHits(KalFitResult& KRes, TrkT0 const& t0) {
    const Tracker& tracker = getTrackerOrThrow();
    TrkDef const& tdef = *KRes._tdef;
// compute the propagaion velocity
    double flt0 = tdef.helix().zFlight(0.0);
    double mom = TrkMomCalculator::vecMom(tdef.helix(),bField(),flt0).mag();
    double vflt = tdef.particle().beta(mom)*CLHEP::c_light;
    unsigned nind = tdef.strawHitIndices().size();
    for(unsigned iind=0;iind<nind;iind++){
      size_t istraw = tdef.strawHitIndices()[iind]._index;
      const StrawHit& strawhit(tdef.strawHitCollection()->at(istraw));
      const Straw& straw = tracker.getStraw(strawhit.strawIndex());
      double fltlen = tdef.helix().zFlight(straw.getMidPoint().z());
    // estimate arrival time at the wire
      TrkT0 hitt0(t0);
      hitt0._t0 += (fltlen-flt0)/vflt;
    // create the hit object.  Start with the 1st additional error for anealing
      TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,hitt0,fltlen,_hiterr.front(),_maxdriftpull);
      assert(trkhit != 0);
    // set the initial ambiguity based on the input
      trkhit->setAmbig(tdef.strawHitIndices()[iind]._ambig);
    // refine the flightlength, as otherwise hits in the same plane are at exactly the same flt, which can cause problems
      TrkErrCode pstat = trkhit->updatePoca(&tdef.helix());
      if(pstat.failure()){
        trkhit->setActivity(false);
      }
      KRes._hits.push_back(trkhit);
    }
 // sort the hits by flightlength
    std::sort(KRes._hits.begin(),KRes._hits.end(),fltlencomp(tdef.fitdir().fitDirection()));
  }

//-----------------------------------------------------------------------------
// '*' in front of the hit drift radius: the hit drift sign is not defined 
//     and the drift ambiguity has been set to 0
// '?': the drift sign determined by the resolver is different from the MC truth
//-----------------------------------------------------------------------------
  void KalFitHack::printHits(KalFitResult& KRes, const char* Caller) {
    const KalRep* Trk  = KRes._krep;

    if (Trk != 0) {

      printf("[KalFitHack::printHits] BEGIN called from %s iherr:%i \n",Caller,_annealingStep);
      printf("---------------------------------------------------------------------------------");
      printf("-----------------------------------------------------\n");
      //      printf("%s",Prefix);
      printf("  TrkID       Address    N  NA      P       pT     costh    T0      T0Err   Omega");
      printf("      D0       Z0      Phi0   TanDip    Chi2    FCons\n");
      printf("---------------------------------------------------------------------------------");
      printf("-----------------------------------------------------\n");

      Hep3Vector trk_mom;

      double h1_fltlen = Trk->firstHit()->kalHit()->hitOnTrack()->fltLen() - 10;
      trk_mom          = Trk->momentum(h1_fltlen);
      double mom       = trk_mom.mag();
      double pt        = trk_mom.perp();
      double costh     = trk_mom.cosTheta();
      double chi2      = Trk->chisq();

      int    nhits(0);

      const TrkHotList* hots = Trk->hotList();
      for (TrkHotList::hot_iterator ihot=hots->begin(); ihot != hots->end(); ++ihot) {
	nhits++;
      }

      int    nact      = Trk->nActive();
      double t0        = Trk->t0().t0();
      double t0err     = Trk->t0().t0Err();
//-----------------------------------------------------------------------------
// in all cases define momentum at lowest Z - ideally, at the tracker entrance
//-----------------------------------------------------------------------------
      double s1     = Trk->firstHit()->kalHit()->hitOnTrack()->fltLen();
      double s2     = Trk->lastHit ()->kalHit()->hitOnTrack()->fltLen();
      double s      = std::min(s1,s2);

      double d0     = Trk->helix(s).d0();
      double z0     = Trk->helix(s).z0();
      double phi0   = Trk->helix(s).phi0();
      double omega  = Trk->helix(s).omega();
      double tandip = Trk->helix(s).tanDip();

      double fit_consistency = Trk->chisqConsistency().consistency();
      int q         = Trk->charge();
      
      printf("%5i %16p %3i %3i %8.3f %7.3f %8.4f %7.3f %7.4f",
	     -1,
	     Trk,
	     nhits,
	     nact,
	     q*mom,pt,costh,t0,t0err
	     );

      printf(" %8.5f %8.3f %8.3f %8.4f %7.4f",
	     omega,d0,z0,phi0,tandip
	     );
      printf(" %8.3f %10.3e\n",
	     chi2,
	     fit_consistency);
    }
    
    //-----------------------------------------------------------------------------
    // print detailed information about the track hits
    //-----------------------------------------------------------------------------
    //    const TrkHotList* hot_list = Trk->hotList();
    int nhits = KRes._hits.size();
    printf("--------------------------------------------------------------------");
    printf("----------------------------------------------------------------");
    printf("--------------------------------------------\n");
    printf(" ih  SInd U A     len         x        y        z      HitT    HitDt");
    printf(" Ch Pl  L  W     T0       Xs      Ys        Zs      resid  sigres");
    printf("   Rdrift   mcdoca  totErr hitErr  t0Err penErr extErr\n");
    printf("--------------------------------------------------------------------");
    printf("----------------------------------------------------------------");
    printf("--------------------------------------------\n");

    mu2e::TrkStrawHit     *hit;
    Hep3Vector            pos;
    const mu2e::StrawHit  *sh;
    const mu2e::Straw     *straw;
    int                   ihit, volume_id, nstraws;
    double                len;
    HepPoint              plen;

    ihit = 0;
    for (int it=0; it<nhits; ++it) {
      hit   =  KRes._hits.at(it);
      sh    = &hit->strawHit();
      straw = &hit->straw();

      hit->hitPosition(pos);

      len   = hit->fltLen();
      plen  = Trk->position(len);
//-----------------------------------------------------------------------------
// find MC truth DOCA in a given straw
// start from finding the right vector of StepPointMC's
//-----------------------------------------------------------------------------
      nstraws = fListOfMCStrawHits->size();

      const mu2e::StepPointMC* step(0);

      for (int i=0; i<nstraws; i++) {
	const mu2e::PtrStepPointMCVector&  mcptr(fListOfMCStrawHits->at(i));
	step = &(*mcptr.at(0));
	volume_id = step->volumeId();
 	if (volume_id == straw->index().asInt()) {
//-----------------------------------------------------------------------------
// step found - use the first one in the straw
//-----------------------------------------------------------------------------
 	  break;
 	}
      }

      double mcdoca = -99.0;

      if (step) {
	const Hep3Vector* v1 = &straw->getMidPoint();
	HepPoint p1(v1->x(),v1->y(),v1->z());

	const Hep3Vector* v2 = &step->position();
	HepPoint    p2(v2->x(),v2->y(),v2->z());

	TrkLineTraj trstraw(p1,straw->getDirection()  ,0.,0.);
	TrkLineTraj trstep (p2,step->momentum().unit(),0.,0.);

	TrkPoca poca(trstep, 0., trstraw, 0.);
    
	mcdoca = poca.doca();
      }

      ihit += 1;
      printf("%3i %5i %1i %1i %9.3f %8.3f %8.3f %9.3f %8.3f %7.3f",
	     ihit,
	     straw->index().asInt(), 
	     hit->isUsable(),
	     hit->isActive(),
	     len,
	     //	     hit->hitRms(),
	     plen.x(),plen.y(),plen.z(),
	     sh->time(), sh->dt()
	     );

      printf(" %2i %2i %2i %2i",
	     straw->id().getDevice(),
	     straw->id().getSector(),
	     straw->id().getLayer(),
	     straw->id().getStraw()
	     );

      printf(" %8.3f",hit->hitT0().t0());

      double res, sigres;
      hit->resid(res, sigres, true);

      printf(" %8.3f %8.3f %9.3f %7.3f %7.3f",
	     pos.x(),
	     pos.y(),
	     pos.z(),
	     res,
	     sigres
	     );

      if      (hit->ambig()       == 0) printf(" * %6.3f",hit->driftRadius());
      else if (hit->ambig()*mcdoca > 0) printf("   %6.3f",hit->driftRadius()*hit->ambig());
      else                              printf(" ? %6.3f",hit->driftRadius()*hit->ambig());

      printf("  %7.3f  %6.3f %6.3f %6.3f %6.3f %6.3f\n",		 
	     mcdoca, 
	     hit->totalErr(),
	     hit->hitErr(),
	     hit->t0Err(),
	     hit->penaltyErr(),
	     hit->extErr()
	     );
    }
  }
 

//-----------------------------------------------------------------------------
  void KalFitHack::makeMaterials(KalFitResult& KRes) {
    TrkDef const& tdef = *KRes._tdef;
    for(std::vector<TrkStrawHit*>::iterator ihit=KRes._hits.begin();ihit!=KRes._hits.end();ihit++){
      TrkStrawHit* trkhit = *ihit;
      // create wall and gas intersection objects from each straw hit (active or not)
      DetIntersection wallinter;
      wallinter.delem = 0;
      wallinter.pathlen = trkhit->fltLen();
      DetIntersection gasinter;
      gasinter.delem = 0;
      gasinter.pathlen = trkhit->fltLen();
      if(trkhit->wallElem().reIntersect(&tdef.helix(),wallinter))
	KRes._detinter.push_back(wallinter);
      if(trkhit->gasElem().reIntersect(&tdef.helix(),gasinter))
	KRes._detinter.push_back(gasinter);
    }
  }

//-----------------------------------------------------------------------------
  bool KalFitHack::weedHits(KalFitResult& KRes, int Iteration, int Final) {
    // Loop over HoTs and find HoT with largest contribution to chi2.  If this value
    // is greater than some cut value, deactivate that HoT and reFit
    bool retval(false);
    double worst = -1.;
    TrkStrawHit* worstHot = 0;
    for (std::vector<TrkStrawHit*>::iterator iter = KRes._hits.begin(); iter != KRes._hits.end(); ++iter){
      TrkStrawHit* iHot = *iter;
      if (iHot->isActive()) {
        double resid, residErr;
        if(iHot->resid(resid, residErr, true)){
          double value = fabs(resid/residErr);
          if (value > _maxhitchi && value > worst) {
            worst    = value;
            worstHot = iHot;
          }
        }
      }
    }
    if(0 != worstHot){
      retval = true;
      worstHot->setActivity(false);
      worstHot->setUsability(5); // positive usability allows hot to be re-enabled later

      _ambigresolver[Iteration]->resolveTrk(KRes,Final);

      KRes.fit();
      KRes._krep->addHistory(KRes._fit, "HitWeed");
      // Recursively iterate
      KRes._nweediter++;
      if (KRes._fit.success() && KRes._nweediter < _maxweed ) {
        retval |= weedHits(KRes,Iteration,Final);
      }
    }
    return retval;
  }
  

//-----------------------------------------------------------------------------
// assume that this is a call during the last iteration
//-----------------------------------------------------------------------------
  bool KalFitHack::unweedHits(KalFitResult& KRes, double maxchi) {
    // Loop over inactive HoTs and find the one with the smallest contribution to chi2.  If this value
    // is less than some cut value, reactivate that HoT and reFit
    bool       retval(false);
    int  const final (1);
    double     best = 1.e12;
    int        last_iteration = _hiterr.size()-1;

    TrkStrawHit* bestHot = 0;
    for (std::vector<TrkStrawHit*>::iterator iter = KRes._hits.begin(); iter != KRes._hits.end(); ++iter){
      TrkStrawHit* iHot = *iter;
      if (!iHot->isActive()) {
        double resid, residErr;
        if(iHot->resid(resid, residErr, true)){
          double chival = fabs(resid/residErr);
  // test both for a good chisquared and for the drift radius to be physical
          if (chival < maxchi && iHot->physicalDrift(maxchi) && chival < best) {
            best    = chival;
            bestHot = iHot;
          }
        }
      }
    }
    if (0 != bestHot) {
      retval = true;
      bestHot->setActivity(true);
      bestHot->setUsability(4);
//-----------------------------------------------------------------------------
// update drift signs before the fit again - more doublets could've been recovered
// one more place for the ambiguity resolver to be called from
// also - from addHits
// here we call fitter not invoking the fitIteration
//-----------------------------------------------------------------------------
//      KRes._decisionMode = _decisionMode;
      _ambigresolver[last_iteration]->resolveTrk(KRes,final);

      KRes.fit();
      KRes._krep->addHistory(KRes._fit, "HitUnWeed");
      // Recursively iterate
      KRes._nunweediter++;
      if (KRes._fit.success() && KRes._nunweediter < _maxweed  ) {
        retval |= unweedHits(KRes,maxchi);
      }
    }
    return retval;
  }
  
  BField const&  KalFitHack::bField() const {
    if(_bfield == 0){
      GeomHandle<BFieldConfig> bfconf;
      if(_fieldcorr){
// create a wrapper around the mu2e field 
	_bfield = new BaBarMu2eField();	
      } else {
// create a fixed field using the nominal value
	GeomHandle<BFieldConfig> bfconf;
	_bfield=new BFieldFixed(bfconf->getDSUniformValue());
	assert(_bfield != 0);
      }
    }
    return *_bfield;
  }
   
  const TrkVolume* KalFitHack::trkVolume(trkDirection trkdir) const {
    //FIXME!!!!
    return 0;
  }

//-----------------------------------------------------------------------------
// time initialization
//-----------------------------------------------------------------------------
  void KalFitHack::initCaloT0(CalTimePeak* TPeak, TrkDef const& tdef, TrkT0& t0) {
//    2014-11-24 gianipez and Pasha removed time offset between caloriemter and tracker

    // get flight distance of z=0
    double t0flt = tdef.helix().zFlight(0.0);
    // estimate the momentum at that point using the helix parameters.  This is
    // assumed constant for this crude estimate
    double mom = TrkMomCalculator::vecMom(tdef.helix(),bField(),t0flt).mag();
    // compute the particle velocity
    double vflt = tdef.particle().beta(mom)*CLHEP::c_light;
//-----------------------------------------------------------------------------
// Calculate the path length of the particle from the middle of the Tracker to the 
// calorimeter, TPeak->Z() is calculated wrt the tracker center 
//-----------------------------------------------------------------------------
    double path = TPeak->ClusterZ()/tdef.helix().sinDip();

    t0._t0 = TPeak->ClusterT0() + _dtoffset - path/vflt;
    
    //Set dummy error value
    t0._t0err = 1.;
  }


//-----------------------------------------------------------------------------
  void KalFitHack::initT0(TrkDef const& tdef, TrkT0& t0) {
    using namespace boost::accumulators;
// make an array of all the hit times, correcting for propagation delay
    const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    unsigned nind = tdef.strawHitIndices().size();
    std::vector<double> times;
    times.reserve(nind);
    // get flight distance of z=0
    double t0flt = tdef.helix().zFlight(0.0);
    // estimate the momentum at that point using the helix parameters.  This is
    // assumed constant for this crude estimate
    double mom = TrkMomCalculator::vecMom(tdef.helix(),bField(),t0flt).mag();
    // compute the particle velocity
    double vflt = tdef.particle().beta(mom)*CLHEP::c_light;
    // for crude estimates, we only need 1 d2t function
    D2T d2t;
    static CLHEP::Hep3Vector zdir(0.0,0.0,1.0);
    // loop over strawhits
    for(unsigned iind=0;iind<nind;iind++){
      size_t istraw = tdef.strawHitIndices()[iind]._index;
      const StrawHit& strawhit(tdef.strawHitCollection()->at(istraw));
      const Straw& straw = tracker.getStraw(strawhit.strawIndex());
      // compute the flightlength to this hit from z=0 (can be negative)
      double hflt = tdef.helix().zFlight(straw.getMidPoint().z()) - t0flt;
      // Use this to estimate the time for the track to reaches this hit from z=0
      double tprop = hflt/vflt;
      // estimate signal propagation time on the wire assuming the middle (average)
      double vwire = tcal->SignalVelocity(straw.index());
      double teprop = straw.getHalfLength()/vwire;
      // correct the measured time for these effects: this gives the aveage time the particle passed this straw, WRT
      // when the track crossed Z=0
    // assume the average drift time is half the maximum drift distance.  This is a poor approximation, but good enough for now
      if(iind==0)tcal->DistanceToTime(straw.index(),0.5*straw.getRadius(),zdir,d2t);
      double htime = strawhit.time() - tprop - teprop - d2t._tdrift;
      times.push_back(htime);
    }
    // find the median time
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > med;
    med = std::for_each( times.begin(), times.end(), med );
    t0._t0 = extract_result<tag::median>(med);
    accumulator_set<double, stats<tag::min> >  min;
    accumulator_set<double, stats<tag::max> > max;
    min = std::for_each( times.begin(), times.end(), min );
    max = std::for_each( times.begin(), times.end(), max );
    double tmin = extract_result<tag::min>(min);
    double tmax = extract_result<tag::max>(max);
    // estimate the error using the range
    t0._t0err = (tmax-tmin)/sqrt(12*nind);
  }

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  void KalFitHack::updateCalT0(KalFitResult& KRes, CalTimePeak* TPeak) {
//    2014-11-24 gianipez and Pasha removed time offset between caloriemter and tracker

    TrkT0 t0;
    double mom, vflt, path, t0flt, flt0(0.0);
    bool converged = TrkHelixUtils::findZFltlen(KRes._krep->traj(),0.0,flt0);
    
    //get helix from kalrep
    HelixTraj trkHel(KRes._krep->helix(flt0).params(),KRes._krep->helix(flt0).covariance());
    
					// get flight distance of z=0
    t0flt = trkHel.zFlight(0.0);
    
    if (converged) {
//-----------------------------------------------------------------------------
// estimate the momentum at that point using the helix parameters.  
// This is assumed constant for this crude estimate
// compute the particle velocity
//-----------------------------------------------------------------------------
      mom  = TrkMomCalculator::vecMom(trkHel,bField(),t0flt).mag();
      vflt = KRes._tdef->particle().beta(mom)*CLHEP::c_light;
//-----------------------------------------------------------------------------
// path length of the particle from the middle of the Tracker to the  calorimeter
// set dummy error value
//-----------------------------------------------------------------------------
      path      = TPeak->ClusterZ()/trkHel.sinDip();
      t0._t0    = TPeak->ClusterT0() + _dtoffset - path/vflt;
      t0._t0err = 1.;
      
      KRes._krep->setT0(t0,flt0);
      updateHitTimes(KRes);
    }
  }
  

//-----------------------------------------------------------------------------
  bool KalFitHack::updateT0(KalFitResult& KRes) {
    using namespace boost::accumulators;
    bool retval(false);
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    KalRep* krep = KRes._krep;
// need to have a valid fit
    if(krep->fitValid()) {
// find the global fltlen associated with z=0. 
      double flt0(0.0);
      bool converged = TrkHelixUtils::findZFltlen(krep->traj(),0.0,flt0);
      if(converged){
	std::vector<double> hitt0; // store t0, to allow outlyer removal
	std::vector<double> hitt0err;
	hitt0.reserve(KRes._hits.size());
	hitt0err.reserve(KRes._hits.size());
	// loop over the hits
	for(std::vector<TrkStrawHit*>::iterator ihit= KRes._hits.begin();ihit != KRes._hits.end(); ihit++){
	  TrkStrawHit* hit = *ihit;
	  if(hit->isActive() && hit->poca()!= 0 && hit->poca()->status().success()){
	    // find the residual, exluding this hits measurement
	    double resid,residerr;
	    if(krep->resid(hit,resid,residerr,true)){
	      // convert this to a distance to the wire
	      double doca = (resid + hit->driftRadius()*hit->ambig());
	      if(hit->ambig() == 0)
		doca = fabs(doca);
	      else
		doca *= hit->ambig();
	      // restrict the range, symmetrically to avoid bias
	      double rad = hit->straw().getRadius();
	      if(doca > _mint0doca && doca < rad-_mint0doca){
		// translate the DOCA into a time
		D2T d2t;
		tcal->DistanceToTime(hit->straw().index(),doca,krep->traj().direction(hit->fltLen()),d2t);
		// subtracting hitT0 makes this WRT the previous track t0
		hitt0.push_back(hit->time() - d2t._tdrift - hit->signalTime() - hit->hitT0()._t0);
		// assume residual error dominates
		hitt0err.push_back(residerr/d2t._vdrift);
	      }
	    }
	  }
	}
	if(hitt0.size() >1){
	  TrkT0 t0;
	  // find the median
	  accumulator_set<double, stats<tag::median(with_p_square_quantile) > > med;
	  med = std::for_each( hitt0.begin(), hitt0.end(), med );
	  t0._t0 = extract_result<tag::median>(med);
	  // iterate an outlier search and linear fit until the set of used hits doesn't change
	  bool changed(true);
	  std::vector<bool> used(hitt0.size(),true);
	  unsigned niter(0);
	  while(changed && niter < 10){
	    niter++;
	    changed = false;
	    accumulator_set<double,stats<tag::weighted_variance>,double > wmean;
	    for(unsigned ihit=0;ihit<hitt0.size();ihit++){
	      bool useit = fabs(hitt0[ihit]-t0._t0) < _t0nsig*hitt0err[ihit];
	      changed |= useit != used[ihit];
	      used[ihit] = useit;
	      if(useit){
		wmean(hitt0[ihit], weight=1.0/(hitt0err[ihit]*hitt0err[ihit]));
	      }
	    }
	    unsigned nused = extract_result<tag::count>(wmean);
	    if(nused > 1){
	      t0._t0 = extract_result<tag::weighted_mean>(wmean);
	      t0._t0err = sqrt(extract_result<tag::weighted_variance>(wmean)/nused);
	    } else {
	      break;
	    }
	  }
	  // reset t0
	  if(!changed){
	    // put in t0 from the track.
	    t0._t0 += krep->t0()._t0;
	    krep->setT0(t0,flt0);
	    updateHitTimes(KRes);
	    retval = true;
	  }
	}
      }
    }
    return retval;
  }

//-----------------------------------------------------------------------------
  void KalFitHack::updateHitTimes(KalFitResult& KRes) {
  // compute the time the track came closest to the wire for each hit, starting from t0 and working out.
  // this function allows for momentum change along the track.
  // find the bounding hits on either side of this
    std::sort(KRes._hits.begin(),KRes._hits.end(),fltlencomp(KRes._tdef->fitdir().fitDirection()));
    std::vector<TrkStrawHit*>::iterator ihigh;
    std::vector<TrkStrawHit*>::reverse_iterator ilow;
    findBoundingHits(KRes._hits,KRes._krep->flt0(),ilow,ihigh);
    // reset all the hit times
    double hflt = KRes._krep->flt0();
    TrkT0 hitt0 = KRes._krep->t0();
    for(std::vector<TrkStrawHit*>::iterator ihit= ihigh;ihit != KRes._hits.end(); ++ihit){
      TrkStrawHit* hit = *ihit;
// particle momentum at this point, using the full fit
      double mom = KRes._krep->momentum(hit->fltLen()).mag();
// relativistic velocity from that
      double beta = KRes._tdef->particle().beta(mom);
// particle transit time to this hit from the reference
      double tflt = (hit->fltLen()-hflt)/(beta*CLHEP::c_light);
// update the time in the TrkT0 object
      hitt0._t0 += tflt;
      (*ihit)->updateHitT0(hitt0);
// update the reference flightlength
      hflt = hit->fltLen();
    }
//     if (_debug > 1) {
//       printf("[KalFitHack::updateHitTimes] moving forward\n");
//       printHits(KRes,"updateTimes_001");
//     }

// now the same, moving backwards
    hflt = KRes._krep->flt0();
    hitt0 = KRes._krep->t0();
    for(std::vector<TrkStrawHit*>::reverse_iterator ihit= ilow;ihit != KRes._hits.rend(); ++ihit){
      TrkStrawHit* hit = *ihit;
      double mom = KRes._krep->momentum(hit->fltLen()).mag();
      double beta = KRes._tdef->particle().beta(mom);
      double tflt = (hit->fltLen()-hflt)/(beta*CLHEP::c_light);
      hitt0._t0 += tflt;
      (*ihit)->updateHitT0(hitt0);
      hflt = hit->fltLen();
    }

//     if (_debug > 1) {
//       printf("[KalFitHack::updateHitTimes] moving backwards\n");
//       printHits(KRes,"updateTimes_002");
//     }
  }

//-----------------------------------------------------------------------------
  void KalFitHack::findBoundingHits(std::vector<TrkStrawHit*>& hits,double flt0,
    std::vector<TrkStrawHit*>::reverse_iterator& ilow,
    std::vector<TrkStrawHit*>::iterator& ihigh) {
    ilow = hits.rbegin();
    ihigh = hits.begin();
    while(ilow != hits.rend() && (*ilow)->fltLen() > flt0 )++ilow;
    while(ihigh != hits.end() && (*ihigh)->fltLen() < flt0 )++ihigh;
  }

//-----------------------------------------------------------------------------
// stolen from AmbigResolver class - 
//-----------------------------------------------------------------------------
//   const TrkSimpTraj* KalFitHack::findTraj(std::vector<TrkStrawHit*> const& Hits, 
// 					  const KalRep*                    Krep) const {
    
//     typedef std::vector<KalSite*>::const_iterator KSI;

//     const TrkSimpTraj* retval(0);
// // if the fit is valid, use the full fit result.  Otherwise, use the reference traj
//     if(!Krep->fitValid()){
//       double locdist;
//       retval = Krep->referenceTraj()->localTrajectory(Hits[0]->fltLen(),locdist);
//     } else {
// // find the range of these hits in the KalRep site vector
//       std::vector<KalSite*> const& sites = Krep->siteList();
//       KSI first = sites.end();
//       KSI last  = sites.begin();
//       for(auto ihit = Hits.begin();ihit != Hits.end();++ihit){
// 	if((*ihit)->isActive()){
// 	  const KalHit* hitsite = Krep->findHotSite(*ihit);
// 	  // find the index to this site
// 	  KSI ifnd = std::find(sites.begin(),sites.end(),hitsite);
// 	  if(ifnd != sites.end()) {
// 	    if(ifnd < first){
// 	      first = ifnd;
// 	    }
// 	    if(ifnd > last){
// 	      last = ifnd;
// 	    }
// 	  }
// 	}
//       }
// // create a trajectory from the fit which excludes this set of hits
//       static TrkSimpTraj* straj = Krep->seed()->clone();
//       if(Krep->smoothedTraj(first,last,straj)){
// 	retval = straj;
//       } else {
// 	double locdist;
// 	retval = Krep->referenceTraj()->localTrajectory(Hits[0]->fltLen(),locdist);
//       }
//     }
//     return retval;
//   }

} // namespace



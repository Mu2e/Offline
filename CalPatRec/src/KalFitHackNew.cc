///////////////////////////////////////////////////////////////////////////////
// $Id: KalFitHackNew.cc,v 1.4 2014/04/08 04:25:46 murat Exp $
// $Author: murat $
// $Date: 2014/04/08 04:25:46 $
//
// the following data members have to be initialized at construction time by the
// code creating the KalPatRec object:
//
// KalFitHackNew::_tracker
// KalFitHackNew::_tcal
///////////////////////////////////////////////////////////////////////////////
// framework
#include "fhiclcpp/ParameterSet.h"
// the following has to come before other BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "CalPatRec/inc/KalFitHackNew.hh"

#include "BTrkData/inc/TrkStrawHit.hh"
#include "TrkReco/inc/PanelAmbigResolver.hh"
#include "TrkReco/inc/PocaAmbigResolver.hh"
#include "TrkReco/inc/HitAmbigResolver.hh"
#include "TrkReco/inc/FixedAmbigResolver.hh"
#include "TrkReco/inc/DoubletAmbigResolver.hh"
#include "Mu2eBTrk/inc/BaBarMu2eField.hh"
#include "Mu2eBTrk/inc/Mu2eDetectorModel.hh"
//geometry
#include "BTrkHelper/inc/BTrkHelper.hh"
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
// #include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
// BaBar
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
//09 - 26 - 2013 gianipez added the following include file
#include "BTrk/TrkBase/HelixParams.hh"
//----------------------------------------
#include "BTrk/KalmanTrack/KalMaterial.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/BaBar/ErrLog.hh"
#include "BTrk/BField/BFieldFixed.hh"
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "BTrk/DetectorModel/DetMaterial.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
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

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {
// comparison functor for ordering hits
  struct fltlencomp : public binary_function<TrkHit*, TrkHit*, bool> {
    fltlencomp(TrkFitDirection::FitDirection fdir=TrkFitDirection::downstream) : _fdir(fdir) {}
    bool operator()(TrkHit* x, TrkHit* y) {
      return _fdir == TrkFitDirection::downstream ? x->fltLen() < y->fltLen() : y->fltLen() < x->fltLen() ;
    }
    TrkFitDirection::FitDirection _fdir;
  };

// struct for finding materials
  struct StrawFlight {
    StrawIndex _index;  // straw being tested
    double _flt; // flight where trajectory comes near this straw
// construct from pair
    StrawFlight(StrawIndex strawind, double flt) : _index(strawind), _flt(flt) {}
  };

// comparison operators understand that the same straw could be hit twice, so the flight lengths need
// to be similar befoew we consider these 'the same'
  struct StrawFlightComp : public binary_function<StrawFlight, StrawFlight, bool> {
    double _maxdiff; // maximum flight difference; below this, consider 2 intersections 'the same'
    StrawFlightComp(double maxdiff) : _maxdiff(maxdiff) {}
    bool operator () (StrawFlight const& a, StrawFlight const& b) { return a._index < b._index ||
    ( a._index == b._index && a._flt < b._flt && fabs(a._flt-b._flt)>=_maxdiff);}
  };

// construct from a parameter set
  KalFitHackNew::KalFitHackNew(fhicl::ParameterSet const& pset) :
// KalFitHackNew parameters
    _debugLevel   (pset.get<int>     ("debugLevel")),
    _minnstraws   (pset.get<unsigned>("minnstraws")),
    _maxmatfltdiff(pset.get<double>  ("MaximumMaterialFlightDifference")), // mm separation in flightlength
    _weedhits(pset.get<vector<bool>>("weedhits")),
    _maxhitchi(pset.get<double>     ("maxhitchi")),
    _maxweed(pset.get<unsigned>       ("maxweed")),
    _hiterr(pset.get<vector<double>>  ("hiterr")),
    _maxdriftpull(pset.get<double>    ("maxDriftPull")),
    // t0 parameters
	    //    _initt0(pset.get<bool>            ("initT0",true)),
    _updateT0(pset.get<bool>          ("updateT0",true)),
    _updateT0Mode(pset.get<int>       ("updateT0Mode")),
    _minHitDrift(pset.get<double>     ("HitMinDrift")),
    fRdriftMinusDocaTol(pset.get<double>("RdriftMinusDocaTol")),
    _t0tol(pset.get< vector<double> >   ("t0Tolerance")),
    //    _t0errfac(pset.get<double>          ("t0ErrorFactor")),
    _mint0doca(pset.get<double>         ("minT0DOCA")),
    _t0nsig(pset.get<double>            ("t0window")),
    _dtoffset(pset.get<double>          ("dtOffset")),
    fScaleErrDoublet(pset.get<double>   ("scaleErrDoublet")),
    //    fUseDoublets(0),
    fMinDriftDoublet  (pset.get<double>      ("minDriftDoublet")),
    fDeltaDriftDoublet(pset.get<double>      ("deltaDriftDoublet")),
    _maxDoubletChi2   (pset.get<double>      ("maxDoubletChi2")),
    fSigmaSlope       (pset.get<double>      ("sigmaSlope")),
    fMakeStrawHitModuleLabel(pset.get<string>("makeStrawHitModuleLabel")),
    //
    _removefailed(pset.get<bool>             ("removeFailedFits")),
    _ambigstrategy(pset.get<vector<int>>     ("ambiguityStrategy")),
    _addmaterial(pset.get<vector<bool>>      ("AddMaterial")),
    _bfield(0),
    _nIter(0)
  {
    _darPset = new fhicl::ParameterSet(pset.get<fhicl::ParameterSet>("DoubletAmbigResolver",fhicl::ParameterSet()));

// set KalContext parameters
    _disttol      = pset.get<double>("IterationTolerance",0.1);
    _intertol     = pset.get<double>("IntersectionTolerance",100.0);
    _maxiter      = pset.get<long>("MaxIterations",10);
    _maxinter     = pset.get<long>("MaxIntersections",0);
    _matcorr      = pset.get<bool>("materialCorrection",true);
    _fieldcorr    = pset.get<bool>("fieldCorrection",false);
    _smearfactor  = pset.get<double>("SeedSmear",1.0e6);
    _sitethresh   = pset.get<double>("SiteMomThreshold",0.2);
    _momthresh    = pset.get<double>("MomThreshold",10.0);
    _mingap       = pset.get<double>("mingap",0.1);
    _minfltlen    = pset.get<double>("MinFltLen",0.1);
    _minmom       = pset.get<double>("MinMom",10.0);
    _fltepsilon   = pset.get<double>("FltEpsilon",0.001);
    _divergeflt   = pset.get<double>("DivergeFlt");
    _mindot       = pset.get<double>("MinDot",0.0);
    _maxmomdiff   = pset.get<double>("MaxMomDiff",0.5);
    _momfac       = pset.get<double>("MomFactor",0.0);
    _maxpardif[0] = _maxpardif[1] = pset.get<double>("MaxParameterDifference",1.0);
    // DOF counting subdivision is illogical, FIXME!!!!
    _mindof       = pset.get<double>("MinNDOF");
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
        << "mu2e::KalFitHackNew: no hit errors specified" << endl;
    }

    if (_hiterr.size() != _ambigstrategy.size()) {
      throw cet::exception("RECO")
        << "mu2e::KalFitHackNew: inconsistent ambiguity resolution" << endl;
    }

    if (_hiterr.size() != _t0tol.size()) {
      throw cet::exception("RECO")
        <<"mu2e::KalFitHackNew: inconsistent ambiguity resolution" << endl;
    }

    if (_hiterr.size() != _weedhits.size()) {
      throw cet::exception("RECO")
        <<"mu2e::KalFitHackNew: inconsistent _weedhits size" << endl;
    }
    if (_hiterr.size() != _addmaterial.size()) {
      throw cet::exception("RECO") 
	<< "mu2e::KalFit: inconsistent ambiguity resolution AddMaterial" << endl;
    } 
    AmbigResolver* ar;
    double         err;
    int            n, final(0);

    n = _ambigstrategy.size();
    for(int i=0; i<n; ++i) {
      err = _hiterr[i];
      switch (_ambigstrategy[i]) {
      case kFixedAmbig: default:
        ar = new FixedAmbigResolver(pset,err);
        break;
      case kHitAmbig:
        ar = new HitAmbigResolver(pset,err);
        break;
      case kPanelAmbig:
        ar = new PanelAmbig::PanelAmbigResolver(pset,err,i);
        break;
      case kPocaAmbig:
        ar = new PocaAmbigResolver(pset,err);
        break;
      case kDoubletAmbig: // 4
        ar = new DoubletAmbigResolver(*_darPset,err,i,final);
        break;
      }
      _ambigresolver.push_back(ar);
    }
//-----------------------------------------------------------------------------
// 2016-01-27 P.Murat: for DoubletAmbigResolver need one more instance with 
// 'final' set to 1 - now need to figure how to use it
// assume that KalFitHackNew is using the DoubletAmbigResolver only
//-----------------------------------------------------------------------------
    err = _hiterr[n-1];
    ar    = new DoubletAmbigResolver(*_darPset,err,n,final);
    _ambigresolver.push_back(ar);
//-----------------------------------------------------------------------------
// print routine
//-----------------------------------------------------------------------------
    _mcTruth = pset.get <int >("mcTruth"); 

    if (_mcTruth != 0) {
      fhicl::ParameterSet ps = pset.get<fhicl::ParameterSet>("mcUtils");
      _mcUtils = art::make_tool  <McUtilsToolBase>(ps);
    }
    else               _mcUtils = std::make_unique<McUtilsToolBase>();
  }

//-----------------------------------------------------------------------------
  KalFitHackNew::~KalFitHackNew(){
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
  void KalFitHackNew::makeTrack(KalFitResultNew& Result) {

    Result.fit = TrkErrCode(TrkErrCode::fail);
//-----------------------------------------------------------------------------
// first, find t0
//-----------------------------------------------------------------------------
    const CaloCluster* cluster = Result.caloCluster;

    initT0(Result);
//-----------------------------------------------------------------------------
// knowing t0, create the BaBar hit list, and fill. 
// The BaBar list takes ownership
// This will go away when we cleanup the BaBar hit storage, FIXME!!!
// Find the wall and gas material description objects for these hits
// TrkStrawHit inherits from TrkHit
//-----------------------------------------------------------------------------
    vector<DetIntersection> detinter;       // material intersections, used by the KalRep
    vector<TrkHit*>    hotlist;

    makeHits(Result,hotlist);
    if (_matcorr) makeMaterials(hotlist,*Result.helixTraj,detinter);

// create Kalman representation

    double flt0 = Result.helixTraj->zFlight(0.0);
    Result.krep = new KalRep(*Result.helixTraj, hotlist, detinter, *this, Result.tpart,Result.t0,flt0);

    hotlist.clear();         // clear hotlist, as from now on the hits are owned by Result._krep

    assert(Result.krep != 0);

// initialize krep t0; eventually, this should be in the constructor, FIXME!!!

    // double flt0 = Result.helixTraj()->zFlight(0.0);
    // Result._krep->setT0(Result._t0,flt0);

    if (_debugLevel > 0) printHits(Result,"makeTrack_001");
//-----------------------------------------------------------------------------
// FIT. 10-07-2013 giani added the following line. It updates the hit times
//                 following changes in the t0 value
//-----------------------------------------------------------------------------
    if (cluster) updateHitTimes(Result.krep);
//-----------------------------------------------------------------------------
// 09 - 26 - 2013 giani: include the calorimeter information when it is avaiable
//-----------------------------------------------------------------------------
      fitTrack(Result);

    if (_removefailed) Result.removeFailed();
  }

//-----------------------------------------------------------------------------
// addHits is called _only_once_ from CalPatRec_module::produce after 
// the last iteration
//-----------------------------------------------------------------------------
  void KalFitHackNew::addHits(KalFitResultNew& KFRes, double MaxChi) {

    const char* oname = "[KalFitHackNew::addHits]";

                                        // there must be a valid Kalman fit to add hits to
    int  activity(0);
					// last iteration
    int    iteration = _hiterr.size();
    double hit_error = _hiterr[iteration-1];

    KalRep* krep = KFRes.krep;

    if ((krep != 0) && KFRes.fit.success()){
      Mu2eDetectorModel const& detmodel{ art::ServiceHandle<BTrkHelper>()->detectorModel() };
      std::vector<TrkHit*>::iterator         ihigh;
      std::vector<TrkHit*>::reverse_iterator ilow;
//-----------------------------------------------------------------------------
// use the reference trajectory, as that's what all the existing hits do
//-----------------------------------------------------------------------------
      const TrkDifPieceTraj* reftraj = krep->referenceTraj();

      vector<TrkHit*>& hits = krep->hitVector();

      int debug_addhits = _debugLevel/100;

      if (debug_addhits > 0) {
        printf("%s  shId   A      sec   panel       ambig   res       hitRMS      drift      res/hitRMS\n",oname);
      }
      int                iambig;
      double             hflt, mom, beta, tflt;
      const TrkStrawHit* nearhit;
      TrkStrawHit*       trkhit;

      int nadd = KFRes.missingHits.size();
      for (int i=0; i<nadd; i++) {
        size_t istraw = KFRes.missingHits[i].index;
        const StrawHit& strawhit(KFRes.shcol->at(istraw));
        const Straw& straw = _tracker->getStraw(strawhit.strawId());
//-----------------------------------------------------------------------------
// estimate  initial flightlength
//-----------------------------------------------------------------------------
        hflt = 0;
        TrkHelixUtils::findZFltlen(*reftraj,straw.getMidPoint().z(),hflt);
//-----------------------------------------------------------------------------
// find the bounding sites near this hit, and extrapolate to get the hit t0
// hits are supposed to be sorted 
//-----------------------------------------------------------------------------
        std::sort(hits.begin(),hits.end(),fltlencomp(KFRes.fdir.fitDirection()));
        findBoundingHits(hits,hflt,ilow,ihigh);

        if(ihigh != hits.end()) nearhit = static_cast<TrkStrawHit*> (*ihigh);
        else                    nearhit = static_cast<TrkStrawHit*> (*ilow );

        TrkT0 hitt0 = nearhit->hitT0();

        mom  = krep->momentum(nearhit->fltLen()).mag();
        beta = krep->particleType().beta(mom);
        tflt = (hflt-nearhit->fltLen())/(beta*CLHEP::c_light);

                                        // update the time in the TrkT0 object and create a new
                                        // hit object.  Assume we're at the last iteration over added error
        hitt0._t0 += tflt;
        trkhit     = new TrkStrawHit(strawhit,straw,istraw,hitt0,hflt,hit_error,_maxdriftpull,1.,_mint0doca);
        assert(trkhit != 0);
                                        // 3 means "Added by addHits"
        trkhit->setFlag(3);
	if (KFRes.missingHits[i].doca > 0) iambig =  1;
	else                               iambig = -1;
					// can set ambiguity only for deactivated hit
	trkhit->setActivity(false);
	trkhit->setAmbig(iambig);
	trkhit->setActivity(true);
                                        // set by constructor, don't need to do it twice

	//        trkhit->setAmbigUpdate(false);

                                        // add the hit to the track and the fit, from now on track owns it
        krep->addHit(trkhit);
                                        // create intersections for the material of this hit
                                        // and add those to the track

        const DetStrawElem* strawelem = detmodel.strawElem(trkhit->straw());
        DetIntersection strawinter;
        strawinter.delem   = strawelem;
        strawinter.pathlen = trkhit->fltLen();
        strawinter.thit    = trkhit;
        if(strawelem->reIntersect(reftraj,strawinter)) {
          krep->addInter(strawinter);
	}
//-----------------------------------------------------------------------------
// hit is added with _iamb=0, which means that the position uncertainty
// includes the drift distance
// by design, findMissingHits looks for hits close to the trajectory - why trying to be smart here?
// what happens, if I turn this off?
//-----------------------------------------------------------------------------
        double dr     = trkhit->residual();
	double sig_dr = trkhit->hitRms();
        double chi    = dr/sig_dr;

        activity = 1;
// if it's outside limits, deactivate the HOT
        if ((fabs(chi) > MaxChi) || (!trkhit->isPhysical( MaxChi))) {
          trkhit->setActivity(false);
          activity = 0;
        }

        if (debug_addhits > 0) {
          printf("%s %5i %3i  %6i  %6i  %6i  %10.3f  %10.3f %10.3f  %10.3f \n",oname,
                 straw.index().asInt(),
                 activity,
                 straw.id().getPlane(),
                 straw.id().getPanel(),
                 trkhit->ambig(),
                 trkhit->residual(),
                 trkhit->hitRms(),
                 trkhit->driftRadius(),
                 chi);
        }
      }
//-----------------------------------------------------------------------------
// sort hits by flightlength (or in Z, which is the same)
//-----------------------------------------------------------------------------
      std::sort(hits.begin(),hits.end(),fltlencomp(KFRes.fdir.fitDirection()));
//---------------------------------------------------------------------------
// refit the track one more time with minimal external errors
//
// 2015 - 03 - 09 Gainipez added the following line for forcing the fiITeration procedure
// to use findAndUseDoublets
// 2015 - 02 - 27 Gianipez added the loop for including the external errors
// 2015-04-03 P.Murat: not sure I understand why there are different number of iterations
//                     in different cases - Giani?
// 2015-04-10        : perform one iteration, -1 means 'use the smallest external error defined'
//------------------------------------------------------------------------------------------
      fitIteration(KFRes, -1);

      KFRes.krep->addHistory(KFRes.fit,"AddHits");
    }
  }


//-----------------------------------------------------------------------------
// Dave's "annealing" procedure: refit the track 'KRes' multiple times 
//                               gradually reducing the external hit errors
//
// assume this is not the final fit, 
// the final fit will be performed later, after the hit pickup
//
// final=1 in fitIteration means 'make decision about the hit drift directions no matter what'
//      =0                 means 'reduce external hit error down to zero only for hits which 
//                                drift direction is well defined'
//-----------------------------------------------------------------------------
  void KalFitHackNew::fitTrack(KalFitResultNew& KRes) {
    //    int not_final(0);

    int n = _hiterr.size();
    for (int i=0; i<n; ++i) {
      fitIteration(KRes,i);

      if (! KRes.fit.success()) break; //commented by gianipez
     }

    if(KRes.krep != 0) KRes.krep->addHistory(KRes.fit,"KalFitHackNew");
  }

//-----------------------------------------------------------------------------
// one step of the track fit
// update external hit errors.  This isn't strictly necessary on the 1st iteration.
//-----------------------------------------------------------------------------
  void KalFitHackNew::fitIteration(KalFitResultNew& KRes, int Iteration) {

    // update the external hit errors.  This isn't strictly necessary on the 1st iteration.
    TrkHitVector* thv   = &(KRes.krep->hitVector());
    for (auto itsh=thv->begin();itsh!=thv->end(); ++itsh){
      (*itsh)->setTemperature(_hiterr[Iteration]);
    }

    unsigned     niter(0);
    bool         changed(true);
    bool         fit_success;
    char         msg[100];
    double       oldt0;

    if (_debugLevel >0) {
      printf("------------------------------------------------------------------------------------------\n");
      printf("[KalFitHackNew::fitIteration] BEGIN Iteration:%i \n", Iteration);
      printf("------------------------------------------------------------------------------------------\n");
    }
//-----------------------------------------------------------------------------------
// 2015 -02 -17 G. Pezzullo: loop over the hits and assign a smaller external error
// for the doublets
// 2015-04-03: IHErr = -1: special value, invoke last iteration .. 
// 2016-05-08 P.Murat: there are Niter+1 ambig resolvers, so last iteration = _hiterr.size()-1
//-----------------------------------------------------------------------------------
    if (Iteration == -1) Iteration = _hiterr.size()-1;

    _annealingStep     = Iteration;
//--------------------------------------------------------------------------------
// 2015-02-19 G.Pezzu: re-search multiplets using updated fit results
// 2015-03-25 P.Murat: *TODO* I don't think this call is needed, the one in the
//                     loop should be sufficient - check !
//-----------------------------------------------------------------------------
    KalRep* krep =  KRes.krep;
    KRes.nt0iter = 0;
    KRes.fit     = TrkErrCode::succeed;
    oldt0        = krep->t0()._t0;

    while (KRes.fit.success() && changed && niter < maxIterations()) {
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
      _ambigresolver[_annealingStep]->resolveTrk(krep);
//-----------------------------------------------------------------------------
// perform the track fit
//-----------------------------------------------------------------------------
      krep->resetFit();
      KRes.fit = krep->fit();

      fit_success = KRes.fit.success();

      if (_debugLevel > 0) {
        sprintf(msg,"KalFitHackNew::fitIteration::002 niter = %2i fit_success = %i",
                niter,fit_success);
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
        if (KRes.caloCluster != NULL) updateCalT0(KRes);
        else                          updateT0   (krep);
      }
      else if (_updateT0Mode == 1) {
//-----------------------------------------------------------------------------
// mode = 1: when iterating, don't look back at the cluster T0,
//           in this mode the cluster T0 is used only to seed the process
//-----------------------------------------------------------------------------
        updateT0(krep);
      }

//-----------------------------------------------------------------------------
// look at T0 changes only during the final fit
//-----------------------------------------------------------------------------
      if (KRes.fitType == 1) changed |= fabs(krep->t0()._t0-oldt0) > _t0tol[_annealingStep];
      oldt0    = krep->t0()._t0;
//-----------------------------------------------------------------------------
// drop outliers. weedHits() calls KalRep::fit(), so the fit success code may change
//-----------------------------------------------------------------------------
      if(_weedhits[Iteration]) {
        KRes.nweediter = 0;
        changed        |= weedHits(KRes,Iteration);
        fit_success     = KRes.fit.success();
      }
//-----------------------------------------------------------------------------
// find missing materials - latest borrowed from Dave
//-----------------------------------------------------------------------------
      if (_addmaterial[Iteration]) {
        changed |= addMaterial(krep) > 0;
      }

      niter++;
    }

    if (!krep->fitCurrent()) {
      KRes.fit = krep->fit();
    }
//-----------------------------------------------------------------------------
// done iterating, define drift signs with respect to the final trajectory
// 2015-02-17 G. Pezzu: ::resolveTrk() updates drift signs of ALL hits,
// so fix the ambiguity of the doublets after that
//-----------------------------------------------------------------------------
    _nIter += niter;
  }


//-----------------------------------------------------------------------------
// stolen clone of KalFit::addMaterials - need to work towards
//-----------------------------------------------------------------------------
  // this function belongs in TrkDifTraj, FIXME!!!!
  double KalFitHackNew::zFlight(KalRep* krep,double pz) {
// get the helix at the middle of the track
    double loclen;
    double fltlen(0.0);
    const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(krep->referenceTraj()->localTrajectory(fltlen,loclen));
// Iterate
    const HelixTraj* oldtraj;
    unsigned iter(0);
    do {
// remember old traj
      oldtraj = htraj;
// correct the global fltlen for this difference in local trajectory fltlen at this Z position
      fltlen += (htraj->zFlight(pz)-loclen);
      htraj = dynamic_cast<const HelixTraj*>(krep->referenceTraj()->localTrajectory(fltlen,loclen));
    } while(oldtraj != htraj && iter++<10);
    return fltlen;
  }

//-----------------------------------------------------------------------------
  unsigned KalFitHackNew::addMaterial(KalRep* KRep) {
    unsigned retval(0);
    // TTracker geometry
    const Tracker& tracker = getTrackerOrThrow();
    const TTracker& ttracker = dynamic_cast<const TTracker&>(tracker);
    // fetcth the DetectorModel
    Mu2eDetectorModel const& detmodel{ art::ServiceHandle<BTrkHelper>()->detectorModel() };
    // storage of potential straws
    StrawFlightComp strawcomp(_maxmatfltdiff);
    std::set<StrawFlight,StrawFlightComp> matstraws(strawcomp);
    // loop over Planes
    double strawradius = ttracker.strawRadius();
    unsigned nadded(0);
    // for(auto const& plane : ttracker.getPlanes()){
    for ( size_t i=0; i!= ttracker.nPlanes(); ++i){
      const auto& plane = ttracker.getPlane(i);
      int nstraws = plane.getPanel(0).nStraws();
      // get an approximate z position for this plane from the average position of the 1st and last straws
      Hep3Vector s0 = plane.getPanel(0).getStraw(0).getMidPoint();
      Hep3Vector sn = plane.getPanel(0).getStraw(nstraws-1).getMidPoint();
      double pz = 0.5*(s0.z() + sn.z());
      // find the transverse position at this z using the reference trajectory
      double flt = zFlight(KRep,pz);
      HepPoint pos = KRep->referenceTraj()->position(flt);
      Hep3Vector posv(pos.x(),pos.y(),pos.z());
      // see if this position is in the active region.  Double the straw radius to be generous
      double rho = posv.perp();
      double rmin = s0.perp()-2*strawradius;
      double rmax = sn.perp()+2*strawradius;
      if(rho > rmin && rho < rmax){
	// loop over panels
        for(auto const& panel : plane.getPanels()){
	  // get the straw direction for this panel
          Hep3Vector sdir = panel.getStraw(0).getDirection();
	  // get the transverse direction to this and z
          static Hep3Vector zdir(0,0,1.0);
          Hep3Vector pdir = sdir.cross(zdir);
     //  project the position along this
          double prho = posv.dot(pdir);
      // test for acceptance of this panel
          if(prho > rmin && prho < rmax) {
          // translate the transverse position into a rough straw number
            int istraw = (int)rint(nstraws*(prho-s0.perp())/(sn.perp()-s0.perp()));
            // take a few straws around this
            for(int is = max(0,istraw-3); is<min(nstraws,istraw+3); ++is){
              matstraws.insert(StrawFlight(panel.getStraw(is).index(),flt));
              ++nadded;
            }
          }
        }
      }
    }
// Now test if the Kalman rep hits these straws
    if(_debugLevel > 2) std::cout << "Found " << matstraws.size() << " unique possible straws " << " out of " << nadded << std::endl;
    for(auto strawflt : matstraws){
    // hack
      Straw const& straw = tracker.getStraw(strawflt._index);
      const DetStrawElem* strawelem = detmodel.strawElem(straw.id());
      DetIntersection strawinter;
      strawinter.delem = strawelem;
      strawinter.pathlen = strawflt._flt;
      if(strawelem->reIntersect(KRep->referenceTraj(),strawinter)){
// If the rep already has a material site for this element, skip it
        std::vector<const KalMaterial*> kmats;
        KRep->findMaterialSites(strawelem,kmats);
        if(_debugLevel > 2) std::cout << "found intersection with straw " << strawelem->straw()->index() << " with "
        << kmats.size() << " materials " << std::endl;
// test material isn't on the track
        bool hasmat(false);
        for(auto kmat : kmats ){
          const DetStrawElem* kelem = dynamic_cast<const DetStrawElem*>(kmat->detIntersection().delem);
          if(kelem != 0){
            StrawFlight ksflt(kelem->straw()->index(),kmat->globalLength());
            if(_debugLevel>2)std::cout << " comparing flights " << kmat->globalLength() << " and " << strawflt._flt << std::endl;
            if(!strawcomp.operator()(strawflt,ksflt)){
              if(_debugLevel>2)std::cout << "operator returned false!!" << std::endl;
              // this straw is already on the track: stop
              hasmat = true;
              break;
            }
          }
        }
        if (kmats.size() == 0 || !hasmat) {
          if(_debugLevel > 2) std::cout << "Adding material element" << std::endl;
          // this straw doesn't have an entry in the Kalman fit: add it`
          DetIntersection detinter(strawelem, KRep->referenceTraj(),strawflt._flt);
          KRep->addInter(detinter);
          ++retval;
        }
      }
    }
    if (_debugLevel > 1) std::cout << "Added " << retval << " new material sites" << std::endl;
    return retval;
  }

//-----------------------------------------------------------------------------
// create list of track candidate hits, sort hits in flightlength
// this should be the only place where sorting is done
//-----------------------------------------------------------------------------
  void KalFitHackNew::makeHits(KalFitResultNew& KRes, TrkHitVector& ListOfHits) {

    TrkT0 hitt0; 

    const HelixTraj* hel = KRes.helixTraj;
    double     flt0 = hel->zFlight(0.0);
    double     mom  = TrkMomCalculator::vecMom(*hel,bField(),flt0).mag();
    double     vflt = KRes.tpart.beta(mom)*CLHEP::c_light;

    unsigned   nind = KRes.strawHitIndices()->size();

    for (unsigned iind=0; iind<nind; iind++) {
      size_t istraw = KRes.strawHitIndices()->at(iind);             //[iind];
      const StrawHit& strawhit(KRes.shcol->at(istraw));
      const Straw& straw = _tracker->getStraw(strawhit.strawId());
      double fltlen      = hel->zFlight(straw.getMidPoint().z());
    // estimate arrival time at the wire
      hitt0._t0          = KRes.t0._t0 + (fltlen-flt0)/vflt;
      // create the hit object.  Start with the 1st additional error for anealing
      TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,hitt0,fltlen,_hiterr.front(),_maxdriftpull,1.,_mint0doca);
      assert(trkhit != 0);
					// set hit flag and the initial ambiguity to null
      trkhit->setFlag(0);
      trkhit->setAmbig(0);
    // refine the flightlength, as otherwise hits in the same plane are at exactly the same flt, which can cause problems
      TrkErrCode pstat = trkhit->updatePoca(hel);
      if (pstat.failure()) {
        trkhit->setActivity(false);
      }
      ListOfHits.push_back(trkhit);
    }
					// sort the hits by flightlength
    std::sort(ListOfHits.begin(),ListOfHits.end(),fltlencomp(KRes.fdir.fitDirection()));
  }

//-----------------------------------------------------------------------------
// '*' in front of the hit drift radius: the hit drift sign is not defined
//     and the drift ambiguity has been set to 0
// '?': the drift sign determined by the resolver is different from the MC truth
//-----------------------------------------------------------------------------
  void KalFitHackNew::printHits(KalFitResultNew& KRes, const char* Caller) {
    const KalRep* Trk  = KRes.krep;

    if (Trk == NULL)  return;

    printf("[KalFitHackNew::printHits] BEGIN called from %s _annealingStep:%i \n",Caller,_annealingStep);
    printf("---------------------------------------------------------------------------------");
    printf("-----------------------------------------------------\n");
      //      printf("%s",Prefix);
    printf("  TrkID       Address    N  NA      P       pT     costh    T0      T0Err   Omega");
    printf("      D0       Z0      Phi0   TanDip    Chi2    FCons\n");
    printf("---------------------------------------------------------------------------------");
    printf("-----------------------------------------------------\n");

    Hep3Vector trk_mom;

    double h1_fltlen = Trk->firstHit()->kalHit()->hit()->fltLen() - 10;
    trk_mom          = Trk->momentum(h1_fltlen);
    double mom       = trk_mom.mag();
    double pt        = trk_mom.perp();
    double costh     = trk_mom.cosTheta();
    double chi2      = Trk->chisq();
    int    nhits     = Trk->hitVector().size();
    int    nact      = Trk->nActive();
    double t0        = Trk->t0().t0();
    double t0err     = Trk->t0().t0Err();
//-----------------------------------------------------------------------------
// in all cases define momentum at lowest Z - ideally, at the tracker entrance
//-----------------------------------------------------------------------------
    double s1     = Trk->firstHit()->kalHit()->hit()->fltLen();
    double s2     = Trk->lastHit ()->kalHit()->hit()->fltLen();
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

    printf(" %8.5f %8.3f %8.3f %8.4f %7.4f",omega,d0,z0,phi0,tandip
	   );
    printf(" %8.3f %10.3e\n",chi2,fit_consistency);
//-----------------------------------------------------------------------------
// print detailed information about the track hits
//-----------------------------------------------------------------------------
    printf("--------------------------------------------------------------------");
    printf("----------------------------------------------------------------");
    printf("--------------------------------------------\n");
    printf(" ih  SInd Flag      A     len         x        y        z      HitT     HitDt");
    printf(" Pl Pn  L  W     T0       Xs      Ys        Zs      resid  sigres");
    printf("   Rdrift   mcdoca  totErr hitErr  t0Err penErr extErr\n");
    printf("--------------------------------------------------------------------");
    printf("----------------------------------------------------------------");
    printf("--------------------------------------------\n");

    mu2e::TrkStrawHit     *hit;
    Hep3Vector            pos;
    const mu2e::StrawHit  *sh;
    const mu2e::Straw     *straw;
    int                   ihit;
    double                len;
    HepPoint              plen;

    ihit = 0;
    for (int it=0; it<nhits; ++it) {
      hit   = static_cast<TrkStrawHit*> (KRes.krep->hitVector().at(it));
      sh    = &hit->strawHit();
      straw = &hit->straw();

      hit->hitPosition(pos);

      len   = hit->fltLen();
      plen  = Trk->position(len);

      double mcdoca = _mcUtils->mcDoca(KRes.event,KRes.shDigiLabel.data(),straw);

      ihit += 1;
      printf("%3i %5i 0x%08x %1i %9.3f %8.3f %8.3f %9.3f %8.3f %7.3f",
             ihit,
             straw->index().asInt(),
             hit->hitFlag(),
             hit->isActive(),
             len,
             //      hit->hitRms(),
             plen.x(),plen.y(),plen.z(),
             sh->time(), sh->dt()
             );

      printf(" %2i %2i %2i %2i",
             straw->id().getPlane(),
             straw->id().getPanel(),
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
             hit->temperature()*hit->driftVelocity()
             );
    }
  }


//-----------------------------------------------------------------------------
  void KalFitHackNew::makeMaterials(TrkHitVector& ListOfHits, const HelixTraj& Helix, vector<DetIntersection>& ListOfIntersections) {
  // fetcth the DetectorModel
    DetIntersection strawinter;

    Mu2eDetectorModel const& detmodel{ art::ServiceHandle<BTrkHelper>()->detectorModel() };

    for(std::vector<TrkHit*>::iterator ihit=ListOfHits.begin();ihit!=ListOfHits.end();ihit++){
      TrkStrawHit* trkhit = static_cast<TrkStrawHit*>(*ihit);
      const DetStrawElem* strawelem = detmodel.strawElem(trkhit->straw());
      strawinter.delem   = strawelem;
      strawinter.pathlen = trkhit->fltLen();
      strawinter.thit    = trkhit;
      if (strawelem->reIntersect(&Helix,strawinter)) {
        ListOfIntersections.push_back(strawinter);
      }
    }
  }

//-----------------------------------------------------------------------------
  bool KalFitHackNew::weedHits(KalFitResultNew& KRes, int Iteration) {
    // Loop over HoTs and find HoT with largest contribution to chi2.  If this value
    // is greater than some cut value, deactivate that HoT and reFit

    bool         retval(false);
    double       value, worst(-1.);
    TrkStrawHit* worstHot(NULL);

    KalRep* krep = KRes.krep;

    TrkHitVector& hits = krep->hitVector();

    for (std::vector<TrkHit*>::iterator iter = hits.begin(); iter != hits.end(); ++iter){
      TrkStrawHit* iHot = static_cast<TrkStrawHit*> (*iter);
      if (iHot->isActive()) {
        double resid, residErr;
        if (iHot->resid(resid, residErr, true)) {
	  if (KRes.fitType == 0) {
//-----------------------------------------------------------------------------
// seed fit, not drift sign ambiguities
//-----------------------------------------------------------------------------
	    value = fabs(resid/residErr);
	    if (value > _maxhitchi && value > worst) {
	      worst    = value;
	      worstHot = iHot;
	    }
	  }
	  else if (KRes.fitType == 1) {
//-----------------------------------------------------------------------------
// final fit, logic depends on the drift ambiguity
// attempt to consider separately cases when the drift sign is defined and is not
//-----------------------------------------------------------------------------
	    if (iHot->ambig() != 0) {
	      value = fabs(resid/residErr);
	      if (value > _maxhitchi && value > worst) {
		worst    = value;
		worstHot = iHot;
	      }
	    }
	    else {
//-----------------------------------------------------------------------------
// ambig = 0, resid is the distance to the wire .... compare 2 radii
// if the hit is OK, one of the radii should be giving the right solution
// however, what is the error to use here? - the fit was using large error
//-----------------------------------------------------------------------------
	      double rdrift = iHot->driftRadius();
	      double d1     = fabs(resid+rdrift);
	      double d2     = fabs(resid-rdrift);

	      double dr = (d1 <= d2) ? d1 : d2;

	      double err = iHot->totalErr();
	      if (err > 1.5) err = 1.5 ;   // D= 5 mm , uniformly distributed

	      value = fabs(dr/err);
	      if (value > _maxhitchi && value > worst) {
		worst    = value;
		worstHot = iHot;
	      }
	    }
	  }
	}
      }
    }

    int iter;

    if (Iteration == -1) iter = _ambigresolver.size()-1;
    else                 iter = Iteration;

    if (worstHot) {
      retval = true;
      worstHot->setActivity(false);
      worstHot->setFlag(5);             // positive usability allows hot to be re-enabled later

      _ambigresolver[iter]->resolveTrk(krep);

      KRes.fit = krep->fit();
      krep->addHistory(KRes.fit, "HitWeed");

      if (_debugLevel > 0) {
	char msg[200];
        sprintf(msg,"KalFitHackNew::weedHits Iteration = %2i success = %i",Iteration,KRes.fit.success());
        printHits(KRes,msg);
      }
//-----------------------------------------------------------------------------
// Recursively iterate
//-----------------------------------------------------------------------------
      KRes.nweediter++;
      if (KRes.fit.success() && KRes.nweediter < _maxweed ) {
        retval |= weedHits(KRes,iter);
      }
    }
    return retval;
  }


//-----------------------------------------------------------------------------
// assume that this is a call during the last iteration
//-----------------------------------------------------------------------------
  bool KalFitHackNew::unweedHits(KalFitResultNew& KRes, double maxchi) {
    // Loop over inactive HoTs and find the one with the smallest contribution to chi2.  If this value
    // is less than some cut value, reactivate that HoT and reFit
    bool       retval(false);
//    int  const final (1);
    double     best = 1.e12;
    int        last_iteration = _hiterr.size()-1;

    KalRep* krep         = KRes.krep;
    TrkHitVector& hits   = krep->hitVector();
    TrkStrawHit* bestHot = 0;

    for (std::vector<TrkHit*>::iterator iter = hits.begin(); iter != hits.end(); ++iter){
      TrkStrawHit* iHot = static_cast<TrkStrawHit*> (*iter);
      if (!iHot->isActive()) {
        double resid, residErr;
        if(iHot->resid(resid, residErr, true)){
          double chival = fabs(resid/residErr);
  // test both for a good chisquared and for the drift radius to be physical
          if (chival < maxchi && iHot->isPhysical(maxchi) && chival < best) {
            best    = chival;
            bestHot = iHot;
          }
        }
      }
    }
    if (0 != bestHot) {
      retval = true;
      bestHot->setActivity(true);
      bestHot->setFlag(4);
//-----------------------------------------------------------------------------
// update drift signs before the fit again - more doublets could've been recovered
// one more place for the ambiguity resolver to be called from
// also - from addHits
// here we call fitter not invoking the fitIteration
//-----------------------------------------------------------------------------
      _ambigresolver[last_iteration]->resolveTrk(krep);

      KRes.fit = krep->fit();
      krep->addHistory(KRes.fit, "HitUnWeed");

      if (_debugLevel > 0) {
	char msg[200];
        sprintf(msg,"KalFitHackNew::unweedHits Iteration = %2i success = %i",last_iteration,KRes.fit.success());
        printHits(KRes,msg);
      }
      // Recursively iterate
      KRes.nunweediter++;
      if (KRes.fit.success() && bestHot->isActive()) {
        retval |= unweedHits(KRes,maxchi);
      }
    }
    return retval;
  }


  BField const&  KalFitHackNew::bField() const {
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

  //-----------------------------------------------------------------------------
  const TrkVolume* KalFitHackNew::trkVolume(trkDirection trkdir) const {
    //FIXME!!!!
    return 0;
  }

//----------------------------------------------------------------------------
// if a calorimeter cluster is present, estimate T0 using the cluster time and 
// the existing parameterization of the trajectory
// in case there is no cluster, estimate T0 using the hit times and the trajectory 
// parameterization and assuming that the drift time is defined by R(straw)/2
//-----------------------------------------------------------------------------
  void KalFitHackNew::initT0(KalFitResultNew& KRes) {
    using namespace boost::accumulators;

    static CLHEP::Hep3Vector zdir(0.0,0.0,1.0);
    const HelixTraj* hel = KRes.helixTraj;
    double t0flt         = hel->zFlight(0.0);
    double mom           = TrkMomCalculator::vecMom(*hel,bField(),t0flt).mag();
    double vflt          = KRes.tpart.beta(mom)*CLHEP::c_light;

    const CaloCluster* cl = KRes.caloCluster;
    if (cl) {
//-----------------------------------------------------------------------------
// calculate the path length of the particle from the middle of the Tracker to the
// calorimeter, cl->Z() is calculated wrt the tracker center
// _dtoffset : global time offset between the tracker and the calorimeter, 
//             think of an average cable delay
//-----------------------------------------------------------------------------
      Hep3Vector gpos = _calorimeter->geomUtil().diskToMu2e(cl->diskId(),cl->cog3Vector());
      Hep3Vector tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);
      double path     = tpos.z()/hel->sinDip();
      KRes.t0._t0    = cl->time() + _dtoffset - path/vflt;
      KRes.t0._t0err = 1.;
    }
    else {
//-----------------------------------------------------------------------------
// cluster is not present
//-----------------------------------------------------------------------------
      int nind = KRes.strawHitIndices()->size();
      std::vector<double> times(nind);
					// for crude estimates, we only need 1 d2t function
      D2T d2t;
      for (int iind=0; iind<nind; iind++) {

	size_t istraw = KRes.strawHitIndices()->at(iind);
	const StrawHit& strawhit(KRes.shcol->at(istraw));
	const Straw& straw = _tracker->getStraw(strawhit.strawId());
//-----------------------------------------------------------------------------
// compute the flightlength to this hit from z=0 (can be negative)
// and use this to estimate the time for the track to reach this hit from z=0
//-----------------------------------------------------------------------------
	double hflt  = hel->zFlight(straw.getMidPoint().z()) - t0flt;
	double tprop = hflt/vflt;
	// estimate signal propagation time on the wire assuming the middle (average)
	double vwire = _tcal->SignalVelocity(straw.index());
	double teprop = straw.getHalfLength()/vwire;
	// correct the measured time for these effects: this gives the aveage time the particle passed this straw, WRT
	// when the track crossed Z=0
	// assume the average drift time is half the maximum drift distance.  This is a poor approximation, but good enough for now
	if(iind==0)_tcal->DistanceToTime(straw.index(),0.5*straw.getRadius(),zdir,d2t);
	double htime = strawhit.time() - tprop - teprop - d2t._tdrift;
	times.push_back(htime);
      }
      // find the median time
      accumulator_set<double, stats<tag::median(with_p_square_quantile)>> med;
      med = std::for_each( times.begin(), times.end(), med );
      KRes.t0._t0 = extract_result<tag::median>(med);

      accumulator_set<double, stats<tag::min>> min;
      accumulator_set<double, stats<tag::max>> max;
      min = std::for_each( times.begin(), times.end(), min);
      max = std::for_each( times.begin(), times.end(), max);
      double tmin = extract_result<tag::min>(min);
      double tmax = extract_result<tag::max>(max);
					// estimate the error using the range
      KRes.t0._t0err = (tmax-tmin)/sqrt(12*nind);
    }
  }


//-----------------------------------------------------------------------------
// track T0 is the time when the particle crossed the Z=0 plane
// 2014-11-24 gianipez and Pasha removed time offset between caloriemter and tracker
//-----------------------------------------------------------------------------
  void KalFitHackNew::updateCalT0(KalFitResultNew& KRes) {

    TrkT0 t0;
    double mom, vflt, path, t0flt, flt0(0.0);

    KalRep* krep   = KRes.krep;
    bool converged = TrkHelixUtils::findZFltlen(krep->traj(),0.0,flt0);

    // get helix from kalrep

    const CaloCluster* cl = KRes.caloCluster;

    HelixTraj trkHel(krep->helix(flt0).params(),krep->helix(flt0).covariance());

                                        // get flight distance of z=0
    t0flt = trkHel.zFlight(0.0);

    if (converged) {
//-----------------------------------------------------------------------------
// estimate the momentum at that point using the helix parameters.
// This is assumed constant for this crude estimate
// compute the particle velocity
//-----------------------------------------------------------------------------
      mom  = TrkMomCalculator::vecMom(trkHel,bField(),t0flt).mag();
      vflt = krep->particleType().beta(mom)*CLHEP::c_light;
//-----------------------------------------------------------------------------
// path length of the particle from the middle of the Tracker to the  calorimeter
// set dummy error value
//-----------------------------------------------------------------------------
      Hep3Vector gpos = _calorimeter->geomUtil().diskToMu2e(cl->diskId(),cl->cog3Vector());
      Hep3Vector tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);

      path      = tpos.z()/trkHel.sinDip();
      t0._t0    = cl->time() + _dtoffset - path/vflt;
      t0._t0err = 1.;

      krep->setT0(t0,flt0);
      updateHitTimes(krep);
    }
  }


//-----------------------------------------------------------------------------
  bool KalFitHackNew::updateT0(KalRep* KRep) {
    using namespace boost::accumulators;
    //    bool retval(false);

    // TrkHitVector& hits = KRep->hitVector();

    // int nhits   = hits.size();
    // 					// need to have a valid fit
    // if (KRep->fitValid()) {
    // 					// find the global fltlen associated with z=0.
    //   double flt0(0.0);
    //   bool   converged = TrkHelixUtils::findZFltlen(KRep->traj(),0.0,flt0);
    //   if (converged) {
    //     std::vector<double> hitt0   (nhits);          // store t0, to allow outlyer removal
    //     std::vector<double> hitt0err(nhits);
    //     // loop over the hits
    //     for(std::vector<TrkHit*>::iterator ihit= hits.begin();ihit != hits.end(); ihit++) {
    //       TrkStrawHit* hit = static_cast<TrkStrawHit*>(*ihit);
    //       if(hit->isActive() && hit->hasResidual()){
    //         // find the residual, exluding this hits measurement
    //         double resid,residerr;
    //         if(KRep->resid(hit,resid,residerr,true)){
    //           // convert this to a distance to the wire
    //           double doca = (resid + hit->driftRadius()*hit->ambig());
    //           if(hit->ambig() == 0)
    //             doca = fabs(doca);
    //           else
    //             doca *= hit->ambig();
    //           // restrict the range, symmetrically to avoid bias
    //           double rad = hit->straw().getRadius();
    //           if(doca > _mint0doca && doca < rad-_mint0doca){
    //             // translate the DOCA into a time
    //             D2T d2t;
    //             _tcal->DistanceToTime(hit->straw().index(),doca,KRep->traj().direction(hit->fltLen()),d2t);
    //             // subtracting hitT0 makes this WRT the previous track t0
    //             hitt0.push_back(hit->time() - d2t._tdrift - hit->signalTime() - hit->hitT0()._t0);
    //             // assume residual error dominates
    //             hitt0err.push_back(residerr/d2t._vdrift);
    //           }
    //         }
    //       }
    //     }
    //     if (hitt0.size() > 1) {
    //       TrkT0 t0;
    //       // find the median
    //       accumulator_set<double, stats<tag::median(with_p_square_quantile) > > med;
    //       med = std::for_each( hitt0.begin(), hitt0.end(), med );
    //       t0._t0 = extract_result<tag::median>(med);
    //       // iterate an outlier search and linear fit until the set of used hits doesn't change
    //       bool changed(true);
    //       std::vector<bool> used(hitt0.size(),true);
    //       unsigned niter(0);
    //       while(changed && niter < 10){
    //         niter++;
    //         changed = false;
    //         accumulator_set<double,stats<tag::weighted_variance>,double> wmean;
    //         for (unsigned ihit=0; ihit<hitt0.size(); ihit++) {
    //           bool useit = fabs(hitt0[ihit]-t0._t0) < _t0nsig*hitt0err[ihit];
    //           changed   |= useit != used[ihit];
    //           used[ihit] = useit;
    //           if(useit){
    //             wmean(hitt0[ihit], weight=1.0/(hitt0err[ihit]*hitt0err[ihit]));
    //           }
    //         }
    //         unsigned nused = extract_result<tag::count>(wmean);
    //         if (nused > 1) {
    //           t0._t0    = extract_result<tag::weighted_mean>(wmean);
    //           t0._t0err = sqrt(extract_result<tag::weighted_variance>(wmean)/nused);
    //         } 
    // 	    else {
    //           break;
    //         }
    //       }

    TrkHitVector *thv = &(KRep->hitVector());
    bool retval(false);
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
// need to have a valid fit
    if(KRep->fitValid()){
// find the global fltlen associated with z=0.
      double flt0(0.0);
      bool converged = TrkHelixUtils::findZFltlen(KRep->traj(),0.0,flt0);
      if(converged){
        std::vector<double> hitt0; // store t0, to allow outlyer removal
        std::vector<double> hitt0err;
        size_t nhits = KRep->hitVector().size();
        hitt0.reserve(nhits);
        hitt0err.reserve(nhits);
        // loop over the hits
        for(auto ihit=thv->begin(); ihit != thv->end(); ihit++){
          TrkStrawHit* hit = static_cast<TrkStrawHit*> (*ihit);
          if(hit->isActive() && hit->hasResidual()){
            // find the residual, exluding this hits measurement
            double resid,residerr;
            if(KRep->resid(hit,resid,residerr,true)){
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
                _tcal->DistanceToTime(hit->straw().index(),doca,KRep->traj().direction(hit->fltLen()),d2t);
                // subtracting hitT0 makes this WRT the previous track t0
                hitt0.push_back(hit->time() - d2t._tdrift - hit->signalTime() - hit->hitT0()._t0);
                // assume residual error dominates
                hitt0err.push_back(residerr/d2t._vdrift);
              }
            }
          }
        }
        if (hitt0.size() > 1) {
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
            accumulator_set<double,stats<tag::weighted_variance>,double> wmean;
            for (unsigned ihit=0; ihit<hitt0.size(); ihit++) {
              bool useit = fabs(hitt0[ihit]-t0._t0) < _t0nsig*hitt0err[ihit];
              changed   |= useit != used[ihit];
              used[ihit] = useit;
              if(useit){
                wmean(hitt0[ihit], weight=1.0/(hitt0err[ihit]*hitt0err[ihit]));
              }
            }
            unsigned nused = extract_result<tag::count>(wmean);
            if (nused > 1) {
              t0._t0    = extract_result<tag::weighted_mean>(wmean);
              t0._t0err = sqrt(extract_result<tag::weighted_variance>(wmean)/nused);
            } 
	    else {
              break;
            }
          }
          // reset t0
          if (!changed) {
            // put in t0 from the track.
            t0._t0 += KRep->t0()._t0;
            KRep->setT0(t0,flt0);
            updateHitTimes(KRep);
            retval = true;
          }
        }
      }
    }
    return retval;
  }

//-----------------------------------------------------------------------------
  void KalFitHackNew::updateHitTimes(KalRep* KRep) {

  // compute the time the track came closest to the wire for each hit, starting from t0 and working out.
  // this function allows for momentum change along the track.
  // find the bounding hits on either side of this

    TrkHitVector& hits = KRep->hitVector();

    std::vector<TrkHit*>::iterator         ihigh;
    std::vector<TrkHit*>::reverse_iterator ilow ;

    findBoundingHits(hits,KRep->flt0(),ilow,ihigh);

					// reset all the hit times
    double hflt = KRep->flt0();
    TrkT0 hitt0 = KRep->t0();
    for(std::vector<TrkHit*>::iterator ihit= ihigh;ihit != hits.end(); ++ihit) {
      TrkStrawHit* hit = static_cast<TrkStrawHit*> (*ihit);
// particle momentum at this point, using the full fit
      double mom = KRep->momentum(hit->fltLen()).mag();
// relativistic velocity from that
      double beta = KRep->particleType().beta(mom);
// particle transit time to this hit from the reference
      double tflt = (hit->fltLen()-hflt)/(beta*CLHEP::c_light);
// update the time in the TrkT0 object
      hitt0._t0 += tflt;
      hit->setHitT0(hitt0);
// update the reference flightlength
      hflt = hit->fltLen();
    }

// now the same, moving backwards
    hflt  = KRep->flt0();
    hitt0 = KRep->t0();
    for(std::vector<TrkHit*>::reverse_iterator ihit= ilow; ihit != hits.rend(); ++ihit){
      TrkStrawHit* hit = static_cast<TrkStrawHit*> (*ihit);
      double mom  = KRep->momentum(hit->fltLen()).mag();
      double beta = KRep->particleType().beta(mom);
      double tflt = (hit->fltLen()-hflt)/(beta*CLHEP::c_light);
      hitt0._t0  += tflt;
      hit->setHitT0(hitt0);
      hflt = hit->fltLen();
    }
  }

//-----------------------------------------------------------------------------
  void KalFitHackNew::findBoundingHits(vector<TrkHit*>&                   hits,
				       double                                  flt0,
				       vector<TrkHit*>::reverse_iterator& ilow,
				       vector<TrkHit*>::iterator&         ihigh) {
    ilow  = hits.rbegin();
    ihigh = hits.begin ();
    while (ilow  != hits.rend() && (*ilow )->fltLen() > flt0) ++ilow ;
    while (ihigh != hits.end () && (*ihigh)->fltLen() < flt0) ++ihigh;
  }

}



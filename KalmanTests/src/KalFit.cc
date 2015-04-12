//
// Class to perform BaBar Kalman fit
//
// $Id: KalFit.cc,v 1.43 2014/08/22 16:10:41 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2014/08/22 16:10:41 $
//

// the following has to come before other BaBar includes
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/KalFit.hh"
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
#include "KalmanTrack/KalBend.hh"
#include "TrkBase/HelixTraj.hh"
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
  KalFit::KalFit(fhicl::ParameterSet const& pset) :
// KalFit parameters
    _debug(pset.get<int>("debugLevel",0)),
    _weedhits(pset.get<bool>("weedhits",true)),
    _maxhitchi(pset.get<double>("maxhitchi",4.0)),
    _maxweed(pset.get<unsigned>("maxweed",10)),
    _herr(pset.get< vector<double> >("hiterr")),
    _maxdriftpull(pset.get<double>("maxDriftPull",10)),
    // t0 parameters
    _initt0(pset.get<bool>("initT0",true)),
    _updatet0(pset.get<bool>("updateT0",true)),
    _t0tol(pset.get< vector<double> >("t0Tolerance")),
    _t0errfac(pset.get<double>("t0ErrorFactor",1.2)),
    _mint0doca(pset.get<double>("minT0DOCA",-0.2)),
    _t0nsig(pset.get<double>("t0window",2.5)),
    //
    _removefailed(pset.get<bool>("RemoveFailedFits",true)),
    _minnstraws(pset.get<unsigned>("minnstraws",15)),
    _ambigstrategy(pset.get< vector<int> >("ambiguityStrategy")),
    _bfield(0)
  {
    // 2015-04-12 P.Murat add doublet ambig resolver
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
    _divergeflt = pset.get<double>("DivergeFlt",1.0e3);
    _mindot = pset.get<double>("MinDot",0.0);
    _maxmomdiff = pset.get<double>("MaxMomDiff",0.5);
    _momfac = pset.get<double>("MomFactor",0.0);
    _maxpardif[0] = _maxpardif[1] = pset.get<double>("MaxParameterDifference",1.0);

    // DOF counting subdivision is illogical, FIXME!!!!
    _mindof[0] = _mindof[2] = pset.get<double>("MinNDOF",10);
    _mindof[1] = 0;
    _bintconfig._maxRange = pset.get<double>("BFieldIntMaxRange",1.0e5); // 100 m
    _bintconfig._intTolerance = pset.get<double>("BFieldIntTol",0.01); // 10 KeV
    _bintconfig._intPathMin = pset.get<double>("BFieldIntMin",20.0); // 20 mm
    _bintconfig._divTolerance = pset.get<double>("BFieldIntDivTol",0.05); // 50 KeV
    _bintconfig._divPathMin = pset.get<double>("BFieldIntDivMin",50.0); // 50 mm
    _bintconfig._divStepCeiling = pset.get<double>("BFieldIntDivMax",500.0); // 500 mm
    // field integral errors
    double perr = pset.get<double>("BendCorrErrFrac",0.0); // fractional accuracy of trajectory
    double berr = pset.get<double>("BFieldMapErr",0.0); // mapping and interpolation error
//    KalBend::setErrors(perr,berr);
    // make sure we have at least one entry for additional errors
    if(_herr.size() <= 0) throw cet::exception("RECO")<<"mu2e::KalFit: no hit errors specified" << endl;
    if(_herr.size() != _ambigstrategy.size()) throw cet::exception("RECO")<<"mu2e::KalFit: inconsistent ambiguity resolution" << endl;
    if(_herr.size() != _t0tol.size()) throw cet::exception("RECO")<<"mu2e::KalFit: inconsistent ambiguity resolution" << endl;
    // construct the ambiguity resolvers

    AmbigResolver* ar;

    _useDoublets = 0;

    for(size_t i=0; i<_ambigstrategy.size(); ++i){
      switch (_ambigstrategy[i]) {
      case fixedambig: default:
	ar = new FixedAmbigResolver(pset);
	break;
      case hitambig:
	ar = new HitAmbigResolver(pset);
	break;
      case panelambig:
	ar = new PanelAmbigResolver(pset);
	break;
      case pocaambig:
	ar = new PocaAmbigResolver(pset);
	break;
      case doubletambig: // 4
 	ar = new DoubletAmbigResolver(*_darPset,i);
	_useDoublets = 1;
 	break;
      }

      ar->setExterr(_herr[i]);
      _ambigresolver.push_back(ar);
    }

//     for(size_t iambig=0;iambig<_ambigstrategy.size();++iambig){
//       switch (_ambigstrategy[iambig] ){
// 	case fixedambig: default:
// 	  _ambigresolver.push_back(new FixedAmbigResolver(pset));
// 	  break;
// 	case hitambig:
// 	  _ambigresolver.push_back(new HitAmbigResolver(pset));
// 	  break;
// 	case panelambig:
// 	  _ambigresolver.push_back(new PanelAmbigResolver(pset));
// 	  break;
// 	case pocaambig:
// 	  _ambigresolver.push_back(new PocaAmbigResolver(pset));
// 	  break;
//       }
//     }
  }

  KalFit::~KalFit(){
    for(size_t iambig=0;iambig<_ambigresolver.size();++iambig){
      delete _ambigresolver[iambig];
    }
    delete _bfield;
  }

  void KalFit::makeTrack(KalFitResult& kres) {
    kres._fit = TrkErrCode(TrkErrCode::fail);
// test if fitable
    if(fitable(kres._tdef)){
// first, find t0
      TrkT0 t0;
      if(_initt0)
	initT0(kres._tdef, t0);
      else
	t0 = kres._tdef.t0();
// create the hits
      makeHits(kres, t0);
// Create the BaBar hit list, and fill it with these hits.  The BaBar list takes ownership
// This will go away when we cleanup the BaBar hit storage, FIXME!!!
      TrkHotListFull* hotlist = new TrkHotListFull();
      for(std::vector<TrkStrawHit*>::iterator ihit=kres._hits.begin();ihit!=kres._hits.end();ihit++){
        TrkStrawHit* trkhit = *ihit;
	hotlist->append(trkhit);
	if (_debug>0) { trkhit->print(std::cout); }
      }
// Find the wall and gas material description objects for these hits
      if(_matcorr)makeMaterials(kres);
// create Kalman rep
      kres._krep = new KalRep(kres._tdef.helix(), hotlist, kres._detinter, *this, kres._tdef.particle());
      assert(kres._krep != 0);
// initialize krep t0; eventually, this should be in the constructor, FIXME!!!
      double flt0 = kres._tdef.helix().zFlight(0.0);
      kres._krep->setT0(t0,flt0);
// now fit
      fitTrack(kres);
      if(_removefailed)kres.removeFailed();
    }
  }

  void KalFit::addHits(KalFitResult& kres,const StrawHitCollection* straws, std::vector<hitIndex> indices, double maxchi) {
// there must be a valid Kalman fit to add hits to
    if(kres._krep != 0 && kres._fit.success()){
      ConditionsHandle<TrackerCalibrations> tcal("ignored");
      const Tracker& tracker = getTrackerOrThrow();
      std::vector<TrkStrawHit*>::iterator ihigh;
      std::vector<TrkStrawHit*>::reverse_iterator ilow;
// use the reference trajectory, as that's what all the existing hits do
      const TrkDifPieceTraj* reftraj = kres._krep->referenceTraj();
      for(unsigned iind=0;iind<indices.size(); ++iind){
	size_t istraw = indices[iind]._index;
	const StrawHit& strawhit(straws->at(istraw));
	const Straw& straw = tracker.getStraw(strawhit.strawIndex());
// estimate  initial flightlength
	double hflt(0.0);
	TrkHelixUtils::findZFltlen(*reftraj,straw.getMidPoint().z(),hflt);
// find the bounding sites near this hit, and extrapolate to get the hit t0
	std::sort(kres._hits.begin(),kres._hits.end(),fltlencomp(kres._tdef.fitdir().fitDirection()));
	findBoundingHits(kres._hits,hflt,ilow,ihigh);
	const TrkStrawHit* nearhit;
	if(ihigh != kres._hits.end())
	  nearhit = *ihigh;
	else
	  nearhit = *ilow;
	TrkT0 hitt0 = nearhit->hitT0();
	double mom = kres._krep->momentum(nearhit->fltLen()).mag();
	double beta = kres._tdef.particle().beta(mom);
	double tflt = (hflt-nearhit->fltLen())/(beta*CLHEP::c_light);
// update the time in the TrkT0 object
	hitt0._t0 += tflt;
// create the hit object.  Assume we're at the last iteration over added error
	TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,hitt0,hflt,_herr.back(),_maxdriftpull);
	assert(trkhit != 0);
// allow the hit to update its own ambiguity for now: eventually we should get the resolver to do this, FIXME!!!
	trkhit->setAmbigUpdate(true);
// must be initialy active for KalRep to process correctly
	trkhit->setActivity(true);
// flag the added hit
	trkhit->setUsability(3);
// add the hit to the track and the fit
	kres._krep->addHot(trkhit);
	kres._hits.push_back(trkhit);
// create intersections for the material of this hit and add those to the track
	DetIntersection wallinter;
	if(trkhit->wallElem().reIntersect(reftraj,wallinter))
	  kres._krep->addInter(wallinter);	
	DetIntersection gasinter;
	if(trkhit->gasElem().reIntersect(reftraj,gasinter))
	  kres._krep->addInter(gasinter);
// check the raw residual: This call works because the HOT isn't yet processed as part of the fit.
        double chi = fabs(trkhit->residual()/trkhit->hitRms());
//if it's outside limits, deactivate the HOT
	if(chi > maxchi || !trkhit->physicalDrift(maxchi))	
	  trkhit->setActivity(false);
// now that we've got the residual, we can turn of auto-ambiguity resolution
	trkhit->setAmbigUpdate(false);
      }
// refit the last iteration of the track
      fitIteration(kres,_herr.size()-1);
      kres._krep->addHistory(kres._fit,"AddHits");
    }
  }

  void KalFit::fitTrack(KalFitResult& kres){
    // loop over external hit errors, ambiguity assignment, t0 toleratnce
    for(size_t iherr=0;iherr < _herr.size(); ++iherr){
      fitIteration(kres,iherr);
      if(! kres._fit.success())break;
    }
    if(kres._krep != 0) kres._krep->addHistory(kres._fit,"KalFit");
  }

  void KalFit::fitIteration(KalFitResult& kres, size_t iherr){
    // update the external hit errors.  This isn't strictly necessary on the 1st iteration.
    for(std::vector<TrkStrawHit*>::iterator itsh = kres._hits.begin(); itsh != kres._hits.end(); ++itsh){
      (*itsh)->setExtErr(_herr[iherr]);
    }
    // update t0, and propagate it to the hits
    double oldt0 = kres._krep->t0()._t0;
    kres._nt0iter = 0;
    unsigned niter(0);
    bool changed(true);
					// 2015-04-12 P.Murat
    if (_useDoublets) {
      kres._decisionMode = _decisionMode;
    }
    kres._fit = TrkErrCode::succeed;
    while(kres._fit.success() && changed && niter < maxIterations()){
      changed = false;
      _ambigresolver[iherr]->resolveTrk(kres);
      kres._krep->resetFit();
      kres.fit();
      if(! kres._fit.success())break;
      if(_updatet0){
	updateT0(kres);
	changed |= fabs(kres._krep->t0()._t0-oldt0) > _t0tol[iherr];
	oldt0 = kres._krep->t0()._t0;
      }
      // drop outlyers
      if(_weedhits){
	kres._nweediter = 0;
	changed |= weedHits(kres);
      }
      niter++;
    }
    kres._ninter = kres._krep->intersections();
  }

  bool
  KalFit::fitable(TrkDef const& tdef){
    return tdef.strawHitIndices().size() >= _minnstraws;
  }
  
  void
  KalFit::makeHits(KalFitResult& kres,TrkT0 const& t0) {
    const Tracker& tracker = getTrackerOrThrow();
    TrkDef const& tdef = kres._tdef;
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
      TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,hitt0,fltlen,_herr.front(),_maxdriftpull);
      assert(trkhit != 0);
    // set the initial ambiguity based on the input
      trkhit->setAmbig(tdef.strawHitIndices()[iind]._ambig);
    // refine the flightlength, as otherwise hits in the same plane are at exactly the same flt, which can cause problems
      TrkErrCode pstat = trkhit->updatePoca(&tdef.helix());
      if(pstat.failure()){
        trkhit->setActivity(false);
      }
      kres._hits.push_back(trkhit);
    }
 // sort the hits by flightlength
    std::sort(kres._hits.begin(),kres._hits.end(),fltlencomp(tdef.fitdir().fitDirection()));
  }

  void
  KalFit::makeMaterials(KalFitResult& kres) {
    TrkDef const& tdef = kres._tdef;
    for(std::vector<TrkStrawHit*>::iterator ihit=kres._hits.begin();ihit!=kres._hits.end();ihit++){
      TrkStrawHit* trkhit = *ihit;
      // create wall and gas intersection objects from each straw hit (active or not)
      DetIntersection wallinter;
      wallinter.delem = 0;
      wallinter.pathlen = trkhit->fltLen();
      DetIntersection gasinter;
      gasinter.delem = 0;
      gasinter.pathlen = trkhit->fltLen();
      if(trkhit->wallElem().reIntersect(&tdef.helix(),wallinter))
	kres._detinter.push_back(wallinter);
      if(trkhit->gasElem().reIntersect(&tdef.helix(),gasinter))
	kres._detinter.push_back(gasinter);
    }
  }

  bool
  KalFit::weedHits(KalFitResult& kres) {
    // Loop over HoTs and find HoT with largest contribution to chi2.  If this value
    // is greater than some cut value, deactivate that HoT and reFit
    bool retval(false);
    double worst = -1.;
    TrkStrawHit* worstHot = 0;
    for (std::vector<TrkStrawHit*>::iterator iter = kres._hits.begin(); iter != kres._hits.end(); ++iter){
      TrkStrawHit* iHot = *iter;
      if (iHot->isActive()) {
        double resid, residErr;
        if(iHot->resid(resid, residErr, true)){
          double value = fabs(resid/residErr);
          if (value > _maxhitchi && value > worst) {
            worst = value;
            worstHot = iHot;
          }
        }
      }
    }
    if(0 != worstHot){
      retval = true;
      worstHot->setActivity(false);
      worstHot->setUsability(5); // positive usability allows hot to be re-enabled later
      kres.fit();
      kres._krep->addHistory(kres._fit, "HitWeed");
      // Recursively iterate
      kres._nweediter++;
      if (kres._fit.success() && kres._nweediter < _maxweed ) {
        retval |= weedHits(kres);
      }
    }
    return retval;
  }
  
  bool
  KalFit::unweedHits(KalFitResult& kres, double maxchi) {
    // Loop over inactive HoTs and find the one with the smallest contribution to chi2.  If this value
    // is less than some cut value, reactivate that HoT and reFit
    bool retval(false);
    double best = 1.e12;
    TrkStrawHit* bestHot = 0;
    for (std::vector<TrkStrawHit*>::iterator iter = kres._hits.begin(); iter != kres._hits.end(); ++iter){
      TrkStrawHit* iHot = *iter;
      if (!iHot->isActive()) {
        double resid, residErr;
        if(iHot->resid(resid, residErr, true)){
          double chival = fabs(resid/residErr);
  // test both for a good chisquared and for the drift radius to be physical
          if (chival < maxchi && iHot->physicalDrift(maxchi) && chival < best) {
            best = chival;
            bestHot = iHot;
          }
        }
      }
    }
    if(0 != bestHot){
      retval = true;
      bestHot->setActivity(true);
      bestHot->setUsability(4);

      if (_useDoublets == 1) {
      // 2015-04-12 P.Murat: mode doublets could be recovered, for now - a hack: 
      // assume that unweedHits is called only once from TrkpatRec in the very end
	kres._decisionMode = _decisionMode;
	int last = _herr.size()-1;
	_ambigresolver[last]->resolveTrk(kres);
      }

      kres.fit();
      kres._krep->addHistory(kres._fit, "HitUnWeed");
      // Recursively iterate
      kres._nunweediter++;
      if (kres._fit.success() && kres._nunweediter < _maxweed  ) {
        retval |= unweedHits(kres,maxchi);
      }
    }
    return retval;
  }

  BField const&
  KalFit::bField() const {
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
   
  const TrkVolume* 
  KalFit::trkVolume(trkDirection trkdir) const {
    //FIXME!!!!
    return 0;
  }

  void
  KalFit::initT0(TrkDef const& tdef, TrkT0& t0) {
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

  bool
  KalFit::updateT0(KalFitResult& kres){
    using namespace boost::accumulators;
    bool retval(false);
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    KalRep* krep = kres._krep;
// need to have a valid fit
    if(krep->fitValid()){
// find the global fltlen associated with z=0. 
      double flt0(0.0);
      bool converged = TrkHelixUtils::findZFltlen(krep->traj(),0.0,flt0);
      if(converged){
	std::vector<double> hitt0; // store t0, to allow outlyer removal
	std::vector<double> hitt0err;
	hitt0.reserve(kres._hits.size());
	hitt0err.reserve(kres._hits.size());
	// loop over the hits
	for(std::vector<TrkStrawHit*>::iterator ihit= kres._hits.begin();ihit != kres._hits.end(); ihit++){
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
	    updateHitTimes(kres);
	    retval = true;
	  }
	}
      }
    }
    return retval;
  }

  void
  KalFit::updateHitTimes(KalFitResult& kres) {
  // compute the time the track came closest to the wire for each hit, starting from t0 and working out.
  // this function allows for momentum change along the track.
  // find the bounding hits on either side of this
    std::sort(kres._hits.begin(),kres._hits.end(),fltlencomp(kres._tdef.fitdir().fitDirection()));
    std::vector<TrkStrawHit*>::iterator ihigh;
    std::vector<TrkStrawHit*>::reverse_iterator ilow;
    findBoundingHits(kres._hits,kres._krep->flt0(),ilow,ihigh);
    // reset all the hit times
    double hflt = kres._krep->flt0();
    TrkT0 hitt0 = kres._krep->t0();
    for(std::vector<TrkStrawHit*>::iterator ihit= ihigh;ihit != kres._hits.end(); ++ihit){
      TrkStrawHit* hit = *ihit;
// particle momentum at this point, using the full fit
      double mom = kres._krep->momentum(hit->fltLen()).mag();
// relativistic velocity from that
      double beta = kres._tdef.particle().beta(mom);
// particle transit time to this hit from the reference
      double tflt = (hit->fltLen()-hflt)/(beta*CLHEP::c_light);
// update the time in the TrkT0 object
      hitt0._t0 += tflt;
      (*ihit)->updateHitT0(hitt0);
// update the reference flightlength
      hflt = hit->fltLen();
    }
// now the same, moving backwards
    hflt = kres._krep->flt0();
    hitt0 = kres._krep->t0();
    for(std::vector<TrkStrawHit*>::reverse_iterator ihit= ilow;ihit != kres._hits.rend(); ++ihit){
      TrkStrawHit* hit = *ihit;
      double mom = kres._krep->momentum(hit->fltLen()).mag();
      double beta = kres._tdef.particle().beta(mom);
      double tflt = (hit->fltLen()-hflt)/(beta*CLHEP::c_light);
      hitt0._t0 += tflt;
      (*ihit)->updateHitT0(hitt0);
      hflt = hit->fltLen();
    }
  }

  void
  KalFit::findBoundingHits(std::vector<TrkStrawHit*>& hits,double flt0,
    std::vector<TrkStrawHit*>::reverse_iterator& ilow,
    std::vector<TrkStrawHit*>::iterator& ihigh) {
    ilow = hits.rbegin();
    ihigh = hits.begin();
    while(ilow != hits.rend() && (*ilow)->fltLen() > flt0 )++ilow;
    while(ihigh != hits.end() && (*ihigh)->fltLen() < flt0 )++ihigh;
  }

}

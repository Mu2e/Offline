//
// Class to perform BaBar Kalman fit
// Original author: Dave Brown LBNL 2012
//
// $Id: KalFit.cc,v 1.43 2014/08/22 16:10:41 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/22 16:10:41 $
//
#include "TrkReco/inc/KalFit.hh"
#include "TrkReco/inc/PanelAmbigResolver.hh"
#include "TrkReco/inc/PocaAmbigResolver.hh"
#include "TrkReco/inc/HitAmbigResolver.hh"
#include "TrkReco/inc/FixedAmbigResolver.hh"
#include "TrkReco/inc/DoubletAmbigResolver.hh"
#include "TrkReco/inc/TrkUtilities.hh"
#include "Mu2eBTrk/inc/BaBarMu2eField.hh"
#include "Mu2eBTrk/inc/Mu2eDetectorModel.hh"
//geometry
#include "BTrkHelper/inc/BTrkHelper.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
// tracker
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/KalmanTrack/KalBend.hh"
#include "BTrk/KalmanTrack/KalMaterial.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/BaBar/ErrLog.hh"
#include "BTrk/BField/BFieldFixed.hh"
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "BTrk/DetectorModel/DetMaterial.hh"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
//CLHEP
#include "CLHEP/Vector/ThreeVector.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <set>

using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::HepVector;
using CLHEP::HepSymMatrix;

namespace mu2e
{
// comparison functor for ordering hits.  This should operate on TrkHit, FIXME!
  struct fcomp : public binary_function<TrkHit*, TrkHit*, bool> {
    bool operator()(TrkHit* x, TrkHit* y) {
      return x->fltLen() < y->fltLen();
    }
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
  KalFit::KalFit(fhicl::ParameterSet const& pset) :
// KalFit parameters
    _debug(pset.get<int>("debugLevel",3)),
    _maxhitchi(pset.get<double>("maxhitchi",3.5)),
    _maxpull(pset.get<double>("maxPull",5)),
    // t0 parameters
    _initt0(pset.get<bool>("initT0",true)),
    _useTrkCaloHit(pset.get<bool>("useTrkCaloHit",false)),
    _updatet0(pset.get<bool>("updateT0",true)),
    _t0tol(pset.get< vector<double> >("t0Tolerance")),
    _t0errfac(pset.get<double>("t0ErrorFactor",1.2)),
    _mint0doca(pset.get<double>("minT0DOCA",-0.2)),
    _t0nsig(pset.get<double>("t0window",2.5)),
    _dtoffset(pset.get<double>("dtOffset")),
    _strHitW(pset.get<double>("strawHitT0Weight")),
    _calHitW(pset.get<double>("caloHitT0Weight")),
    //
    _minnstraws(pset.get<unsigned>("minnstraws",15)),
    _maxmatfltdiff(pset.get<double>("MaximumMaterialFlightDifference",1000.0)), // mm separation in flightlength
    _weedhits(pset.get<vector<bool> >("weedhits")),
    _herr(pset.get< vector<double> >("hiterr")),
    _ambigstrategy(pset.get< vector<int> >("ambiguityStrategy")),
    _addmaterial(pset.get<vector<bool> >("AddMaterial")),
    _resolveAfterWeeding(pset.get<bool>("ResolveAfterWeeding",false)),
    _exup((extent)pset.get<int>("UpstreamExtent",noextension)),
    _exdown((extent)pset.get<int>("DownstreamExtent",noextension)),
    _bfield(0)
  {
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
    _mingap = pset.get<double>("mingap",1.0);
    _minfltlen = pset.get<double>("MinFltLen",0.1);
    _minmom = pset.get<double>("MinMom",10.0);
    _fltepsilon = pset.get<double>("FltEpsilon",0.001);
    _divergeflt = pset.get<double>("DivergeFlt",1.0e3);
    _mindot = pset.get<double>("MinDot",0.0);
    _maxmomdiff = pset.get<double>("MaxMomDiff",0.5);
    _momfac = pset.get<double>("MomFactor",0.0);
    _maxpardif[0] = _maxpardif[1] = pset.get<double>("MaxParameterDifference",1.0);

    _mindof = pset.get<double>("MinNDOF",10);
    // this config belongs in the BField integrator, FIXME!!!
    _bintconfig._maxRange = pset.get<double>("BFieldIntMaxRange",1.0e5); // 100 m
    _bintconfig._intTolerance = pset.get<double>("BFieldIntTol",0.01); // 10 KeV
    _bintconfig._intPathMin = pset.get<double>("BFieldIntMin",20.0); // 20 mm
    _bintconfig._divTolerance = pset.get<double>("BFieldIntDivTol",0.05); // 50 KeV
    _bintconfig._divPathMin = pset.get<double>("BFieldIntDivMin",50.0); // 50 mm
    _bintconfig._divStepCeiling = pset.get<double>("BFieldIntDivMax",500.0); // 500 mm
    // field integral errors.  This is commented out as it hasn't been shown to improve the fit
    //    double perr = pset.get<double>("BendCorrErrFrac",0.0); // fractional accuracy of trajectory
    //    double berr = pset.get<double>("BFieldMapErr",0.0); // mapping and interpolation error
//    KalBend::setErrors(perr,berr);
    // make sure we have at least one entry for additional errors
    if(_herr.size() <= 0) throw cet::exception("RECO")<<"mu2e::KalFit: no hit errors specified" << endl;
    if(_herr.size() != _ambigstrategy.size()) throw cet::exception("RECO")<<"mu2e::KalFit: inconsistent ambiguity resolution hiterr" << endl;
    if(_herr.size() != _t0tol.size()) throw cet::exception("RECO")<<"mu2e::KalFit: inconsistent ambiguity resolution t0" << endl;
    if(_herr.size() != _weedhits.size()) throw cet::exception("RECO")<<"mu2e::KalFit: inconsistent ambiguity resolution WeedHits" << endl;
    if(_herr.size() != _addmaterial.size()) throw cet::exception("RECO")<<"mu2e::KalFit: inconsistent ambiguity resolution AddMaterial" << endl;
    // Search for explicit resolver parameter sets.  These may not be used
    fhicl::ParameterSet const& fixedPset = pset.get<fhicl::ParameterSet>("FixedAmbigResolver",fhicl::ParameterSet());
    fhicl::ParameterSet const& hitPset = pset.get<fhicl::ParameterSet>("HitAmbigResolver",fhicl::ParameterSet());
    fhicl::ParameterSet const& panelPset = pset.get<fhicl::ParameterSet>("PanelAmbigResolver",fhicl::ParameterSet());
    fhicl::ParameterSet const& pocaPset = pset.get<fhicl::ParameterSet>("POCAAmbigResolver",fhicl::ParameterSet());
    fhicl::ParameterSet const& doubletPset = pset.get<fhicl::ParameterSet>("DoubletAmbigResolver",fhicl::ParameterSet());
// construct the explicit ambiguity resolvers, 1 instance per iteration
    size_t niter = _ambigstrategy.size();
    for(size_t iter=0; iter<niter; ++iter) {
      int Final = iter==niter-1 ? 1 : 0;
      AmbigResolver* ar(0);
      switch (_ambigstrategy[iter]) {
      case fixedambig:
        ar = new FixedAmbigResolver(fixedPset,_herr[iter]);
        break;
      case hitambig:
        ar = new HitAmbigResolver(hitPset,_herr[iter]);
        break;
      case panelambig:
        ar = new PanelAmbig::PanelAmbigResolver(panelPset,_herr[iter],iter);
        break;
      case pocaambig:
        ar = new PocaAmbigResolver(pocaPset,_herr[iter]);
        break;
      case doubletambig: // 4
        ar = new DoubletAmbigResolver(doubletPset,_herr[iter],iter,Final);
        break;
      default:
        break;
      }
      if(ar != 0)
        _ambigresolver.push_back(ar);
      else
        throw cet::exception("RECO")<<"mu2e::KalFit: unknown ambiguity resolver " << _ambigstrategy[iter] << " for iteration " << iter << endl;
    }
  }

  KalFit::~KalFit(){
    for(size_t iambig=0;iambig<_ambigresolver.size();++iambig){
      delete _ambigresolver[iambig];
    }
    delete _bfield;
  }

//-----------------------------------------------------------------------------
// create the track (KalRep) from a track seed
//-----------------------------------------------------------------------------
  void KalFit::makeTrack(const StrawHitCollection* shcol, KalSeed const& kseed, KalRep*& krep) {
// test if fitable
    if(fitable(kseed)){
      // find the segment at the 0 flight
      double flt0 = kseed.flt0();
      auto kseg = kseed.segments().end();
      for(auto iseg= kseed.segments().begin(); iseg != kseed.segments().end(); ++iseg){
	if(iseg->fmin() <= flt0 && iseg->fmax() > flt0){
	  kseg = iseg;
	  break;
	}
      }
      if(kseg == kseed.segments().end()){
	std::cout << "Helix segment range doesn't cover flt0" << std::endl;
	kseg = kseed.segments().begin();
      }
      // create a trajectory from the seed. This shoudl be a general utility function that
      // can work with multi-segment seeds FIXME!
      // create CLHEP objects from seed native members.  This will
      // go away when we switch to SMatrix FIXME!!!
      HepVector pvec(5,0);
      HepSymMatrix pcov(5,0);
      kseg->helix().hepVector(pvec);
      kseg->covar().symMatrix(pcov);
      // Create the traj from these
      HelixTraj htraj(pvec,pcov);
      // create the hits
      TrkStrawHitVector tshv;
      makeTrkStrawHits(shcol, htraj, kseed.hits(), tshv);
      
   // Find the wall and gas material description objects for these hits
      std::vector<DetIntersection> detinter;
      if(_matcorr)makeMaterials(tshv,htraj,detinter);
   // Create the BaBar hit list, and fill it with these hits.  The BaBar list takes ownership
      // We should use the TrkHit vector everywhere, FIXME!
      std::vector<TrkHit*> thv(0);
      for(auto ihit = tshv.begin(); ihit != tshv.end(); ++ihit){
        thv.push_back(*ihit);
        if (_debug>2) { (*ihit)->print(std::cout); }
      }
      if (_useTrkCaloHit){    //use the TrkCaloHit to initialize the t0?
	//create the TrkCaloHit
	TrkCaloHit* tch(0);
	makeTrkCaloHit(kseed, tch);
	if (tch != 0) thv.push_back(tch);
      }
 


      TrkT0 t0(kseed.t0());
      // create Kalman rep
      krep = new KalRep(htraj, thv, detinter, *this, kseed.particle(), t0, flt0);
      assert(krep != 0);
      if(_initt0){
	initT0(krep);
      }
// initialize history list
      krep->addHistory(TrkErrCode(),"KalFit creation");
// now fit
      TrkErrCode fitstat = fitTrack(krep);
      krep->addHistory(fitstat,"KalFit fit");
// extend the fit
      if(fitstat.success()){
	fitstat = extendFit(krep);
	krep->addHistory(fitstat,"KalFit extension");
      }
    }
  }

  void KalFit::addHits(KalRep* krep,const StrawHitCollection* shcol, std::vector<StrawHitIndex> indices, double maxchi) {
    //2017-05-02: Gianipez. In this function inten
  // fetcth the DetectorModel
   Mu2eDetectorModel const& detmodel{ art::ServiceHandle<BTrkHelper>()->detectorModel() };
// there must be a valid Kalman fit to add hits to
    if(krep != 0 && indices.size() > 0 && krep->fitStatus().success()){
      //      TrkHitVector thv;
      //      thv = krep->hitVector();
      ConditionsHandle<TrackerCalibrations> tcal("ignored");
      const Tracker& tracker = getTrackerOrThrow();
      TrkHitVector::iterator ihigh;
      TrkHitVector::reverse_iterator ilow;
// use the reference trajectory, as that's what all the existing hits do
      const TrkDifPieceTraj* reftraj = krep->referenceTraj();
      for(unsigned iind=0;iind<indices.size(); ++iind){
        size_t istraw = indices[iind];
        const StrawHit& strawhit(shcol->at(istraw));
        const Straw& straw = tracker.getStraw(strawhit.strawIndex());
// estimate  initial flightlength
        double hflt(0.0);
        TrkHelixUtils::findZFltlen(*reftraj,straw.getMidPoint().z(),hflt);
// find the bounding sites near this hit, and extrapolate to get the hit t0
        findBoundingHits(krep,hflt,ilow,ihigh);
        const TrkHit* nearhit;
        if(ihigh != krep->hitVector().end())
          nearhit = *ihigh;
        else
          nearhit = *ilow;
        TrkT0 hitt0 = nearhit->hitT0();
        double mom = krep->momentum(nearhit->fltLen()).mag();
        double beta = krep->particleType().beta(mom);
        double tflt = (hflt-nearhit->fltLen())/(beta*CLHEP::c_light);
// update the time in the TrkT0 object
        hitt0._t0 += tflt;
// create the hit object.  Assume we're at the last iteration over added error
        TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,hitt0,hflt,_herr.back(),
					      _maxpull,_strHitW, _mint0doca);
        assert(trkhit != 0);
// allow the hit to update its own ambiguity for now: eventually we should get the resolver to do this, FIXME!!!
        trkhit->setAmbigUpdate(true);
        trkhit->setFlag(TrkHit::addedHit);
// must be initialy active for KalRep to process correctly
        trkhit->setActivity(true);
// add the hit to the track
        krep->addHit(trkhit);
// check the raw residual: This call works because the HOT isn't yet processed as part of the fit.
        double chi = fabs(trkhit->residual()/trkhit->hitRms());
//if it's outside limits, deactivate the HOT
        if(chi > maxchi || (!trkhit->isPhysical(maxchi)))
          trkhit->setActivity(false);
// now that we've got the residual, we can turn of auto-ambiguity resolution
        trkhit->setAmbigUpdate(false);
   // find the DetElem associated this straw
        const DetStrawElem* strawelem = detmodel.strawElem(trkhit->straw());
// see if this KalRep already has a KalMaterial with this element: if not, add it
        bool hasmat(false);
        std::vector<const KalMaterial*> kmats;
        krep->findMaterialSites(strawelem,kmats);
// if this is a reflecting track the same material can appear multiple times: check the flight lengths
        if(kmats.size() > 0){
          for(auto kmat: kmats) {
            if( fabs( kmat->globalLength() - trkhit->fltLen()) < _maxmatfltdiff){
              hasmat = true;
              break;
            }
          }
        }
        if(!hasmat){
          // create intersection object for this element; it includes all materials
          DetIntersection strawinter(strawelem, krep->referenceTraj(),trkhit->fltLen());
          strawinter.thit = trkhit;
          // compute initial intersection: this gets updated each fit iteration
          strawelem->reIntersect(krep->referenceTraj(),strawinter);
          krep->addInter(strawinter);
        }
      }
// refit the last iteration of the track
      TrkErrCode fitstat = fitIteration(krep,_herr.size()-1);
      krep->addHistory(fitstat,"AddHits");
    }
  }
//
  TrkErrCode KalFit::fitTrack(KalRep* krep) {
    // loop over external hit errors, ambiguity assignment, t0 toleratnce
    TrkErrCode fitstat;
    for(size_t iherr=0;iherr < _herr.size(); ++iherr) {
      fitstat = fitIteration(krep,iherr);
      if(_debug > 0) cout << "Iteration " << iherr 
      << " NDOF = " << krep->nDof() 
      << " Fit Status = " <<  fitstat << endl;
      if(!fitstat.success())break;
    }
    return fitstat;
  }

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
  TrkErrCode KalFit::fitIteration(KalRep* krep, size_t iter) {
    // update the external hit errors.  This isn't strictly necessary on the 1st iteration.
    TrkHitVector* thv   = &(krep->hitVector());
    for (auto itsh=thv->begin();itsh!=thv->end(); ++itsh){
      (*itsh)->setTemperature(_herr[iter]);
    }
    // update t0, and propagate it to the hits
    double oldt0 = krep->t0()._t0;
    unsigned niter(0);
    bool changed(true);
    TrkErrCode retval = TrkErrCode::succeed;
    while(retval.success() && changed && ++niter < maxIterations()){
//-----------------------------------------------------------------------------
// convention: resolve drift signs before the fit with respect to the trajectory
// determined at the previous iteration
//-----------------------------------------------------------------------------
      changed = _ambigresolver[iter]->resolveTrk(krep);
      // force a refit
      krep->resetFit();
      retval = krep->fit();
      if(! retval.success())break;
      // updates
      if(_updatet0){
	updateT0(krep);
        changed |= fabs(krep->t0()._t0-oldt0) > _t0tol[iter];
      }
      // drop outliers
      if(_weedhits[iter]){
        changed |= weedHits(krep,iter);
	changed |= unweedBestHit(krep,_maxhitchi);
      }
      // find missing materials
      unsigned nmat(0);
      if(_addmaterial[iter]){
	nmat = addMaterial(krep);
        changed |= nmat>0;
      }
      if(_debug > 1) std::cout << "Inner iteration " << niter << " changed = "
	<< changed << " t0 old " << oldt0 << " new " << krep->t0()._t0 
	<< " nmat = " << nmat << endl;
      oldt0 = krep->t0()._t0;
    }
    if(_debug > 1)
      std::cout << "Outer iteration " << iter << " stopped after "
      << niter << " iterations" << std::endl;
// make sure the fit is current
    if(!krep->fitCurrent())
      retval = krep->fit();
    return retval;
  }

  bool
  KalFit::fitable(TrkDef const& tdef){
    return tdef.strawHitIndices().size() >= _minnstraws;
  }

  bool
  KalFit::fitable(KalSeed const& kseed){
    return kseed.segments().size() > 0 && kseed.hits().size() >= _minnstraws;
  }

  void
  KalFit::makeTrkStrawHits(const StrawHitCollection* shcol, HelixTraj const& htraj,
			   std::vector<TrkStrawHitSeed>const& hseeds, TrkStrawHitVector& tshv ) {
    const Tracker& tracker = getTrackerOrThrow();
    // compute particle velocity to 
    for(auto ths : hseeds ){
      // create a TrkStrawHit from this seed.
      size_t index = ths.index();
      const StrawHit& strawhit(shcol->at(index));
      const Straw& straw = tracker.getStraw(strawhit.strawIndex());
      TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,ths.index(),ths.t0(),ths.trkLen(),
					    _herr.front(),_maxpull,_strHitW,_mint0doca);
      assert(trkhit != 0);
      // set the initial ambiguity
      trkhit->setAmbig(ths.ambig());
      // refine the flightlength, as otherwise hits in the same plane are at exactly the same flt, which can cause problems
      TrkErrCode pstat = trkhit->updatePoca(&htraj);
      if(pstat.failure()){
        trkhit->setActivity(false);
      }
      tshv.push_back(trkhit);
    }
 // sort the hits by flightlength
    std::sort(tshv.begin(),tshv.end(),fcomp());
  }

  void 
  KalFit::makeTrkCaloHit  (KalSeed const& kseed, TrkCaloHit *tch){
    if (kseed.caloCluster().get() != 0){
      HitT0 ht0;
      ht0._t0    = kseed.caloCluster()->time();
      ht0._t0err = 0.5;//dummy error FIXME!
      
      double fltlen(0);//dummy value FIXME! maybe I can use the dip angle from the kseed?
      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      Hep3Vector          cog = ch->geomUtil().mu2eToTracker(ch->geomUtil().diskToMu2e(
		   kseed.caloCluster()->diskId(), kseed.caloCluster()->cog3Vector())); 
      
      Hep3Vector const& clusterAxis = Hep3Vector(0, 0, 1);//FIX ME!
      double      crystalHalfLength = ch->caloInfo().crystalHalfLength();
      tch = new TrkCaloHit(*kseed.caloCluster().get(), cog, crystalHalfLength, clusterAxis, ht0, fltlen, _calHitW, _dtoffset);
    }
  }


  void
  KalFit::makeMaterials(TrkStrawHitVector const& tshv, HelixTraj const& htraj,std::vector<DetIntersection>& detinter) {
  // fetcth the DetectorModel
    Mu2eDetectorModel const& detmodel{ art::ServiceHandle<BTrkHelper>()->detectorModel() };
    // loop over strawhits and extract the straws
    for (auto trkhit : tshv) {
   // find the DetElem associated this straw
      const DetStrawElem* strawelem = detmodel.strawElem(trkhit->straw());
      // create intersection object for this element; it includes all materials
      DetIntersection strawinter;
      strawinter.delem = strawelem;
      strawinter.pathlen = trkhit->fltLen();
      strawinter.thit = trkhit;
      // compute initial intersection: this gets updated each fit iteration
      strawelem->reIntersect(&htraj,strawinter);
      detinter.push_back(strawinter);
    }
  }

  unsigned KalFit::addMaterial(KalRep* krep) {
    _debug>2 && std::cout << __func__ << " called " << std::endl;
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
    for(auto const& plane : ttracker.getPlanes()){
       _debug>2 && std::cout << __func__ << " plane " << plane.id() << std::endl;
      // if (!(plane.exists())) continue;
      // # of straws in a panel
      int nstraws = plane.getPanel(0).nStraws();
       _debug>2 && std::cout << __func__ << " nstraws " << nstraws << std::endl;
// get an approximate z position for this plane from the average position of the 1st and last straws
      // Hep3Vector s0 = plane.getPanel(0).getLayer(0).getStraw(0).getMidPoint();
      // plane id is id of 0th straw
      Hep3Vector s0 = plane.getPanel(0).getStraw(StrawId2(plane.id())).getMidPoint();
       _debug>2 && std::cout << __func__ << " s0 via panel " << s0 << std::endl;
      // funky convention for straw numbering in a layer FIXME!!!!
      // Hep3Vector sn = plane.getPanel(0).getLayer(1).getStraw(2*plane.getPanel(0).getLayer(1).nStraws()-1).getMidPoint();
      Hep3Vector sn = plane.getPanel(0).getStraw(nstraws-1).getMidPoint();
      _debug>2 && std::cout << __func__ << " sn via panel " << sn << std::endl;
      double pz = 0.5*(s0.z() + sn.z());
      _debug>2 && std::cout << __func__ << " an approximate z position for this plane " << plane.id() << " " << pz << std::endl;
// find the transverse position at this z using the reference trajectory
      double flt = krep->referenceTraj()->zFlight(pz);
      HepPoint pos = krep->referenceTraj()->position(flt);
      Hep3Vector posv(pos.x(),pos.y(),pos.z());
// see if this position is in the active region.  Double the straw radius to be generous
      double rho = posv.perp();
      double rmin = s0.perp()-2*strawradius;
      double rmax = sn.perp()+2*strawradius;
      if(rho > rmin && rho < rmax){
  // loop over panels
        for(auto const& panel : plane.getPanels()){
          if (_debug>2) {
            std::cout << __func__ << " panel " << panel.id() << std::endl;

            std::cout << __func__ << " printing all straws in layer 0 " << std::endl;
            std::cout << __func__ << " ";
            for (const auto straw_p : panel.getStrawPointers() ) {
              Straw const&       straw(*straw_p);
              StrawId sid = straw.id();
              if ( sid.getLayer() != 0 ) continue;
              std::cout.width(10);
              std::cout  << sid << ", ";
            }
            std::cout << std::endl;

            std::cout << __func__ << " printing all straws in layer 1 " << std::endl;
            std::cout << __func__ << " ";
            for (const auto straw_p : panel.getStrawPointers() ) {
              Straw const&       straw(*straw_p);
              StrawId sid = straw.id();
              if ( sid.getLayer() != 1 ) continue;
              std::cout.width(10);
              std::cout  << sid << ", ";
            }
            std::cout << std::endl;
          }
      // get the straw direction for this panel
          // Hep3Vector sdir = panel.getLayer(0).getStraw(0).getDirection();
          Hep3Vector sdir = panel.getStraw(0).getDirection();
      // get the transverse direction to this and z
          static Hep3Vector zdir(0,0,1.0);
          Hep3Vector pdir = sdir.cross(zdir);
     //  project the position along this
          double prho = posv.dot(pdir);
      // test for acceptance of this panel
          if(prho > rmin && prho < rmax) {
          // translate the transverse position into a rough straw number
          // nstraws is the number of straws in the panel
            int istraw = (int)rint(nstraws*(prho-s0.perp())/(sn.perp()-s0.perp()));
            // take a few straws around this
            for(int is = max(0,istraw-2); is<min(nstraws-1,istraw+2); ++is){
              _debug>2 && std::cout << __func__ << " taking a few straws, istraw, is " << istraw << ", " << is << std::endl;
              // must do this twice due to intrusion of layer on hierarchy FIXME!!!
              // std::cout << __func__ << " straw id l0 " << panel.getLayer(0).getStraw(is).id() << std::endl;
              // std::cout << __func__ << " straw id l1 " << panel.getLayer(1).getStraw(is).id() << std::endl;
               _debug>2 && std::cout << __func__ << " straw id l0 " << panel.getStraw(is).id() << std::endl;
              if ( panel.getStraw(is).id().getLayer()==0) {
                 _debug>2 && std::cout << __func__ << " straw id l0 by id "
                          << panel.getStraw(StrawId2(panel.id().asUint16()+is)).id() << std::endl;
              } else {
                 _debug>2 && std::cout << __func__ << " straw id l1 by id "
                          << panel.getStraw(StrawId2(panel.id().asUint16()+is)).id() << std::endl;
              }
              matstraws.insert(StrawFlight(panel.getStraw(is).index(),flt));
              ++nadded;
            }
          }
        }
      }
    }
// Now test if the Kalman rep hits these straws
    if(_debug>2)std::cout << "Found " << matstraws.size() << " unique possible straws " << " out of " << nadded << std::endl;
    for(auto const& strawflt : matstraws){
      const DetStrawElem* strawelem = detmodel.strawElem(strawflt._index);
      DetIntersection strawinter;
      strawinter.delem = strawelem;
      strawinter.pathlen = strawflt._flt;
      if(strawelem->reIntersect(krep->referenceTraj(),strawinter)){
// If the rep already has a material site for this element, skip it
        std::vector<const KalMaterial*> kmats;
        krep->findMaterialSites(strawelem,kmats);
        if(_debug>2)std::cout << "found intersection with straw " << strawelem->straw()->index() << " with "
        << kmats.size() << " materials " << std::endl;
// test material isn't on the track
        bool hasmat(false);
        for(auto kmat : kmats ){
          const DetStrawElem* kelem = dynamic_cast<const DetStrawElem*>(kmat->detIntersection().delem);
          if(kelem != 0){
            StrawFlight ksflt(kelem->straw()->index(),kmat->globalLength());
            if(_debug>2)std::cout << " comparing flights " << kmat->globalLength() << " and " << strawflt._flt << std::endl;
            if(!strawcomp.operator()(strawflt,ksflt)){
              if(_debug>2)std::cout << "operator returned false!!" << std::endl;
              // this straw is already on the track: stop
              hasmat = true;
              break;
            }
          }
        }
        if(kmats.size() == 0 || !hasmat) {
          if(_debug>2)std::cout << "Adding material element" << std::endl;
          // this straw doesn't have an entry in the Kalman fit: add it`
          DetIntersection detinter(strawelem, krep->referenceTraj(),strawflt._flt);
          krep->addInter(detinter);
          ++retval;
        }
      }
    }
    if(_debug>1)std::cout << "Added " << retval << " new material sites" << std::endl;
    return retval;
  }

  bool
  KalFit::weedHits(KalRep* krep, size_t iter) {
    // Loop over HoTs and find HoT with largest contribution to chi2.  If this value
    // is greater than some cut value, deactivate that HoT and reFit
    bool retval(false);
    double worst = -1.;
    //    TrkHit* worsthit = 0;
    TrkStrawHit  *worsthit = 0;
    TrkHitVector *thv      = &(krep->hitVector());

    for (auto ihit=thv->begin();ihit!=thv->end(); ++ihit){
      //      TrkHit* hit = *ihit;
      TrkStrawHit*hit = dynamic_cast<TrkStrawHit*>(*ihit);
      if (hit == 0)     continue;
      if (hit->isActive()) {
        double resid, residErr;
        if(hit->resid(resid, residErr, true)){
          double value = fabs(resid/residErr);
          if (value > _maxhitchi && value > worst) {
            worst = value;
            worsthit = hit;
          }
        }
      }
    }
    if(0 != worsthit){
      retval = true;
      worsthit->setActivity(false);
      worsthit->setFlag(TrkHit::weededHit);
      if (_resolveAfterWeeding) {
//-----------------------------------------------------------------------------
// _resolveAfterWeeding=0 makes changes in the logic fully reversible
//-----------------------------------------------------------------------------
        _ambigresolver[iter]->resolveTrk(krep);
      }
      TrkErrCode fitstat = krep->fit();
      krep->addHistory(fitstat, "HitWeed");
      // Recursively iterate
      if (fitstat.success() ) {
        retval |= weedHits(krep,iter);
      }
    }
    return retval;
  }

  bool
  KalFit::unweedHits(KalRep* krep, double maxchi) {

    bool retval = unweedBestHit(krep, maxchi);
    // if any hits were added, re-analyze ambiguity
    if (retval && _resolveAfterWeeding) {
      // 2015-04-12 P.Murat: '_resolveAfterWeeding' is here to make my changes fully reversible
      // I think, resolving ambiguities before each fit, makes a lot of sense
      //
      // Moved to after iteration: PanelAmbig resolver can change the state of hit resulting in infinte
      // loop if the resolver is called each iteration
      int last = _herr.size()-1;
      _ambigresolver[last]->resolveTrk(krep);
      if(!krep->fitCurrent()){
    // if this changed the track state, refit it
        krep->resetFit();
        TrkErrCode fitstat = krep->fit();
        krep->addHistory(fitstat, "HitUnWeedResolver");
      }
    }
    return retval;
  }

  bool
  KalFit::unweedBestHit(KalRep* krep, double maxchi) {
    // Loop over inactive HoTs and find the one with the smallest contribution to chi2.  If this value
    // is less than some cut value, reactivate that HoT and reFit
    bool      retval(false);
    double    best = 1.e12;
// no need to cast
    TrkHit* besthit = 0;
    const TrkHitVector* thv = &(krep->hitVector());
    for (auto ihit=thv->begin();ihit!=thv->end(); ++ihit){
      TrkHit* hit = *ihit;
      if (!hit->isActive()) {
        double resid, residErr;
        if(hit->resid(resid, residErr, true)){
          double chival = fabs(resid/residErr);
  // test both for a good chisquared and for the drift radius to be physical
          if (chival < maxchi && hit->isPhysical(maxchi) && chival < best) {
            best = chival;
            besthit = hit;
          }
        }
      }
    }
    if(0 != besthit){
      retval = true;
      besthit->setActivity(true);
      besthit->setFlag(TrkHit::unweededHit);
      TrkErrCode fitstat = krep->fit();
      if (fitstat.success() && besthit->isActive() ) {
	krep->addHistory(fitstat, "HitUnWeed");
	// Recursively iterate
        retval |= unweedBestHit(krep, maxchi);
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

// const StrawHitCollection* shcol,TrkParticle const& part,
// 		 TrkT0& t0,std::vector<StrawHitIndex> const& hits,
// 		 HelixTraj const& htraj   ) {

  void
  KalFit::initT0(KalRep*krep) {
    TrkT0 t0;
    using namespace boost::accumulators;
    // make an array of all the hit times, correcting for propagation delay
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    unsigned nind = krep->hitVector().size();
    std::vector<double> times;
    std::vector<double> timesweight;
    times.reserve(nind);
    timesweight.reserve(nind);
    // get flight distance of z=0
    double t0flt = krep->referenceTraj()->zFlight(0);//htraj.zFlight(0.0);
    // estimate the momentum at that point using the helix parameters.  This is
    // assumed constant for this crude estimate
    double loclen;
    double fltlen(0.0);
    const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(krep->referenceTraj()->localTrajectory(fltlen,loclen));
    double mom = TrkMomCalculator::vecMom(*htraj,bField(),t0flt).mag();
    // compute the particle velocity
    double vflt = krep->particleType().beta(mom)*CLHEP::c_light;
    // use the reference trajectory, as that's what all the existing hits do
    const TrkDifPieceTraj* reftraj = (krep->referenceTraj());
    // loop over hits
    double      htime(0);    
    for(auto ith=krep->hitVector().begin(); ith!=krep->hitVector().end(); ++ith){
      (*ith)->trackT0Time(htime, t0flt, reftraj, vflt);
      times.push_back(htime);
      timesweight.push_back((*ith)->t0Weight());
    }

    // find the median time
    accumulator_set<double,stats<tag::weighted_variance >,double >         wmean;
    //fill the accumulator using the weights
    int nhits(times.size());
    for (int i=0; i<nhits; ++i){
      wmean(times.at(i), weight=timesweight.at(i));
    }
    t0._t0    = extract_result<tag::weighted_mean>(wmean);
    t0._t0err = sqrt(extract_result<tag::weighted_variance>(wmean)/nhits);
    //set the new T0
    krep->setT0(t0, t0flt);
  }

  bool
  KalFit::updateT0(KalRep* krep){
    using namespace boost::accumulators;
    TrkHitVector *thv = &(krep->hitVector());
    bool retval(false);
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
// need to have a valid fit
    if(krep->fitValid()){
// find the global fltlen associated with z=0.
      double flt0(0.0);
      bool converged = TrkHelixUtils::findZFltlen(krep->traj(),0.0,flt0);
      if(converged){
        std::vector<double> hitt0; // store t0, to allow outlyer removal
        std::vector<double> hitt0err;
        size_t nhits = krep->hitVector().size();
        hitt0.reserve(nhits);
        hitt0err.reserve(nhits);
        // loop over the hits
        for(auto ihit=thv->begin(); ihit != thv->end(); ihit++){
          TrkHit* hit = *ihit;
          if(hit->isActive() && hit->hasResidual()){
            // find the residual, exluding this hits measurement
            double resid,residerr;
	    double pTime, doca;//propagation-time
	    CLHEP::Hep3Vector trjDir(krep->traj().direction(hit->fltLen()));
            if(krep->resid(hit,resid,residerr,true)){
	      if (hit->signalPropagationTime(pTime, doca, resid, residerr, trjDir)){	      
                // subtracting hitT0 makes this WRT the previous track t0
                hitt0.push_back(hit->time() - pTime - hit->hitT0()._t0);
                // assume residual error dominates
                hitt0err.push_back(residerr);
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
            updateHitTimes(krep);
            retval = true;
          }
        }
      }
    }
    return retval;
  }
  
  void
  KalFit::updateHitTimes(KalRep* krep) {
  // compute the time the track came closest to the wire for each hit, starting from t0 and working out.
  // this function allows for momentum change along the track.
  // find the bounding hits on either side of this
    TrkHitVector *thv = &(krep->hitVector());
    std::sort(thv->begin(),thv->end(),fcomp());
    TrkHitVector::iterator ihigh;
    TrkHitVector::reverse_iterator ilow;
    findBoundingHits(krep, krep->flt0(),ilow,ihigh);
    // reset all the hit times
    double flt0 = krep->flt0();
    TrkT0 hitt0 = krep->t0();
    for(TrkHitVector::iterator ihit= ihigh;ihit != thv->end(); ++ihit){
      TrkHit* hit = *ihit;
      double flt1 = hit->fltLen();
// particle transit time to this hit from the reference
      double tflt = krep->transitTime(flt0, flt1);
// update the time in the TrkT0 object
      hitt0._t0 += tflt;
      //      (*ihit)->updateHitT0(hitt0);
      (*ihit)->setHitT0(hitt0);
// update the reference flightlength
      flt0 = flt1;
    }
// now the same, moving backwards.
    flt0 = krep->flt0();
    hitt0 = krep->t0();
    for(TrkHitVector::reverse_iterator ihit= ilow;ihit != thv->rend(); ++ihit){
      TrkHit* hit = *ihit;
      double flt1 = hit->fltLen();
      double tflt = krep->transitTime(flt0, flt1);
      hitt0._t0 += tflt;
      //      (*ihit)->updateHitT0(hitt0);
      (*ihit)->setHitT0(hitt0);
      flt0 = flt1;
    }
  }

  void
  KalFit::findBoundingHits(KalRep* krep,double flt0,
			   TrkHitVector::reverse_iterator& ilow,
			   TrkHitVector::iterator& ihigh) {
    TrkHitVector* hits = &(krep->hitVector());
    ilow = hits->rbegin();
    ihigh = hits->begin();
    while(ilow != hits->rend() && (*ilow)->fltLen() > flt0 )++ilow;
    while(ihigh != hits->end() && (*ihigh)->fltLen() < flt0 )++ihigh;
  }

// attempt to extend the fit to the specified location
  TrkErrCode KalFit::extendFit(KalRep* krep) {
    TrkErrCode retval;
    // find the downstream and upstream Z positions to extend to
    if(_exdown != noextension){
      double downz = extendZ(_exdown);
    // convert to flightlength using the fit trajectory
      double downflt = krep->pieceTraj().zFlight(downz);
    // actually extend the track
      retval = krep->extendThrough(downflt);
    }
    // same for upstream extension
    if(retval.success() && _exup != noextension){
      double upz = extendZ(_exup);
      double upflt = krep->pieceTraj().zFlight(upz);
      retval = krep->extendThrough(upflt);
    }
    return retval;
  }

  double KalFit::extendZ(extent ex) {
    double retval(0.0);
    if(ex == target){
      GeomHandle<StoppingTarget> target;
      GeomHandle<DetectorSystem> det;
      retval = det->toDetector(target->centerInMu2e()).z() - 0.5*target->cylinderLength();
    } else if(ex == ipa) {
    // the following is wrong FIXME!!
      GeomHandle<StoppingTarget> target;
      GeomHandle<DetectorSystem> det;
      retval = det->toDetector(target->centerInMu2e()).z() +  0.5*target->cylinderLength();
    } else if(ex == tracker) {
      retval = 0.0;
    } else if(ex == calo) {
      GeomHandle<Calorimeter> cg;
      return cg->caloInfo().envelopeZ1();
    }
    return retval;
  }

  void KalFit::findTrkCaloHit(KalRep*krep, TrkCaloHit*tch){
    for(auto ith=krep->hitVector().begin(); ith!=krep->hitVector().end(); ++ith){
      TrkCaloHit* tsh = dynamic_cast<TrkCaloHit*>(*ith);
      if(tsh != 0) {
	tch = tsh;
	break;
      }
    }

  }
  

}

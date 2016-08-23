//
// Class to perform BaBar Kalman fit
// Original author: Dave Brown LBNL 2012
//
// $Id: KalFit.cc,v 1.43 2014/08/22 16:10:41 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/22 16:10:41 $
//

// the following has to come before other BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/KalFit.hh"
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
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
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
  KalFit::KalFit(fhicl::ParameterSet const& pset, TrkFitDirection const& fitdir) :
// KalFit parameters
    _debug(pset.get<int>("debugLevel",0)),
    _maxhitchi(pset.get<double>("maxhitchi",3.5)),
    _maxdriftpull(pset.get<double>("maxDriftPull",10)),
    // t0 parameters
    _initt0(pset.get<bool>("initT0",true)),
    _updatet0(pset.get<bool>("updateT0",true)),
    _t0tol(pset.get< vector<double> >("t0Tolerance")),
    _t0errfac(pset.get<double>("t0ErrorFactor",1.2)),
    _mint0doca(pset.get<double>("minT0DOCA",-0.2)),
    _t0nsig(pset.get<double>("t0window",2.5)),
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
    _fdir(fitdir),
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
// create the track (KalRep) from the track definition. 
//-----------------------------------------------------------------------------
  void KalFit::makeTrack(const StrawHitCollection* shcol, TrkDef& tdef, KalRep*& krep) {
// test if fitable
    if(fitable(tdef)){
// if requested, initialize t0
      if(_initt0)
        initT0(shcol,tdef);
// create the hits
      TrkStrawHitVector tshv;
      makeHits(shcol, tdef, tshv);
// Create the BaBar hit list, and fill it with these hits.  The BaBar list takes ownership
      std::vector<TrkHit*> thv;
      for(auto ihit = tshv.begin(); ihit != tshv.end(); ++ihit){
        thv.push_back(*ihit);
        if (_debug>2) { (*ihit)->print(std::cout); }
      }
// Find the wall and gas material description objects for these hits
      std::vector<DetIntersection> detinter;
      if(_matcorr)makeMaterials(tshv,tdef,detinter);
// create Kalman rep
      krep = new KalRep(tdef.helix(), thv, detinter, *this, tdef.particle());
      assert(krep != 0);
// initialize krep t0; eventually, this should be in the constructor, FIXME!!!
      double flt0 = tdef.helix().zFlight(0.0);
      krep->setT0(tdef.t0(),flt0);
// initialize history list
      krep->addHistory(TrkErrCode(),"KalFit creation");
// now fit
      TrkErrCode fitstat = fitTrack(krep,tshv);
      krep->addHistory(fitstat,"KalFit fit");
// extend the fit
      if(fitstat.success()){
	fitstat = extendFit(krep);
	krep->addHistory(fitstat,"KalFit extension");
      }
    }
  }

  void KalFit::addHits(KalRep* krep,const StrawHitCollection* shcol, std::vector<hitIndex> indices, double maxchi) {
  // fetcth the DetectorModel
   Mu2eDetectorModel const& detmodel{ art::ServiceHandle<BTrkHelper>()->detectorModel() };
// there must be a valid Kalman fit to add hits to
    if(krep != 0 && indices.size() > 0 && krep->fitStatus().success()){
      TrkStrawHitVector tshv;
      convert(krep->hitVector(),tshv);
      ConditionsHandle<TrackerCalibrations> tcal("ignored");
      const Tracker& tracker = getTrackerOrThrow();
      TrkStrawHitVector::iterator ihigh;
      TrkStrawHitVector::reverse_iterator ilow;
// use the reference trajectory, as that's what all the existing hits do
      const TrkDifPieceTraj* reftraj = krep->referenceTraj();
      for(unsigned iind=0;iind<indices.size(); ++iind){
        size_t istraw = indices[iind]._index;
        const StrawHit& strawhit(shcol->at(istraw));
        const Straw& straw = tracker.getStraw(strawhit.strawIndex());
// estimate  initial flightlength
        double hflt(0.0);
        TrkHelixUtils::findZFltlen(*reftraj,straw.getMidPoint().z(),hflt);
// find the bounding sites near this hit, and extrapolate to get the hit t0
        std::sort(tshv.begin(),tshv.end(),fltlencomp(_fdir.fitDirection()));
        findBoundingHits(tshv,hflt,ilow,ihigh);
        const TrkStrawHit* nearhit;
        if(ihigh != tshv.end())
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
        TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,hitt0,hflt,_herr.back(),_maxdriftpull);
        assert(trkhit != 0);
// allow the hit to update its own ambiguity for now: eventually we should get the resolver to do this, FIXME!!!
        trkhit->setAmbigUpdate(true);
        trkhit->setFlag(TrkStrawHit::addedHit);
// must be initialy active for KalRep to process correctly
        trkhit->setActivity(true);
// add the hit to the track
        krep->addHit(trkhit);
// check the raw residual: This call works because the HOT isn't yet processed as part of the fit.
        double chi = fabs(trkhit->residual()/trkhit->hitRms());
//if it's outside limits, deactivate the HOT
        if(chi > maxchi || !trkhit->physicalDrift(maxchi))
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
      TrkErrCode fitstat = fitIteration(krep,tshv,_herr.size()-1);
      krep->addHistory(fitstat,"AddHits");
    }
  }
//
  TrkErrCode KalFit::fitTrack(KalRep* krep,TrkStrawHitVector& tshv) {
    // loop over external hit errors, ambiguity assignment, t0 toleratnce
    TrkErrCode fitstat;
    for(size_t iherr=0;iherr < _herr.size(); ++iherr) {
      fitstat = fitIteration(krep,tshv,iherr);
      if(!fitstat.success())break;
    }
    if(_debug > 0) cout << fitstat << endl;
    return fitstat;
  }

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
  TrkErrCode KalFit::fitIteration(KalRep* krep,TrkStrawHitVector& tshv, size_t iter) {
    // update the external hit errors.  This isn't strictly necessary on the 1st iteration.
    for (auto itsh=tshv.begin();itsh!=tshv.end(); ++itsh){
      (*itsh)->setExtErr(_herr[iter]);
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
        updateT0(krep,tshv);
        changed |= fabs(krep->t0()._t0-oldt0) > _t0tol[iter];
        oldt0 = krep->t0()._t0;
      }
      // drop outliers
      if(_weedhits[iter]){
        changed |= weedHits(krep,tshv,iter);
      }
      // find missing materials
      if(_addmaterial[iter])
        changed |= addMaterial(krep) > 0;
    }
    if(_debug > 1)
      std::cout << "Fit iteration " << iter << " stopped after "
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

  void
  KalFit::makeHits(const StrawHitCollection* shcol, TrkDef const& tdef, TrkStrawHitVector& tshv ) {
    const Tracker& tracker = getTrackerOrThrow();
// compute the propagaion velocity
    double flt0 = tdef.helix().zFlight(0.0);
    double mom = TrkMomCalculator::vecMom(tdef.helix(),bField(),flt0).mag();
    double vflt = tdef.particle().beta(mom)*CLHEP::c_light;
    unsigned nind = tdef.strawHitIndices().size();
    for(unsigned iind=0;iind<nind;iind++){
      size_t istraw = tdef.strawHitIndices()[iind]._index;
      const StrawHit& strawhit(shcol->at(istraw));
      const Straw& straw = tracker.getStraw(strawhit.strawIndex());
      double fltlen = tdef.helix().zFlight(straw.getMidPoint().z());
    // estimate arrival time at the wire
      TrkT0 hitt0(tdef.t0());
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
      tshv.push_back(trkhit);
    }
 // sort the hits by flightlength
    std::sort(tshv.begin(),tshv.end(),fltlencomp(_fdir.fitDirection()));
  }

  void
  KalFit::makeMaterials(TrkStrawHitVector const& tshv, TrkDef const& tdef,std::vector<DetIntersection>& detinter) {
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
      strawelem->reIntersect(&tdef.helix(),strawinter);
      detinter.push_back(strawinter);
    }
  }

  unsigned KalFit::addMaterial(KalRep* krep) {
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
    for(auto plane : ttracker.getPlanes()){
    // crappy access to # of straws in a panel
      int nstraws = 2*plane.getPanel(0).getLayer(0).nStraws();
// get an approximate z position for this plane from the average position of the 1st and last straws
      Hep3Vector s0 = plane.getPanel(0).getLayer(0).getStraw(0).getMidPoint();
      // funky convention for straw numbering in a layer FIXME!!!!
      Hep3Vector sn = plane.getPanel(0).getLayer(1).getStraw(2*plane.getPanel(0).getLayer(1).nStraws()-1).getMidPoint();
      double pz = 0.5*(s0.z() + sn.z());
// find the transverse position at this z using the reference trajectory
      double flt = zFlight(krep,pz);
      HepPoint pos = krep->referenceTraj()->position(flt);
      Hep3Vector posv(pos.x(),pos.y(),pos.z());
// see if this position is in the active region.  Double the straw radius to be generous
      double rho = posv.perp();
      double rmin = s0.perp()-2*strawradius;
      double rmax = sn.perp()+2*strawradius;
      if(rho > rmin && rho < rmax){
  // loop over panels
        for(auto panel : plane.getPanels()){
      // get the straw direction for this panel
          Hep3Vector sdir = panel.getLayer(0).getStraw(0).getDirection();
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
            for(int is = max(0,istraw-2); is<min(nstraws-1,istraw+2); ++is){
            // must do this twice due to intrusion of layer on hierarchy FIXME!!!
              matstraws.insert(StrawFlight(panel.getLayer(0).getStraw(is).index(),flt));
              matstraws.insert(StrawFlight(panel.getLayer(1).getStraw(is).index(),flt));
              nadded += 2;
            }
          }
        }
      }
    }
// Now test if the Kalman rep hits these straws
    if(_debug>2)std::cout << "Found " << matstraws.size() << " unique possible straws " << " out of " << nadded << std::endl;
    for(auto strawflt : matstraws){
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
  KalFit::weedHits(KalRep* krep, TrkStrawHitVector& tshv,size_t iter) {
    // Loop over HoTs and find HoT with largest contribution to chi2.  If this value
    // is greater than some cut value, deactivate that HoT and reFit
    bool retval(false);
    double worst = -1.;
    TrkStrawHit* worsthit = 0;
    for (auto ihit=tshv.begin();ihit!=tshv.end(); ++ihit){
      TrkStrawHit* hit = *ihit;
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
      worsthit->setFlag(TrkStrawHit::weededHit);
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
        retval |= weedHits(krep,tshv,iter);
      }
    }
    return retval;
  }

  bool
  KalFit::unweedHits(KalRep* krep, double maxchi) {
    TrkStrawHitVector tshv;
    convert(krep->hitVector(),tshv);
    bool retval = unweedHits(krep,tshv,maxchi);
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
  KalFit::unweedHits(KalRep* krep, TrkStrawHitVector& tshv, double maxchi) {
    // Loop over inactive HoTs and find the one with the smallest contribution to chi2.  If this value
    // is less than some cut value, reactivate that HoT and reFit
    bool      retval(false);
    double    best = 1.e12;
// no need to cast
    TrkStrawHit* besthit = 0;
    for (auto ihit=tshv.begin();ihit!=tshv.end(); ++ihit){
      TrkStrawHit* hit = *ihit;
      if (!hit->isActive()) {
        double resid, residErr;
        if(hit->resid(resid, residErr, true)){
          double chival = fabs(resid/residErr);
  // test both for a good chisquared and for the drift radius to be physical
          if (chival < maxchi && hit->physicalDrift(maxchi) && chival < best) {
            best = chival;
            besthit = hit;
          }
        }
      }
    }
    if(0 != besthit){
      retval = true;
      besthit->setActivity(true);
      besthit->setFlag(TrkStrawHit::unweededHit);
      TrkErrCode fitstat = krep->fit();
      if (fitstat.success() && besthit->isActive() ) {
	krep->addHistory(fitstat, "HitUnWeed");
	// Recursively iterate
        retval |= unweedHits(krep,tshv,maxchi);
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
  KalFit::initT0(const StrawHitCollection* shcol,TrkDef& tdef) {
    TrkT0 t0 = tdef.t0();
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
      const StrawHit& strawhit(shcol->at(istraw));
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
    tdef.setT0(t0);
  }

  bool
  KalFit::updateT0(KalRep* krep,TrkStrawHitVector& tshv){
    using namespace boost::accumulators;
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
        for(auto ihit=tshv.begin(); ihit != tshv.end(); ihit++){
          TrkStrawHit* hit = *ihit;
          if(hit->isActive() && hit->hasResidual()){
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
            updateHitTimes(krep,tshv);
            retval = true;
          }
        }
      }
    }
    return retval;
  }

  void
  KalFit::updateHitTimes(KalRep* krep,TrkStrawHitVector& tshv) {
  // compute the time the track came closest to the wire for each hit, starting from t0 and working out.
  // this function allows for momentum change along the track.
  // find the bounding hits on either side of this
    std::sort(tshv.begin(),tshv.end(),fltlencomp(_fdir.fitDirection()));
    TrkStrawHitVector::iterator ihigh;
    TrkStrawHitVector::reverse_iterator ilow;
    findBoundingHits(tshv,krep->flt0(),ilow,ihigh);
    // reset all the hit times
    double hflt = krep->flt0();
    TrkT0 hitt0 = krep->t0();
    for(TrkStrawHitVector::iterator ihit= ihigh;ihit != tshv.end(); ++ihit){
      TrkStrawHit* hit = *ihit;
// particle momentum at this point, using the full fit
      double mom = krep->momentum(hit->fltLen()).mag();
// relativistic velocity from that
      double beta = krep->particleType().beta(mom);
// particle transit time to this hit from the reference
      double tflt = (hit->fltLen()-hflt)/(beta*CLHEP::c_light);
// update the time in the TrkT0 object
      hitt0._t0 += tflt;
      (*ihit)->updateHitT0(hitt0);
// update the reference flightlength
      hflt = hit->fltLen();
    }
// now the same, moving backwards
    hflt = krep->flt0();
    hitt0 = krep->t0();
    for(TrkStrawHitVector::reverse_iterator ihit= ilow;ihit != tshv.rend(); ++ihit){
      TrkStrawHit* hit = *ihit;
      double mom = krep->momentum(hit->fltLen()).mag();
      double beta = krep->particleType().beta(mom);
      double tflt = (hit->fltLen()-hflt)/(beta*CLHEP::c_light);
      hitt0._t0 += tflt;
      (*ihit)->updateHitT0(hitt0);
      hflt = hit->fltLen();
    }
  }

  void
  KalFit::findBoundingHits(TrkStrawHitVector& hits,double flt0,
    TrkStrawHitVector::reverse_iterator& ilow,
    TrkStrawHitVector::iterator& ihigh) {
    ilow = hits.rbegin();
    ihigh = hits.begin();
    while(ilow != hits.rend() && (*ilow)->fltLen() > flt0 )++ilow;
    while(ihigh != hits.end() && (*ihigh)->fltLen() < flt0 )++ihigh;
  }

  // this function belongs in TrkDifTraj, FIXME!!!!
  double KalFit::zFlight(KalRep* krep,double pz) {
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

// attempt to extend the fit to the specified location
  TrkErrCode KalFit::extendFit(KalRep* krep) {
    TrkErrCode retval;
    // find the downstream and upstream Z positions to extend to
    if(_exdown != noextension){
      double downz = extendZ(_exdown);
    // convert to flightlength using the fit trajectory
      double downflt = zFlight(krep,downz);
    // actually extend the track
      retval = krep->extendThrough(downflt);
    }
    // same for upstream extension
    if(retval.success() && _exup != noextension){
      double upz = extendZ(_exup);
      double upflt = zFlight(krep,upz);
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
      return cg->caloGeomInfo().envelopeZ1();
    }
    return retval;
  }
}



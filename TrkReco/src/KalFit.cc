//
// Class to perform BaBar Kalman fit
// Original author: Dave Brown LBNL 2012
//
//
#include "Offline/TrkReco/inc/KalFit.hh"
#include "Offline/TrkReco/inc/PanelAmbigResolver.hh"
#include "Offline/TrkReco/inc/HitAmbigResolver.hh"
#include "Offline/TrkReco/inc/FixedAmbigResolver.hh"
#include "Offline/TrkReco/inc/DoubletAmbigResolver.hh"
#include "Offline/TrkReco/inc/TrkUtilities.hh"
#include "Offline/Mu2eBTrk/inc/BaBarMu2eField.hh"
//geometry
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
// conditions
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
// data
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
// tracker
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
// BaBar
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/KalmanTrack/KalBend.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalMaterial.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
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
  struct fcomp {
    bool operator()(TrkHit* x, TrkHit* y) {
      return x->fltLen() < y->fltLen();
    }
  };

  // struct for finding materials
  struct StrawFlight {
    StrawId _id;  // straw being tested
    double _flt; // flight where trajectory comes near this straw
    // construct from pair
    StrawFlight(StrawId strawid, double flt) : _id(strawid), _flt(flt) {}
  };

  // comparison operators understand that the same straw could be hit twice, so the flight lengths need
  // to be similar befoew we consider these 'the same'
  struct StrawFlightComp {
    double _maxdiff; // maximum flight difference; below this, consider 2 intersections 'the same'
    StrawFlightComp(double maxdiff) : _maxdiff(maxdiff) {}
    bool operator () (StrawFlight const& a, StrawFlight const& b) const { return a._id < b._id ||
      ( a._id == b._id && a._flt < b._flt && fabs(a._flt-b._flt)>=_maxdiff);}
  };

  // construct from a parameter set
  KalFit::KalFit(fhicl::ParameterSet const& pset) :
    // KalFit parameters
    _debug(pset.get<int>("debugLevel",0)),
    _maxhitchi(pset.get<double>("maxhitchi",3.5)),
    _maxpull(pset.get<double>("maxPull",5)),
    _maxweed(pset.get<unsigned>("maxweed",10)),
    _maxweedtch(pset.get<unsigned>("maxweedtch",1)),
    // t0 parameters
    _useTrkCaloHit(pset.get<bool>("useTrkCaloHit")),
    _nCaloExtrapolSteps(pset.get<float>("nCaloExtrapolSteps", 100)),
    _caloHitErr(pset.get<double>("caloHitError")),
    _updatet0(pset.get<vector<bool>>("updateT0")),
    _t0tol(pset.get< vector<double> >("t0Tolerance")),
    _t0errfac(pset.get<double>("t0ErrorFactor",1.2)),
    _t0nsig(pset.get<double>("t0window",2.5)),
    _mindocatch(pset.get<double>("mindocatch",-50.)),
    _maxdocatch(pset.get<double>("maxdocatch", 50.)),
    _mindepthtch(pset.get<double>("mindepthtch",-50.)),
    _maxdepthtch(pset.get<double>("maxdepthtch",250.)),
    _maxtchdt(pset.get<double>("maxtchdt", 5.)),//ns
    _mintchenergy(pset.get<double>("mintchEnergy", 10.)),//MeV
    _mintchtrkpath(pset.get<double>("mintchTrkPath", 1.)),//mm
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
    _ttcalc            (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet())),
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

    _printUtils = new TrkPrintUtils(pset.get<fhicl::ParameterSet>("printUtils",fhicl::ParameterSet()));

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
      int Final(0);//      int Final = iter==niter-1 ? 1 : 0;
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

  void KalFit::setCaloGeom(){
    mu2e::GeomHandle<mu2e::Calorimeter> ch;

    _nCaloDisks = ch->nDisks();
    double      crystalLength = ch->caloInfo().getDouble("crystalZLength");
    for (unsigned i=0; i<_nCaloDisks; ++i){
      CLHEP::Hep3Vector pos(ch->disk(i).geomInfo().frontFaceCenter());
      pos = ch->geomUtil().mu2eToTracker(pos);

      _zmincalo[i] = (pos.z());
      _zmaxcalo[i] = (pos.z()+crystalLength);
      _rmincalo[i] = (ch->disk(i).geomInfo().innerEnvelopeR());
      _rmaxcalo[i] = (ch->disk(i).geomInfo().outerEnvelopeR());
    }
  }


  //-----------------------------------------------------------------------------
  // create the track (KalRep) from a track seed
  //-----------------------------------------------------------------------------
  void KalFit::makeTrack(StrawResponse::cptr_t srep,
      Mu2eDetector::cptr_t detmodel,
      KalFitData& kalData){

    // test if fitable
    if(fitable(*kalData.kalSeed)){
      // find the segment at the 0 flight
      double flt0 = kalData.kalSeed->flt0();
      auto kseg = kalData.kalSeed->nearestSegmentFlt(flt0);
      if(kseg->fmin() > kseg->localFlt(flt0) ||
          kseg->fmax() < kseg->localFlt(flt0) ){
        std::cout << "FitType: "<< kalData.fitType<<", number 0f segments = "<<kalData.kalSeed->segments().size()
          <<", Helix segment range doesn't cover flt0 = " << flt0 << std::endl;
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
      kalData.helixTraj = &htraj;
      // create the hits
      TrkStrawHitVector tshv;
      makeTrkStrawHits(srep,kalData, tshv);

      // Find the wall and gas material description objects for these hits
      std::vector<DetIntersection> detinter;
      if(_matcorr)makeMaterials(detmodel, tshv,*kalData.helixTraj,detinter);
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
        makeTrkCaloHit(kalData, tch);
        if (tch != 0) thv.push_back(tch);
      }

      TrkT0 t0(kalData.kalSeed->t0());
      // create Kalman rep
      kalData.krep = new KalRep(htraj, thv, detinter, *this, TrkParticle(TrkParticle::type(kalData.kalSeed->particle())), t0, flt0);
      assert(kalData.krep != 0);

      if (_debug > 0) {
        char msg[100];
        sprintf(msg,"makeTrack_001 annealing step: %2i",_annealingStep);
        _printUtils->printTrack(kalData.event,kalData.krep,"banner+data+hits",msg);
      }

      // initialize history list
      kalData.krep->addHistory(TrkErrCode(),"KalFit creation");
      // now fit
      TrkErrCode fitstat = fitTrack(detmodel,kalData);
      kalData.krep->addHistory(fitstat,"KalFit fit");
      // extend the fit
      if(fitstat.success()){
        fitstat = extendFit(kalData.krep);
        kalData.krep->addHistory(fitstat,"KalFit extension");
      }
    }
  }

  void KalFit::addHits(StrawResponse::cptr_t srep, Mu2eDetector::cptr_t detmodel,
      KalFitData&kalData, double maxchi) {
    //2017-05-02: Gianipez. In this function inten
    // there must be a valid Kalman fit to add hits to
    KalRep* krep = kalData.krep;

    if(kalData.krep != 0 && kalData.missingHits.size() > 0 && krep->fitStatus().success()){
      TrkHitVector::iterator ihigh;
      TrkHitVector::reverse_iterator ilow;
      // use the reference trajectory, as that's what all the existing hits do
      const TrkDifPieceTraj* reftraj = krep->referenceTraj();
      for(unsigned iind=0;iind<kalData.missingHits.size(); ++iind){
        size_t istraw = kalData.missingHits[iind].index;
        const ComboHit& strawhit(kalData.chcol->at(istraw));
        const Straw& straw = _tracker->getStraw(strawhit.strawId());
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

        //can I just do krep_hitt0()???; this whole block (next 5 lines of code) will disappear!
        HitT0  hitt0 = krep_hitT0(krep, nearhit);//nearhit->hitT0();
        double mom   = krep->momentum(nearhit->fltLen()).mag();
        double beta  = krep->particleType().beta(mom);
        double tflt  = (hflt-nearhit->fltLen())/(beta*CLHEP::c_light);
        // update the time in the HitT0 object
        hitt0._t0 += tflt;
        // create the hit object.  Assume we're at the last iteration over added error
        TrkStrawHit* trkhit = new TrkStrawHit(srep,strawhit,*_tracker,istraw,hitt0,hflt,
            _maxpull,_strHitW );
        assert(trkhit != 0);
        trkhit->setTemperature(_herr.back()); // give this hit the final annealing temperature
        trkhit->setFlag(TrkHit::addedHit);
        // guess the ambiguity form the sign of the doca
        int iambig;
        if (kalData.missingHits[iind].doca > 0) iambig =  1;
        else                                    iambig = -1;
        // can set ambiguity only for deactivated hit
        trkhit->setActivity(false);
        trkhit->setAmbig(iambig);
        // must be initialy active for KalRep to process correctly
        trkhit->setActivity(true);
        // set the hit ambiguity.  This is a preliminary value before using the official ambig resolver
        TrkPoca poca(krep->traj(),hflt,*trkhit->hitTraj(),0.0);
        int newamb = poca.doca() > 0 ? 1 : -1;
        trkhit->setAmbig(newamb);
        // add the hit to the track
        krep->addHit(trkhit);
        // check the raw residual: This call works because the HOT isn't yet processed as part of the fit.
        double chi = fabs(trkhit->residual()/trkhit->hitRms());
        //if it's outside limits, deactivate the HOT
        if(chi > maxchi || (!trkhit->isPhysical(maxchi)))
          trkhit->setActivity(false);
        // find the DetElem associated this straw
        const DetStrawElem* strawelem = detmodel->strawElem(trkhit->straw());
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
      TrkErrCode fitstat = fitIteration(detmodel,kalData,_herr.size()-1);
      krep->addHistory(fitstat,"AddHits");
    }
  }
  //
  TrkErrCode KalFit::fitTrack(Mu2eDetector::cptr_t detmodel, KalFitData&kalData) {
    // loop over external hit errors, ambiguity assignment, t0 toleratnce
    TrkErrCode fitstat;
    for(size_t iherr=0;iherr < _herr.size(); ++iherr) {
      fitstat = fitIteration(detmodel,kalData,iherr);
      if(_debug > 0) {
        cout << "Iteration " << iherr
          << " NDOF = " << kalData.krep->nDof()
          << " Fit Status = " <<  fitstat << endl;

        char msg[200];
        sprintf(msg,"KalFit::fitTrack Iteration = %2li success = %i",iherr,fitstat.success());
        _printUtils->printTrack(kalData.event,kalData.krep,"banner+data+hits",msg);
      }
      if(!fitstat.success())break;
    }
    return fitstat;
  }

  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  TrkErrCode KalFit::fitIteration(Mu2eDetector::cptr_t detmodel,
      KalFitData&kalData, int iter) {

    if (iter == -1) iter =  _herr.size()-1;
    _annealingStep = iter;//used in the printHits routine

    // update the external hit errors.  This isn't strictly necessary on the 1st iteration.
    TrkHitVector* thv   = &(kalData.krep->hitVector());
    for (auto itsh=thv->begin();itsh!=thv->end(); ++itsh){
      (*itsh)->setTemperature(_herr[iter]);
    }

    // update t0, and propagate it to the hits
    double oldt0 = kalData.krep->t0()._t0;
    unsigned niter(0);
    bool changed(true);
    TrkErrCode retval = TrkErrCode::succeed;

    KalRep* krep =  kalData.krep;
    bool    flagMaterialAdded(false);

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

      //2019-04-26: Giani&Dave; we want to do the weeding before the call to ::updateT0
      // drop outliers
      if(_weedhits[iter]){
        kalData.nweediter = 0;
        changed |= weedHits(kalData,iter);
        changed |= unweedBestHit(kalData,_maxhitchi);
      }
      // updates
      if(_updatet0[iter]){
        updateT0(kalData, iter);
        changed |= fabs(krep->t0()._t0-oldt0) > _t0tol[iter];
      }
      // find missing materials
      unsigned nmat(0);
      if(_addmaterial[iter]){
        nmat = addMaterial(detmodel,krep);
        changed |= nmat>0;
        if (!flagMaterialAdded) flagMaterialAdded=true;
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
    KalFit::fitable(KalSeed const& kseed){
      return kseed.segments().size() > 0 && kseed.hits().size() >= _minnstraws;
    }

  void
    KalFit::makeTrkStrawHits(StrawResponse::cptr_t srep,
        KalFitData& kalData, TrkStrawHitVector& tshv ) {

      std::vector<TrkStrawHitSeed>const hseeds = kalData.kalSeed->hits();
      HelixTraj const htraj = *kalData.helixTraj;
      // compute particle velocity to
      for(auto ths : hseeds ){
        // create a TrkStrawHit from this seed.
        size_t index = ths.index();
        const ComboHit& strawhit(kalData.chcol->at(index));
        TrkStrawHit* trkhit = new TrkStrawHit(srep,strawhit,*_tracker,ths.index(),ths.t0(),ths.trkLen(),
            _maxpull,_strHitW);
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
    KalFit::makeTrkCaloHit  (KalFitData& kalData, TrkCaloHit *&tch){
      art::Ptr<CaloCluster> const& calo = kalData.kalSeed->caloCluster();
      if (calo.isNonnull()){
        mu2e::GeomHandle<mu2e::Calorimeter> ch;
        Hep3Vector cog = ch->geomUtil().mu2eToTracker(ch->geomUtil().diskFFToMu2e( calo->diskID(), calo->cog3Vector()));
        if(_debug > 0){
          std::cout << "Cluster COG (disk) " << calo->cog3Vector() << std::endl
            << "Cluster COG (Mu2e) " << ch->geomUtil().diskFFToMu2e( calo->diskID(), calo->cog3Vector()) << std::endl
            <<" Cluster COG (Det ) " << cog << std::endl;
        }
        double      crystalLength = ch->caloInfo().getDouble("crystalZLength");
        // estimate fltlen from pitch; take the last segment
        HelixVal const& hval = kalData.kalSeed->segments().back().helix();
        double mom = kalData.kalSeed->segments().back().mom();
        TrkParticle tpart(TrkParticle::type(kalData.kalSeed->particle()));
        double beta = tpart.beta(mom);
        double td = hval.tanDip();
        double sd = td/sqrt(1.0+td*td);
        double fltlen = (cog.z()- hval.z0() + 0.5*crystalLength)/sd;// - kalData.kalSeed->flt0();
        // t0 represents the time the particle reached the sensor; estimate that
        HitT0 ht0 = kalData.kalSeed->t0(); // start with the track t0
        double   flt0 = kalData.kalSeed->flt0();
        double tflt = (fltlen -flt0)/(beta*CLHEP::c_light);
        ht0._t0 += tflt;
        ht0._t0err = _ttcalc.caloClusterTimeErr();
        Hep3Vector clusterAxis = Hep3Vector(0, 0, 1);//FIXME! should come from crystal
        tch = new TrkCaloHit(*kalData.kalSeed->caloCluster().get(), cog, crystalLength, clusterAxis, ht0, fltlen, _calHitW, _caloHitErr, _ttcalc.caloClusterTimeErr(), _ttcalc.trkToCaloTimeOffset());
      }
    }


  void
    KalFit::makeMaterials( Mu2eDetector::cptr_t detmodel,
        TrkStrawHitVector const& tshv, HelixTraj const& htraj,
        std::vector<DetIntersection>& detinter) {
      // loop over strawhits and extract the straws
      for (auto trkhit : tshv) {
        // find the DetElem associated this straw
        const DetStrawElem* strawelem = detmodel->strawElem(trkhit->straw());
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

  unsigned KalFit::addMaterial(Mu2eDetector::cptr_t detmodel, KalRep* krep) {
    unsigned retval(0);
    // Tracker geometry
    const Tracker& tracker = *_tracker;
    // general properties; these should be computed once/job and stored FIXME!
    double strawradius = tracker.strawOuterRadius();
    auto const& frontplane = tracker.planes().front();
    auto const& firstpanel = frontplane.getPanel(0);
    auto const& innerstraw = firstpanel.getStraw(0);
    auto const& outerstraw = firstpanel.getStraw(StrawId::_nstraws-1);
    // compute limits: add some buffer for the finite size of the straw
    auto DStoP = firstpanel.dsToPanel();
    auto innerstraw_origin = DStoP*innerstraw.origin();
    auto outerstraw_origin = DStoP*outerstraw.origin();
    double ymin = innerstraw_origin.y() - strawradius;
    double ymax = outerstraw_origin.y() + strawradius;
    double umax = innerstraw.halfLength() + strawradius;
    // use the outermost straw end to set the max hit radius
    double rmax = outerstraw.wireEnd(StrawEnd::cal).mag() + strawradius;
    double spitch = (StrawId::_nstraws-1)/(ymax-ymin);
    // storage of potential straws
    StrawFlightComp strawcomp(_maxmatfltdiff);
    std::set<StrawFlight,StrawFlightComp> matstraws(strawcomp);
    // loop
    unsigned nadded(0);
    for(auto const& plane : tracker.planes()){
      if(_tracker->planeExists(plane.id())) {
        // get an approximate z position for this plane from the average position of the 1st and last straws
        auto s0 = plane.origin();
        // find the track position at this z using the reference trajectory
        double flt = krep->referenceTraj()->zFlight(s0.z());
        HepPoint pos = krep->referenceTraj()->position(flt);
        Hep3Vector posv(pos.x(),pos.y(),pos.z());
        // loop over panels
        for(auto panel_p : plane.panels()){
          auto const& panel = *panel_p;
          // convert track position into panel coordinates
          auto DStoP = panel.dsToPanel();
          auto pposv = DStoP*posv;
          // see if this point is roughly in the active region of this panel.  Use the z possition as a buffer, to
          // account for the test being performed at the plane center.  Note the radius cut is made in the Mu2e coordinate system
          // this is not a bug!
          double pbuff = fabs(pposv.z());
          if(pposv.y() > ymin - pbuff && pposv.y() < ymax + pbuff && fabs(pposv.x()) < umax && posv.perp() < rmax + pbuff) {
            if(_debug>2)std::cout << "position " << pposv << " in rough acceptance " << std::endl;
            // translate the y position into a rough straw number
            int istraw = (int)rint( (pposv.y()-ymin)*spitch);
            // take a few straws around this.  This value should be configurable FIXME!
            for(int is = max(0,istraw-3); is<min(StrawId::_nstraws-1,istraw+3); ++is){
              if(_debug>3)std::cout << "Adding Straw " << is << " in panel " << panel.id() << std::endl;
              matstraws.insert(StrawFlight(panel.getStraw(is).id(),flt));
              ++nadded;
            }
          } // acceptance
        } // panels
      } // plane exists
    }  // planes
    // Now test if the Kalman rep hits these straws
    if(_debug>2)std::cout << "Found " << matstraws.size() << " unique possible straws " << " out of " << nadded << std::endl;
    unsigned nfound(0);
    for(auto const& strawflt : matstraws){
      const DetStrawElem* strawelem = detmodel->strawElem(strawflt._id);
      DetIntersection strawinter;
      strawinter.delem = strawelem;
      strawinter.pathlen = strawflt._flt;
      if(strawelem->reIntersect(krep->referenceTraj(),strawinter)){
        // If the rep already has a material site for this element, skip it
        std::vector<const KalMaterial*> kmats;
        krep->findMaterialSites(strawelem,kmats);
        if(_debug>2)std::cout << "found intersection with straw " << strawelem->straw()->id() << " with "
          << kmats.size() << " materials " << std::endl;
        // test material isn't on the track
        bool hasmat(false);
        for(auto kmat : kmats ){
          const DetStrawElem* kelem = dynamic_cast<const DetStrawElem*>(kmat->detIntersection().delem);
          if(kelem != 0){
            StrawFlight ksflt(kelem->straw()->id(),kmat->globalLength());
            if(_debug>2)std::cout << " comparing flights " << kmat->globalLength() << " and " << strawflt._flt << std::endl;
            if(!strawcomp.operator()(strawflt,ksflt)){
              if(_debug>2)std::cout << "operator returned false!!" << std::endl;
              // this straw is already on the track: stop
              hasmat = true;
              nfound++;
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
    if(_debug>1)std::cout << "Added " << retval << " new material sites; found " << nfound << " intersections out of " << krep->nActive() << " active hits " << std::endl;
    return retval;
  }

  bool
    KalFit::weedHits(KalFitData& kalData, int iter) {
      // Loop over HoTs and find HoT with largest contribution to chi2.  If this value
      // is greater than some cut value, deactivate that HoT and reFit
      KalRep* krep = kalData.krep;
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

      if (iter == -1) iter = _ambigresolver.size()-1;

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

        if (_debug > 0) {
          char msg[200];
          sprintf(msg,"KalFit::weedHits Iteration = %2i success = %i",iter,fitstat.success());
          _printUtils->printTrack(kalData.event,kalData.krep,"banner+data+hits",msg);
        }

        // Recursively iterate
        kalData.nweediter++;
        if (fitstat.success() && kalData.nweediter < _maxweed) {
          retval |= weedHits(kalData,iter);
        }
      }
      return retval;
    }

  bool
    KalFit::unweedHits(KalFitData& kalData, double maxchi) {
      bool retval = unweedBestHit(kalData, maxchi);
      // if any hits were added, re-analyze ambiguity
      if (retval && _resolveAfterWeeding) {
        KalRep* krep = kalData.krep;
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

          if (_debug > 0) {
            char msg[200];
            sprintf(msg,"KalFit::unweedHits Iteration = %2i success = %i",last,fitstat.success());
            _printUtils->printTrack(kalData.event,kalData.krep,"banner+data+hits",msg);
          }
        }
      }
      return retval;
    }

  bool
    KalFit::unweedBestHit(KalFitData& kalData, double maxchi) {
      // Loop over inactive HoTs and find the one with the smallest contribution to chi2.  If this value
      // is less than some cut value, reactivate that HoT and reFit
      KalRep*   krep = kalData.krep;
      bool      retval(false);
      double    best = 1.e12;
      // no need to cast
      TrkHit*   besthit = 0;
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
        int oldndof = krep->nDof();
        TrkErrCode fitstat = krep->fit();
        krep->addHistory(fitstat, "HitUnWeed");
        if (fitstat.success() && besthit->isActive() && krep->nDof() > oldndof) {
          // Recursively iterate
          retval |= unweedBestHit(kalData, maxchi);
        }
      }
      return retval;
    }

  //--------------------------------------------------------------------------------
  // This function uses the KalRep for searching for the calorimeter disk where
  // the track is supposed to impact
  //--------------------------------------------------------------------------------
  void
    KalFit::findCaloDiskFromTrack(KalFitData& kalData, int& trkToCaloDiskId, double& caloFlt){
      KalRep*krep = kalData.krep;
      const TrkDifPieceTraj* reftraj = krep->referenceTraj();
      float  zExtrapolStep  = (_zmaxcalo[0] - _zmincalo[0])/(float)_nCaloExtrapolSteps;

      //initialize the output values
      trkToCaloDiskId = -1;
      caloFlt         = 0;

      float  fltIn(0), fltOut(0);
      for (unsigned i=0; i<_nCaloDisks; ++i){
        for (unsigned j=0; j<_nCaloExtrapolSteps; ++j){
          double flt(0);
          double zCaloExtrap = _zmincalo[i] + j*zExtrapolStep;
          TrkHelixUtils::findZFltlen(*reftraj, zCaloExtrap, flt);
          HepPoint pos    = krep->referenceTraj()->position(flt);
          float    radius = sqrt(pos.x()*pos.x() + pos.y()*pos.y());
          if ( (radius >= _rmincalo[i]) && (radius <= _rmaxcalo[i])){
            if (trkToCaloDiskId<0) {
              trkToCaloDiskId = i;
              fltIn  = flt;
            }else {
              fltOut = flt;
            }
          }
        }//end loop over the z-steps
        if (trkToCaloDiskId>=0) break;
      }//end loop over tghe disks
      if (trkToCaloDiskId >= 0)  caloFlt = fltOut - fltIn;
    }


  //--------------------------------------------------------------------------------
  // This function loops over the CaloClusterCollection to search for a Cluster
  // compatible with the Track. If the Cluster energy was below the threshold,
  // no Cluster was added in the TimeClusterFinder module
  //--------------------------------------------------------------------------------
  int
    KalFit::addTrkCaloHit( Mu2eDetector::cptr_t detmodel, KalFitData& kalData) {
      int retval(-1);
      //extrapolate the track to the calorimeter region
      //to understand on which disk the track is supposed to impact
      int      trkToCaloDiskId(-1);
      double   trkInCaloFlt(0);
      findCaloDiskFromTrack(kalData, trkToCaloDiskId, trkInCaloFlt);
      if (trkToCaloDiskId >= 0 &&  //the Track doesn't intercept the calorimeter
          fabs(trkInCaloFlt) > _mintchtrkpath) { //FIX ME! should we check the second disk in case the track-path in the first is too small?
        KalRep*  krep = kalData.krep;
        double   minFOM(1e10);
        const CaloCluster*cl(0);
        std::unique_ptr<TrkCaloHit> tchFinal;
        mu2e::GeomHandle<mu2e::Calorimeter> ch;
        double   crystalLength = ch->caloInfo().getDouble("crystalZLength");

        unsigned nClusters = kalData.caloClusterCol->size();
        const TrkDifPieceTraj* reftraj = krep->referenceTraj();
        double   flt0 = krep->flt0();
        double   tflt(0), flt(0);
        if (trkToCaloDiskId>=0){
          //evaluate the flight length at the z of the calorimeter cluster + half crystallength
          TrkHelixUtils::findZFltlen(*reftraj, (_zmincalo[trkToCaloDiskId]+0.5*crystalLength),flt);
          //evaluate the transittime using the full trajectory
          tflt = krep->t0()._t0 + krep->transitTime(flt0, flt);
        }

        for (unsigned icc=0; icc<nClusters; ++icc){
          cl    = &kalData.caloClusterCol->at(icc);
          if (cl->diskID() != trkToCaloDiskId ||
              cl->energyDep() < _mintchenergy) continue;
          // double      hflt(0.0);
          Hep3Vector cog = ch->geomUtil().mu2eToTracker(ch->geomUtil().diskFFToMu2e( cl->diskID(), cl->cog3Vector()));
          double      dt = cl->time() + _ttcalc.trkToCaloTimeOffset() - tflt;

          //check the compatibility of the track and time within a given time window
          if (fabs(dt) > _maxtchdt)        continue;

          HitT0 ht0;
          ht0._t0 = tflt;
          // initial error can't be better than the input error
          ht0._t0err = _ttcalc.caloClusterTimeErr();

          Hep3Vector  clusterAxis   = Hep3Vector(0, 0, 1);//FIXME! should come from crystal
          // construct a temporary TrkCaloHit.  This is just to be able to call POCA
          TrkLineTraj hitTraj(HepPoint(cog.x(), cog.y(), cog.z()),
              clusterAxis, 0.0, crystalLength);
          //evaluate the doca
          TrkPoca poca(krep->traj(),flt,hitTraj,0.5*crystalLength);
          double doca = poca.doca();
          double depth = poca.flt2();
          if( doca  > _mindocatch  && doca  < _maxdocatch &&
              depth > _mindepthtch && depth < _maxdepthtch &&
              fabs(doca) < minFOM) {
            tchFinal.reset(new TrkCaloHit(*cl, cog, crystalLength, clusterAxis,
                  ht0, poca.flt1(),
                  _calHitW, _caloHitErr,
                  _ttcalc.caloClusterTimeErr(), _ttcalc.trkToCaloTimeOffset()));
            minFOM   = doca; // this should be some combination of energy, DOCA, etc FIXME!
            retval = icc;
          }
        }

        if (tchFinal != 0) {

          //add the TrkCaloHit
          krep->addHit(tchFinal.release());

          TrkErrCode fitstat = fitIteration(detmodel,kalData,_herr.size()-1);
          krep->addHistory(fitstat,"AddHits");
        }
      }

      return retval;

    }

  void
    KalFit::fillTchDiag(KalFitData& kalData){
      KalRep* krep = kalData.krep;
      TrkHitVector *thv      = &(krep->hitVector());

      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      double   crystalLength = ch->caloInfo().getDouble("crystalZLength");

      const TrkDifPieceTraj* reftraj = krep->referenceTraj();
      double   flt0 = krep->flt0();

      //evaluate the track-path length in the calorimeter
      int      trkToCaloDiskId(-1);
      double   trkInCaloFlt(0);
      findCaloDiskFromTrack(kalData, trkToCaloDiskId, trkInCaloFlt);

      for (auto ihit=thv->begin();ihit!=thv->end(); ++ihit){
        TrkCaloHit*hit = dynamic_cast<TrkCaloHit*>(*ihit);
        if (hit == 0)     continue;
        if (hit->isActive()) {
          //evaluate the flight length at the z of the calorimeter cluster + half crystallength
          unsigned    diskId = hit->caloCluster().diskID();
          double      flt(0);
          TrkHelixUtils::findZFltlen(*reftraj, (_zmincalo[diskId]+0.5*crystalLength),flt);
          //evaluate the transittime using the full trajectory
          double      tflt = krep->t0()._t0 + krep->transitTime(flt0, flt);
          double      dt   = hit->caloCluster().time() + _ttcalc.trkToCaloTimeOffset() - tflt;

          kalData.diag.diskId   = diskId;
          kalData.diag.depth    = hit->hitLen();
          kalData.diag.dt       = dt;
          kalData.diag.trkPath  = trkToCaloDiskId == (int)diskId ? trkInCaloFlt : -9999;
          kalData.diag.energy   = hit->caloCluster().energyDep();
          kalData.diag.doca     = hit->poca().doca();
        }
      }
    }

  bool
    KalFit::weedTrkCaloHit(KalFitData& kalData, int iter) {
      // check if the TrkCaloHit residuals is within a given limit
      KalRep* krep = kalData.krep;
      bool    retval(false);
      //    double  worst = -1.;
      TrkCaloHit   *worsthit = 0;
      TrkHitVector *thv      = &(krep->hitVector());

      for (auto ihit=thv->begin();ihit!=thv->end(); ++ihit){
        TrkCaloHit*hit = dynamic_cast<TrkCaloHit*>(*ihit);
        if (hit == 0)     continue;
        if (hit->isActive()) {
          //evaluate the doca
          double  doca  = hit->poca().doca();
          //evaluate the crystal depth
          double  depth = hit->hitLen();
          if( (doca  < _mindocatch ) || (doca  > _maxdocatch) ||
              (depth < _mindepthtch) || (depth > _maxdepthtch)) {
            worsthit = hit;
            break;
          }
        }
      }

      if (iter == -1) iter = _ambigresolver.size()-1;

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

        if (_debug > 0) {
          char msg[200];
          sprintf(msg,"KalFit::weedTrkCaloHit Iteration = %2i success = %i",iter,fitstat.success());
          _printUtils->printTrack(kalData.event,kalData.krep,"banner+data+hits",msg);
        }

      }

      kalData.nweedtchiter++;
      return retval;
    }

  BField const&
    KalFit::bField() const {
      if(_bfield == 0){
        if(_fieldcorr){
          // create a wrapper around the mu2e field
          _bfield = new BaBarMu2eField();
        } else {
          // create a fixed field using the field at the tracker origin
          GeomHandle<BFieldManager> bfmgr;
          GeomHandle<DetectorSystem> det;
          // change coordinates to mu2e
          CLHEP::Hep3Vector vpoint(0,0,0);
          CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(vpoint);
          CLHEP::Hep3Vector field = bfmgr->getBField(vpoint_mu2e);
          _bfield=new BFieldFixed(CLHEP::Hep3Vector(0.0,0.0,field.z()));
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


  bool
    KalFit::updateT0(KalFitData& kalData, int iter){
      KalRep* krep = kalData.krep;
      using namespace boost::accumulators;
      TrkHitVector *thv = &(krep->hitVector());
      bool retval(false);
      // need to have a valid fit
      if(krep->fitValid()){
        // find the global fltlen associated with z=0.
        double flt0(0.0);
        bool converged = TrkHelixUtils::findZFltlen(krep->traj(),0.0,flt0);
        if(converged){
          std::vector<double> hitt0, hitt0err; // store t0, to allow outlyer removal
          size_t nhits = krep->hitVector().size();
          hitt0.reserve(nhits);
          hitt0err.reserve(nhits);

          if (iter == -1) iter = _herr.size()-1;
          accumulator_set<double,stats<tag::weighted_variance>,double > wmean;

          // loop over the hits and accumulate t0
          for(auto ihit=thv->begin(); ihit != thv->end(); ihit++){
            TrkHit*      hit   = *ihit;
            bool         trkShAmbigOK(true);
            if (_ambigstrategy[iter] != 0) {
              TrkStrawHit* trkSh = dynamic_cast<TrkStrawHit*>(*ihit);
              if (trkSh !=0){
                if (trkSh->ambig() == 0)
                  trkShAmbigOK = false;
              }
            }
            if(hit->isActive() && trkShAmbigOK) {
              HitT0 st0;
              //      if (hit->signalPropagationTime(st0 )){
              if (hit_time(hit, st0)){//2019-04-22: this function will become TrkHit::time(HitT0& hitT0) in the fututre development. FIXME!
                // subtracting hitT0 makes this WRT the previous track t0
                double    dtHitToTrack = st0._t0 - krep_hitT0(krep, hit)._t0;//FIXME! KalRep should own this function
                wmean(dtHitToTrack, weight=1.0/(st0._t0err*st0._t0err));

              }
            }
            }

            TrkT0 t0; // null t0; this will be the change in t0 from this update
            unsigned nused = extract_result<tag::count>(wmean);
            t0._t0    = extract_result<tag::weighted_mean>(wmean);
            t0._t0err = sqrt(extract_result<tag::weighted_variance>(wmean)/nused);

            // put in t0 from the track.
            t0._t0 += krep->t0()._t0;
            krep->setT0(t0,flt0);
            updateHitTimes(krep);
            retval = true;

          }
        }
        return retval;
      }

      void KalFit::updateHitTimes(KalRep* krep) {
        // compute the time the track came closest to the sensor for each hit, starting from t0 and working out.
        // this function allows for momentum change along the track.
        // find the bounding hits on either side of this
        TrkHitVector *thv = &(krep->hitVector());
        std::sort(thv->begin(),thv->end(),fcomp());
        TrkHitVector::iterator ihigh;
        TrkHitVector::reverse_iterator ilow;
        findBoundingHits(krep, krep->flt0(),ilow,ihigh);
        // reset all the hit times
        double flt0 = krep->flt0();
        HitT0 hitt0 = krep->t0();
        //GIANIPEZ 2019-04-26: update the following loops in the bottom using the function kres_hitt0//FIXME!
        for(TrkHitVector::iterator ihit= ihigh;ihit != thv->end(); ++ihit){
          TrkHit* hit = *ihit;
          double flt1 = hit->fltLen();
          // particle transit time to this hit from the reference
          double tflt = krep->transitTime(flt0, flt1);
          // update the time in the TrkT0 object
          hitt0._t0 += tflt;
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
          (*ihit)->setHitT0(hitt0);
          flt0 = flt1;
        }
      }

      void KalFit::findBoundingHits(KalRep* krep,double flt0,
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
          return cg->caloInfo().getDouble("envelopeZ1");
        }
        return retval;
      }

      HitT0  KalFit::krep_hitT0(KalRep*krep, const TrkHit*hit){
        HitT0  t0;
        double flt0 = krep->flt0();
        double flt1 = hit->fltLen();
        // particle transit time to this hit from the reference
        double tflt = krep->transitTime(flt0, flt1);
        t0._t0    = krep->t0()._t0 + tflt;
        t0._t0err = krep->t0()._t0err;//the error contribution from the *tflt* term is neglected. In the last fit iteration we might be considering adding a contribution from *tlft*
        return t0;
      }


      //--------------------------------------------------------------------------------
      //
      //--------------------------------------------------------------------------------
      bool   KalFit::hit_time(TrkHit*hit, HitT0& hitT0){
        TrkT0 st0;
        if (hit->signalPropagationTime(st0)){
          hitT0._t0    = hit->time() - st0._t0;
          hitT0._t0err = st0._t0err;
          return true;
        }else {
          return false;
        }
      }

    }

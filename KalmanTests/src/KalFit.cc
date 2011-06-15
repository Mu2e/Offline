//
// Class to perform BaBar Kalman fit
//
// $Id: KalFit.cc,v 1.4 2011/06/15 17:52:47 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/06/15 17:52:47 $
//

// the following has to come before other BaBar includes
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/KalFit.hh"
//geometry
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "KalmanTrack/KalHit.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/HelixTraj.hh"
#include "TrkBase/TrkHotListFull.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "TrkBase/TrkPoca.hh"
#include "BaBar/ErrLog.hh"
#include "BField/BFieldFixed.hh"
#include "DetectorModel/DetIntersection.hh"
#include "DetectorModel/DetMaterial.hh"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>


using namespace std; 

namespace mu2e 
{
// statics
// convert speed of light in mm/nsec
  const double KalFit::_vlight = CLHEP::c_light;
// drift velocity should come from a service FIXME!!!  
  const double KalFit::_vdrift = 0.05; // 50 um/nsec
// Material information, BaBar style
  MatDBInfo* KalFit::_matdbinfo(new MatDBInfo);
// comparison functor for ordering hits  
  struct fltlencomp : public binary_function<TrkStrawHit*, TrkStrawHit*, bool> {
          bool operator()(TrkStrawHit* x, TrkStrawHit* y) { return x->fltLen() < y->fltLen(); }
  };
  
// construct from a parameter set  
  KalFit::KalFit(fhicl::ParameterSet const& pset) :
    _wallelem(wall,_matdbinfo),_gaselem(gas,_matdbinfo),
    _debug(pset.get<int>("debugLevel",0)),
    _fieldcorr(pset.get<bool>("fieldCorrections",false)),
    _material(pset.get<bool>("material",false)),
    _ambigflip(pset.get<bool>("ambigflip",false)),
    _weedhits(pset.get<bool>("weedhits",true)),
    _updatet0(pset.get<bool>("updateT0",true)),
    _removefailed(pset.get<bool>("RemoveFailedFits",true)),
    _t0tol(pset.get<double>("t0Tolerance",1.0)),
    _maxhitchi(pset.get<double>("maxhitchi",5.0)),
    _maxiter(pset.get<unsigned>("maxiter",3)),
    _minnstraws(pset.get<unsigned>("minnstraws",20)),
    _maxweed(pset.get<unsigned>("maxweed",10))
  {
      _kalcon = new KalContext;
      _kalcon->setBendSites(_fieldcorr);
      _kalcon->setMaterialSites(_material);
      _kalcon->setForbidAmbigFlips(_ambigflip); //false: free left-rigth ambiguity, true: will be fixed from sim
      _kalcon->setMaxIterations(_maxiter);
      // these are currently fixed, they should be set as parameters and re-optimized FIXME!!!!
      _kalcon->setMaxIntersections(0);
      _kalcon->setMaxDMom(10);
      _kalcon->setSmearFactor(1e6);
      _kalcon->setMinDOF(20,TrkEnums::bothView);
      _kalcon->setMinDOF(20,TrkEnums::xyView);
      _kalcon->setMinDOF(0,TrkEnums::zView);
      _kalcon->setIntersectionTolerance(100);
      _kalcon->setMaxMomDiff(1.0); // 1 MeV
      _kalcon->setTrajBuffer(0.01); // 10um
      _kalcon->setMinGap(0.0); // no minimum separation between sites
      _kalcon->setDefaultType(PdtPid::electron); // by default, fit electrons
  }

  KalFit::~KalFit(){}

  void KalFit::makeTrack(TrkDef const& mytrk,TrkKalFit& myfit) {
// test if fitable
    if(fitable(mytrk)){
// create the hits. This also initializes T0
      makeHits(mytrk,myfit);
// Create the BaBar hit list.  This takes ownership
// Also create a straw hit intersection for each straw hit (active or not)
      std::vector<DetIntersection> detinter;
      TrkHotListFull* hotlist = new TrkHotListFull();
      for(std::vector<TrkStrawHit*>::iterator ihit=myfit._hits.begin();ihit!=myfit._hits.end();ihit++){
        TrkStrawHit* trkhit = *ihit;
        hotlist->append(trkhit);
        double fltlen = trkhit->fltLen();
// note this accounts for the material of both walls
        double wallpath = trkhit->wallPath();
        double gaspath = trkhit->gasPath();
  // offset the paths to avoid stacking elements
        double wlen = fltlen - 0.5*trkhit->straw().getRadius();
        double glen = fltlen + 0.5*trkhit->straw().getRadius();
        detinter.push_back(DetIntersection(&_wallelem,&mytrk.helix(),
          wlen,wlen-wallpath,wlen+wallpath));
        detinter.push_back(DetIntersection(&_gaselem,&mytrk.helix(),
          glen,glen-gaspath,glen+gaspath));
      }
// Create BaBar track and Kalman fit
      myfit._trk = new TrkRecoTrk(_kalcon->defaultType(), 0, 0);
      assert(myfit._trk != 0);
      myfit._trk->setBField(_bfield);
// create Kalman rep
      myfit._krep = new KalRep(mytrk.helix(), hotlist, detinter, myfit._trk, *_kalcon, PdtPid::electron);
      assert(myfit._krep != 0);
      myfit._trk->setRep(myfit._krep);
// fit the track
      myfit.fit();
// update t0, and propagate it to the hits
      double oldt0(-1e8);
      myfit._nt0iter = 0;
      while(_updatet0 && myfit._fit.success() && fabs(myfit._t0.t0()-oldt0) > _t0tol && 
      myfit._nt0iter < _kalcon->maxIterations()){
        oldt0 = myfit._t0.t0();
        if(updateT0(myfit)){
          myfit._krep->resetFit();
          myfit.fit();
          myfit._nt0iter++;
        } else
            break;
// drop outlyers
        if(_weedhits){
          myfit._nweediter = 0;
          weedHits(myfit);
        }
      }
    }
    if(_removefailed)myfit.removeFailed();
  }

  
  bool
  KalFit::fitable(TrkDef const& mytrk){
    return mytrk.strawHitIndices().size() >= _minnstraws;
  }
  
  void
  KalFit::makeHits(TrkDef const& mytrk, TrkKalFit& myfit) {
    const Tracker& tracker = getTrackerOrThrow();
    // find flightlength at z=0
    double flt0 = mytrk.helix().zFlight(0.0);
    unsigned nind = mytrk.strawHitIndices().size();
    double tsum(0.0);
    for(unsigned iind=0;iind<nind;iind++){
      unsigned istraw = mytrk.strawHitIndices()[iind];
      const StrawHit& strawhit(mytrk.strawHitCollection()->at(istraw));
      const Straw& straw = tracker.getStraw(strawhit.strawIndex());
    // compute initial flightlength from helix and hit Z
      double hflt = mytrk.helix().zFlight(straw.getMidPoint().z());
    // estimate the time the track reaches this hit, assuming speed-of-light travel along the helix. Should use actual
    // speed based on momentum and assumed particle species FIXME!!!
      double tprop = (hflt -flt0)/_vlight;
      double hitt0 = mytrk.trkT0().t0() + tprop;
    // subtract the propagation time and the average wire signal delay when computing hit time
    // using vlight = vwire, FIXME!!!!!
      tsum += strawhit.time() - tprop - straw.getHalfLength()/_vlight;
    // create the hit object
      TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,hitt0,fabs(mytrk.trkT0().t0Err()));
      assert(trkhit != 0);
    // refine the flightlength, as otherwise hits in the same plane are at exactly the same flt, which can cause problems
      TrkPoca poca(mytrk.helix(),hflt,*trkhit->hitTraj(),trkhit->hitLen());
      if(poca.status().success()){
        trkhit->setFltLen(poca.flt1());
        trkhit->setHitLen(poca.flt2());
      } else {
        trkhit->setFltLen(hflt);
        trkhit->setActivity(false);
        trkhit->setUsability(-2);
      }
      myfit._hits.push_back(trkhit);
    }
  // if the initial t0error was 0, compute t0 and override
    if(mytrk.trkT0().t0Err() < 0.0 && nind > 0){
  // assuming a flat drift time means correcting by half the maximum drift time
      double tmax = myfit._hits[0]->straw().getRadius()/_vdrift;
      double t0 = tsum/nind - 0.5*tmax;
  // estimate the error using the same assumption
      double t0err = tmax/sqrt(12*nind);
  // save this as the initial T0 value (and the current t0 value)
      myfit.setT00(TrkT0(t0,t0err));
      myfit.setT0(TrkT0(t0,t0err));
    } else {
  // initialize t0 directly from the tracking object
      myfit.setT00(mytrk.trkT0());
      myfit.setT0(mytrk.trkT0());
    }
  // update the hits
    for(std::vector<TrkStrawHit*>::iterator ihit= myfit._hits.begin();ihit != myfit._hits.end(); ihit++){
      double hitt0 = myfit._t0.t0() + ((*ihit)->fltLen() -flt0)/_vlight;
      (*ihit)->updateT0(hitt0,myfit._t0.t0Err());
    }
  // sort the hits by flightlength
    std::sort(myfit._hits.begin(),myfit._hits.end(),fltlencomp());
  }
  
  bool
  KalFit::updateT0(TrkKalFit& myfit){
    bool retval(false);
// need to have a valid fit
    if(myfit._krep->fitValid()){
// find the global fltlen associated with z=0.  This should be a piectraj function, FIXME!!!
      double loclen;
      double flt0 = 0.0;
      double dz(10.0);
      unsigned niter = 0;
      while(fabs(dz) > 1.0 && niter < _kalcon->maxIterations() ) {
        const HelixTraj* helix = dynamic_cast<const HelixTraj*>(myfit._krep->localTrajectory(flt0,loclen));
        flt0 += helix->zFlight(0.0)-loclen;
        dz = myfit._krep->traj().position(flt0).z();
        niter++;
      }
// find hits
      std::vector<double> hitst0; // store t0, to allow outlyer removal
      for(std::vector<TrkStrawHit*>::iterator ihit= myfit._hits.begin();ihit != myfit._hits.end(); ihit++){
        TrkStrawHit* hit = *ihit;
        if(hit->isActive() && hit->poca()!= 0 && hit->poca()->status().success()){
// copy the seed
          static TrkSimpTraj* straj = myfit._krep->seed()->clone();
// find the hit site in the rep
          const KalHit* hitsite = myfit._krep->findHotSite(hit);
// set helix to the local parameters EXCLUDING THIS HIT
          if(hitsite != 0 && myfit._krep->smoothedTraj(hitsite,straj)){
            TrkPoca poca(*straj,hit->fltLen(),*(hit->hitTraj()),hit->hitLen());
            if(poca.status().success()){
              double doca = fabs(poca.doca());
// require a minimum doca to avoid ambiguity bias.  mindoca shoudl be a parameter, FIXME!!!
              static double mindoca(0.4);
              if(doca > mindoca){
// propagation time to this hit from z=0.  This assumes beta=1, FIXME!!!
                double tflt = (hit->fltLen()-flt0)/_vlight;
// drift time of this hit (plus t0)
                double tdrift = hit->time() - tflt;
// t0 = Time difference between the drift time and the DOCA time.  sign of DOCA is irrelevant here.
                double hitt0 = tdrift - doca/_vdrift;
                hitst0.push_back(hitt0);
              }
            }
          }
        }
      }
      if(hitst0.size() >1){
// iterate over outlyer removal.  TrkT0 window should be a parameter, FIXME!!!!
        double nsig(2.5);
        bool changed(true);
        double t0(0.0);
        double t0err(-1.0);

        std::vector<bool> used(hitst0.size(),true);
        unsigned niter(0);
        while(changed && niter < 10){
          niter++;
          unsigned nactive(0);
          double t0sum(0.0);
          double t0sum2(0.0);
          for(unsigned ihit=0;ihit<hitst0.size();ihit++){
            if(used[ihit]){
              nactive++;
              t0sum += hitst0[ihit];
              t0sum2 += hitst0[ihit]*hitst0[ihit];
            }
          }
          t0 = t0sum/nactive;
          double t02 = t0sum2/nactive;
          double t0sig = sqrt(max(t02 - t0*t0,0.0));
// for now, a kludge factor for the t0 error
          t0err = 1.5*t0sig/sqrt(nactive);
          changed = false;
          for(unsigned ihit=0;ihit<hitst0.size();ihit++){
            bool useit = fabs(hitst0[ihit]-t0) < nsig*t0sig;
            changed |= useit != used[ihit];
            used[ihit] = useit;
          }
        }
// reset t0
        myfit._t0.setT0(t0,t0err);
// reset all the hit times
        for(std::vector<TrkStrawHit*>::iterator ihit= myfit._hits.begin();ihit != myfit._hits.end(); ihit++){
          TrkStrawHit* hit = *ihit;
// correct for flightlength.  Again assumes beta=1, FIXME!!!
          double hitt0 = t0 + (hit->fltLen()-flt0)/_vlight;
          hit->updateT0(hitt0,t0err);
        }
        retval = true;
      }
    }
    return retval;
  }

  bool
  KalFit::weedHits(TrkKalFit& myfit) {
    // Loop over HoTs and find HoT with largest contribution to chi2.  If this value
    // is greater than some cut value, deactivate that HoT and reFit
    bool retval(false);
    double worst = -1.;
    TrkHitOnTrk* worstHot = 0;
    TrkHotList* hots = myfit._krep->hotList();
    for (TrkHotList::nc_hot_iterator iHot = hots->begin(); iHot != hots->end(); ++iHot) {
      if (iHot->isActive()) {
        double resid, residErr;
        if(iHot->resid(resid, residErr, true)){
          double value = fabs(resid/residErr);
          if (value > _maxhitchi && value > worst) {
            worst = value;
            worstHot = iHot.get();
          }
        }
      }
    }
    if(0 != worstHot){
      retval = true;
      worstHot->setActivity(false);
      worstHot->setUsability(-5);
      myfit.fit();
      myfit._krep->addHistory(myfit._fit, "HitWeed");
      // Recursively iterate
      myfit._nweediter++;
      if (myfit._fit.success() && myfit._nweediter < _maxweed ) {
        retval |= weedHits(myfit);
      }
    }
    return retval;
  }
}

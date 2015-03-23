//
// Class to perform BaBar Kalman fit
//
// $Id: KalFitHack.cc,v 1.4 2014/04/08 04:25:46 murat Exp $
// $Author: murat $ 
// $Date: 2014/04/08 04:25:46 $
//

// framework
#include "fhiclcpp/ParameterSet.h"
// the following has to come before other BaBar includes
#include "BaBar/BaBar.hh"
#include "CalPatRec/inc/KalFitHack.hh"
#include "KalmanTests/inc/PanelAmbigResolver.hh"
#include "KalmanTests/inc/PocaAmbigResolver.hh"
#include "KalmanTests/inc/HitAmbigResolver.hh"
#include "CalPatRec/inc/HitAmbigResolverHack.hh"
#include "KalmanTests/inc/FixedAmbigResolver.hh"
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
    _herr(pset.get< vector<double> >("hiterr")),
    _maxdriftpull(pset.get<double>("maxDriftPull",10)),
    // t0 parameters
    _initt0(pset.get<bool>("initT0",true)),
    _updatet0(pset.get<bool>("updateT0",true)),
    fMinHitDrift(pset.get<double>("HitMinDrift")),
    fRdriftMinusDocaTol(pset.get<double>("RdriftMinusDocaTol")),
    _daveMode(pset.get<int>("daveMode" ,0)),
    _t0tol(pset.get< vector<double> >("t0Tolerance")),
    _t0errfac(pset.get<double>("t0ErrorFactor",1.2)),
    _mint0doca(pset.get<double>("minT0DOCA",-0.2)),
    _t0nsig(pset.get<double>("t0window",2.5)),
    _dtoffset(pset.get<double>("dtOffset")),
    fScaleErrDoublet(pset.get<double>("scaleErrDoublet")),
    fUseDoublets(0),
    fILoopUseDoublets(pset.get<int>("iLoopUseDoublets")),
    fMinDriftDoublet  (pset.get<double>("minDriftDoublet")),
    fDeltaDriftDoublet(pset.get<double>("deltaDriftDoublet")),
    fSigmaSlope       (pset.get<double>("sigmaSlope")),
    fMakeStrawHitModuleLabel(pset.get<std::string>("makeStrawHitModuleLabel")),
    //
    _removefailed(pset.get<bool>("RemoveFailedFits",true)),
    _minnstraws(pset.get<unsigned>("minnstraws",15)),
    _ambigstrategy(pset.get< vector<int> >("ambiguityStrategy")),
    _bfield(0),
    fNIter(0)
  {
    //    fStopwatch = new TStopwatch();
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
    if(_herr.size() <= 0) throw cet::exception("RECO")<<"mu2e::KalFitHack: no hit errors specified" << endl;
    if(_herr.size() != _ambigstrategy.size()) throw cet::exception("RECO")<<"mu2e::KalFitHack: inconsistent ambiguity resolution" << endl;
    if(_herr.size() != _t0tol.size()) throw cet::exception("RECO")<<"mu2e::KalFitHack: inconsistent ambiguity resolution" << endl;
    // construct the ambiguity resolvers
    for(size_t iambig=0;iambig<_ambigstrategy.size();++iambig){
      switch (_ambigstrategy[iambig] ){
	case fixedambig: default:
	  _ambigresolver.push_back(new FixedAmbigResolver(pset));
	  break;
	case hitambig:
	  _ambigresolver.push_back(new HitAmbigResolver(pset));
	  break;
	case panelambig:
	  _ambigresolver.push_back(new PanelAmbigResolver(pset));
	  break;
	case pocaambig:
	  _ambigresolver.push_back(new PocaAmbigResolver(pset));
	  break;
	case hitambig_hack:
	  //	  _ambigresolver.push_back(new PocaAmbigResolver(pset));
	  _ambigresolver.push_back(new HitAmbigResolverHack(pset,_herr[iambig]));
	  break;
      }
    }
  }

//-----------------------------------------------------------------------------
  KalFitHack::~KalFitHack(){
    for(size_t iambig=0;iambig<_ambigresolver.size();++iambig){
      delete _ambigresolver[iambig];
    }
    delete _bfield;

    //    delete fStopwatch;
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// first step: create list of doublets
//-----------------------------------------------------------------------------
  void KalFitHack::findDoublets (KalRep*                    krep, 
				 std::vector<TrkStrawHit*> *hits, 
				 DoubletCollection         *DCol){

    //first step: create list of doublets
   //    const TrkHotList* hot_list = krep->hotList();
    mu2e::TrkStrawHit* hit;

    int               nhits, station, panel;
    int               oldStation(-1), oldPanel(-1), idlast(0);
    int               trkshsize, shId, layer, istraw;
    
    CLHEP::Hep3Vector wdir, pos, posPanel, wpos[10], tmppos;    
    CLHEP::Hep3Vector tdir, trkpos;
    HepPoint          tpos;

    double            flen, ds, doca, rdrift, phiPanel;
    double            endTrk(0.0);//krep->endFoundRange();

    const Straw       *straw;

    if (_debug > 1){
      printf("[KalFitHack::findDoublets]-------------------------------------------------\n");
      printf("[KalFitHack::findDoublets]  i  shId  ch  panel  il   iw   driftR       doca\n");
      printf("[KalFitHack::findDoublets]-------------------------------------------------\n");
    }
    
    DCol->clear();
    nhits = hits->size();
    
    int multipletIndex(0);
    for (int i=0; i<nhits; ++i) {
//-----------------------------------------------------------------------------
// use active hits only
//-----------------------------------------------------------------------------
      hit       = hits->at(i);
      //      if (hit->isActive() == 0)                             continue;
      straw     = &hit ->straw();
      wdir      = straw->getDirection();
      pos       = straw->getMidPoint();
      station   = straw->id().getDevice();
      panel     = straw->id().getSector();
      shId      = straw->index().asInt();
//-----------------------------------------------------------------------------
// track info 
//-----------------------------------------------------------------------------
      HelixTraj trkHel(krep->helix(endTrk).params(),krep->helix(endTrk).covariance());
      flen   = trkHel.zFlight(pos.z());
      krep->traj().getInfo(flen, tpos, tdir);
//-----------------------------------------------------------------------------
// try to extrapolate helix a bit more accurately
//-----------------------------------------------------------------------------
      ds = (pos.z()-tpos.z())/tdir.z();
      krep->traj().getInfo(flen+ds, tpos, tdir);
      trkpos.set(tpos.x(), tpos.y(), tpos.z());

      //calculate distance of closest approach - from midwire to track
      HepPoint p1(pos.x(),pos.y(),pos.z());
      HepPoint p2(trkpos.x() ,trkpos.y() ,trkpos.z());
      
      TrkLineTraj trstraw(p1, wdir, 0., 0.);
      TrkLineTraj trtrk  (p2, tdir, 0., 0.);
//-----------------------------------------------------------------------------
// distance of closest approach calculated from the track to the hit,
// track trajectory is the first parameter
//-----------------------------------------------------------------------------
      TrkPoca poca  (trtrk,0.,trstraw,0.);
      doca   = poca.doca();
					// get the drift radius
      rdrift  =  hit->driftRadius();

      if (_debug > 1) {
	layer  = straw->id().getLayer();
	istraw = straw->id().getStraw();
	printf("[KalFitHack::findDoublets] %2i  %5i %3i  %4i  %3i %3i %8.3f %8.3f\n",
	       i, shId, station, panel, layer, istraw, rdrift, doca);
      }
//-----------------------------------------------------------------------------
// do not use straw hits with small drift radii
//-----------------------------------------------------------------------------
//      if (rdrift < fMinHitDrift)                           continue;

      if (station != oldStation) { 
//-----------------------------------------------------------------------------
// new chamber : create new doublet candidate
//-----------------------------------------------------------------------------
	DCol->push_back(Doublet(multipletIndex, station, panel, wdir, tdir, trkpos, hit));
	oldStation = station;
	oldPanel   = panel;
	++idlast;
	++multipletIndex;
      } 
      else {
	if (panel == oldPanel) {
//-----------------------------------------------------------------------------
// same chamber, same panel : add one more hit to the last doublet
//-----------------------------------------------------------------------------
	  DCol->at(idlast-1).addStrawHit(tdir, trkpos, hit);
	}
	else {
//-----------------------------------------------------------------------------
// same chamber, different panel : new doublet candidate
//-----------------------------------------------------------------------------
	  DCol->push_back(Doublet(multipletIndex, station, panel, wdir, tdir, trkpos, hit));
	  oldStation = station;
	  oldPanel   = panel;
	  ++idlast;
	  ++multipletIndex;
	}
      }
    }
//-----------------------------------------------------------------------------
// list of doublets is formed, the rest of this routine - debugging only
//-----------------------------------------------------------------------------
    Doublet             *doublet;
    CLHEP::HepRotationZ rot;

    if (_debug > 1){
      printf("[KalFitHack::findDoublets] BEGIN iherr:%i \n", fAnnealingStep);
      printf("-----------------------------------------------------------------------------------------------------------------------------------------\n");
      printf("  i  shId ch pnl lay str      x        y         z      sinphi  tksphi    xtrk     ytrk      ztrk      xr      yr     zr      doca   rdr \n");
      printf("-----------------------------------------------------------------------------------------------------------------------------------------\n");
    }

    int      ndoublets = DCol->size();
    for (int i=0; i<ndoublets; ++i){
      doublet   = &DCol->at(i);
      trkshsize = doublet->fNstrawHits;
//-----------------------------------------------------------------------------
// assume wires are perpendicular to the radial direction to the panel
// this is already ambiguos
// use atan2 to get the right quadrant
//-----------------------------------------------------------------------------
      posPanel = doublet->fHit[0]->straw().getMidPoint();
      phiPanel = atan2(posPanel.y(),posPanel.x());
      rot.set(-phiPanel);
      posPanel = rot*posPanel;
      
      for (int j=0; j<trkshsize; ++j) {
	hit   = doublet->fHit[j];
	straw = &hit->straw();
	shId  = straw->index().asInt();
//-----------------------------------------------------------------------------
// mid-wire position and the wire direction
//-----------------------------------------------------------------------------
	wpos[j] = straw->getMidPoint();
	wdir    = doublet->fShDir;
//-----------------------------------------------------------------------------
// track position and direction
//-----------------------------------------------------------------------------
	trkpos  = doublet->fTrkPos[j];
	tdir    = doublet->fTrkDir[j];

	HepPoint p1(wpos[j].x(),wpos[j].y(),wpos[j].z());
	HepPoint p2(trkpos.x() ,trkpos.y() ,trkpos.z());
	
	TrkLineTraj trstraw(p1, wdir, 0., 0.);
	TrkLineTraj trtrk  (p2, tdir, 0., 0.);
	TrkPoca     poca   (trstraw, 0.,trtrk   , 0.);
	doca   = poca.doca();
	rdrift = hit->driftRadius();
//-----------------------------------------------------------------------------
// rotate into a coordinate system with X axis pointing towards the panel and 
// Y axis pointing in the wire direction
// current channel numbering scheme allows for that
//-----------------------------------------------------------------------------
	wpos[j] = rot*wpos[j];
	trkpos  = rot*trkpos;
	tdir    = rot*tdir;

	if (_debug > 1) {
	  printf(" %2i %5i %2i %3i %3i %3i %9.3f %9.3f %9.3f  %6.3f  %6.3f %9.3f %8.3f %9.3f %8.3f %4.1f %9.3f %6.3f %5.3f\n",
		 i, shId, doublet->fStationId, doublet->fPanelId, 
		 straw->id().getLayer(),
		 straw->id().getStraw(),
		 straw->getMidPoint().x(), straw->getMidPoint().y(), straw->getMidPoint().z(),
		 wdir.y(), 
		 tdir.x()/tdir.z(),
		 trkpos.x(), trkpos.y(), trkpos.z(),
		 wpos[j].x(), wpos[j].y(), wpos[j].z(),
		 doca,
		 rdrift);
	}
      }
    }
 } 

//--------------------------------------------------------------------------------
// given a multiplet, resolve the ambiguity for hit: index0 and index1
//--------------------------------------------------------------------------------
  void KalFitHack::markDoublet(Doublet *doublet, int index0, int index1) {
    mu2e::TrkStrawHit *hit  [2];
    const mu2e::Straw *straw[2];
   
    CLHEP::Hep3Vector spos[2], sposr[2], sdir[2], sdirr[2], wpos[10], posPanel;    
    CLHEP::Hep3Vector tpos[2], tposr[2], tdir[2], tdirr[2];
    
    CLHEP::Hep3Vector wdir, wdir1, wdir2;

    int               layer[2], ibest, inext;
    double            rdrift[2], phiPanel; 
    
    int               shId[2];
    double            trkslope, lineSlopes[4], dxdz[2], chi2[4], doca[4][2];
    double            xdr, dsl, xdsl, sig, chi2min, chi2next;

    double            sflt[2], tflt[2];
    HepPoint          spi[2] , tpi[2], hpos[2];
    Hep3Vector        sdi[2] , tdi[2], u[2];
    TrkPoca           poca[2];
    HepRotationZ      rot;

    wdir = doublet->fShDir;

//-----------------------------------------------------------------------------
// create the array holding the indexes of the straw hit to use
// within a multiplet
//-----------------------------------------------------------------------------
    int               index[2] = {index0, index1};
//-----------------------------------------------------------------------------
// by construction, both hits are in the same panel
//-----------------------------------------------------------------------------
    for (int i=0; i<2; i++) {
      hit   [i] = doublet->fHit[index[i]];
      straw [i] = &hit[i]->straw();
      layer [i] = straw[i]->id().getLayer();
      rdrift[i] = hit[i]->driftRadius();
      shId  [i] = straw[i]->index().asInt();
    }
//-----------------------------------------------------------------------------
// skip doublets with both hits in the same layer
//-----------------------------------------------------------------------------
// 2015-03-22 P.Murat   if (layer[0] == layer[1])                    continue ;
      for (int i=0; i<2; i++) {
	spos [i] = straw[i]->getMidPoint();
	phiPanel = std::atan2(spos[i].y(),spos[i].x());
	rot.set(-phiPanel);

	sposr[i] = rot*spos[i];

	sdir [i] = straw[i]->getDirection();
	sdirr[i] = rot*sdir[i];

	tpos [i] = doublet->fTrkPos[index[i]];
	tposr[i] = rot*tpos[i];

	tdir [i] = doublet->fTrkDir[index[i]];
	tdirr[i] = rot*tdir[i];

	dxdz [i] = tdirr[i].x()/tdirr[i].z();
      }
//-----------------------------------------------------------------------------
// choose the best combination of the drift signs - the one corresponding 
// to the slope closest to that of the track
// also use the dist of closest approach information
//
// 1. coordinate term
//-----------------------------------------------------------------------------
      for (int is=0; is<4; is++) {
	chi2[is] = 0;
      }

      for (int ih=0; ih<2; ih++) {
	  
	HepPoint    p1(sposr[ih].x(),sposr[ih].y(),sposr[ih].z());
	HepPoint    p2(tposr[ih].x(),tposr[ih].y(),tposr[ih].z());

	TrkLineTraj st   (p1,sdirr[ih],0.,0.);
	TrkLineTraj tt   (p2,tdirr[ih],0.,0.);

	poca[ih] = TrkPoca(st,0.,tt,0.);

	sflt[ih] = poca[ih].flt1();
	tflt[ih] = poca[ih].flt2();

	st.getInfo(sflt[ih],spi[ih],sdi[ih]);
	tt.getInfo(tflt[ih],tpi[ih],tdi[ih]);

	u[ih]    = sdi[ih].cross(tdi[ih]).unit();  // direction towards the center

	for (int is=0; is<4; is++) {
//-----------------------------------------------------------------------------
// hit position, given a drift sign
//-----------------------------------------------------------------------------
	  hpos[ih]     = spi[ih]+u[ih]*rdrift[ih]*fSign[is][ih];
	  doca[is][ih] = (hpos[ih]-tpi[ih]).mag();
	  sig          = 2.5; // 1.; // hit[ih]->hitRms();
	  xdr          = doca[is][ih]/sig;
	  chi2[is]    += xdr*xdr;
	}
      }
//-----------------------------------------------------------------------------
// 2. add slope term to chi2
//    for the track dx/dz use average of the two dx/dz slopes 
//    calculated in the two layers corresponding to the doublet hits
//-----------------------------------------------------------------------------
      trkslope  = (dxdz[0]+dxdz[1])/2.;
      findLines(sposr,rdrift,lineSlopes);
    
      for (int is=0; is<4; is++) {
	dsl       = fabs(trkslope-lineSlopes[is]);
	xdsl      = dsl/fSigmaSlope;
	chi2[is] += xdsl*xdsl;
      }
//-----------------------------------------------------------------------------
// now find the best solution
//-----------------------------------------------------------------------------
      ibest    = -1;
      inext    = -1;
      chi2min  = 1.e12;
      chi2next = 1.e12;
    
      for (int is=0; is<4; is++) {
	if (chi2[is] < chi2min) {
	  inext    = ibest;
	  chi2next = chi2min;
	  ibest    = is;
	  chi2min  = chi2[is];
	}
	else if (chi2[is] < chi2next) {
	  inext    = is;
	  chi2next = chi2[is];
	}
      }
//-----------------------------------------------------------------------------
// set best solutions
//-----------------------------------------------------------------------------
      int    os      = fSign[ibest][0]+fSign[ibest][1];
      double ext_err = _herr[fAnnealingStep];

      doublet->fOs      = os;
      doublet->fIBest   = ibest;
      doublet->fINext   = inext;
      doublet->fTrkDxDz = trkslope;
      for (int is=0; is<4; is++) {
	doublet->fDxDz[is] = lineSlopes[is];
	doublet->fChi2[is] = chi2[is];
      }
    
      for (int i=0; i<2; i++) {
	hit[i]->setAmbigUpdate(false);
//-----------------------------------------------------------------------------
// update the straw hit info inside the doublet, however don't rush 
// to resolve the hit sign ambiguities, do it only when completely sure
// this code is executed after the standard HitAmbigResolver, so when not sure, 
// do nothing and default to HitAmbigResolver
//-----------------------------------------------------------------------------
	doublet->fStrawAmbig[index[i]] = fSign[ibest][i];
	if (os == 0) {
	  if (fabs(rdrift[0]+rdrift[1]) > 0.8) {
//-----------------------------------------------------------------------------
// OS doublet reliably resolved, reduce the error
//-----------------------------------------------------------------------------
	    if (rdrift[i] > fMinDriftDoublet) {
//-----------------------------------------------------------------------------
// the hit drift radius is large - reduce the external error
//-----------------------------------------------------------------------------
	      hit[i]->setExtErr(ext_err/fScaleErrDoublet);
	      hit[i]->setAmbig(fSign[ibest][i]);
	    }
	    else {
//-----------------------------------------------------------------------------
// small drift radius : keep the external error large and set the ambiguity 
// to zero to use the wire coordinate
//-----------------------------------------------------------------------------
	      hit[i]->setExtErr(2*rdrift[i]);
	      hit[i]->setAmbig(0);
	    }
	  }
	}
	else {
//-----------------------------------------------------------------------------
// SS doublet
//-----------------------------------------------------------------------------
	  if (fabs(rdrift[0]-rdrift[1]) < fDeltaDriftDoublet) {
//-----------------------------------------------------------------------------
// SS doublet with close radii, scale of uncertainty is defined by the radius
//-----------------------------------------------------------------------------
	    if (fAnnealingStep < fILoopUseDoublets) {
	      double err = fabs(rdrift[0]+rdrift[1])/2.;
	      hit[i]->setExtErr(err);
	      hit[i]->setAmbig(0);
	    }
	  }
	  else {
//-----------------------------------------------------------------------------
// SS doublet, the two radii are different 
//-----------------------------------------------------------------------------
	    if (chi2min < 50) {
//-----------------------------------------------------------------------------
// if the best chi2 is good, the doublet drift signs are determined reliably
//-----------------------------------------------------------------------------
	      if (rdrift[i] > fMinDriftDoublet) {
		hit[i]->setExtErr(ext_err/fScaleErrDoublet);
		hit[i]->setAmbig(fSign[ibest][i]);
	      }
	      else {
//-----------------------------------------------------------------------------
// small radius
//-----------------------------------------------------------------------------
		hit[i]->setExtErr(2*rdrift[i]);
		hit[i]->setAmbig(0);
	      }
	    }
	    else {
//-----------------------------------------------------------------------------
// the best chi2 is large - cant believe anything
//-----------------------------------------------------------------------------
	      if (fAnnealingStep < fILoopUseDoublets) {
		double err = fabs(rdrift[i]);
		hit[i]->setExtErr(err);
		hit[i]->setAmbig(0);
	      }
	    }
	  }
	}
      }
    
      if (_debug > 0) {
	for (int i=0; i<2; i++) {
	  printf(" %2i %5i %2i %3i %2i %2i %8.3f %8.3f %9.3f %6.3f",
		 i, shId[i], doublet->fStationId, doublet->fPanelId, 
		 straw[i]->id().getLayer(),
		 straw[i]->id().getStraw(),
		 spos[i].x(), spos[i].y(), spos[i].z(),
		 wdir.y()
		 );
	  printf(" %6.3f %8.3f %8.3f %9.3f %8.3f %9.3f %6.3f %6.3f",
		 trkslope,
		 tpos[i].x(),tpos[i].y(),tpos[i].z(),
		 sposr[i].x(),tposr[i].x(),doca[ibest][i],rdrift[i]
		 ); 
	  printf(" %2i %6.3f %6.3f %8.2e %6.3f %6.3f %8.2e %6.3f %6.3f %8.2e %6.3f %6.3f %8.2e\n",
		 fSign[ibest][i],
		 lineSlopes[0], doca[0][i], chi2[0], 
		 lineSlopes[1], doca[1][i], chi2[1], 
		 lineSlopes[2], doca[2][i], chi2[2], 
		 lineSlopes[3], doca[3][i], chi2[3]
		 );
	}
      }
      //    }
  }

//---------------------------------------------------------------------------
// loop over the doublets found and mark their ambiguities
//---------------------------------------------------------------------------

  void KalFitHack::markMultiplets (DoubletCollection *DCol) {

    mu2e::TrkStrawHit *hit;
    const mu2e::Straw *straw;
    Doublet           *doublet;
    
    int ndhits;
    int ndoublets  = DCol->size();

    if (_debug > 0) {
      printf("[KalFitHack::markDoublets] BEGIN iherr:%i , ILoopUseDoublets:%2i\n", 
	     fAnnealingStep,fILoopUseDoublets);
      printf("----------------------------------------------------");
      printf("------------------------------------------------------------------------------");
      printf("------------------------------------------------------------------------------------------\n");
      printf("  i  shId ch pnl il is      x       y        z      ");
      printf("  cth   trkth    xtrk     ytrk     ztrk       xr      xtrkr   doca    rdr am  ");
      printf(" s++  doca++   chi2++   s+-  doca+-   chi2+-   s--  doca--   chi2--   s-+  doca-+   chi2-+\n");
      printf("----------------------------------------------------");
      printf("------------------------------------------------------------------------------");
      printf("------------------------------------------------------------------------------------------\n");
    }    

    for (int i=0; i<ndoublets; ++i) {
      doublet = &DCol->at(i);
      ndhits  = doublet->fNstrawHits;//fTrkshcol.size();
      
      if (ndhits < 2) continue;               // goto next doublet
     
      if (ndhits == 2) {
	markDoublet(doublet, 0, 1);
      }
      
      if (ndhits > 2) {
	int      tmpLayerId, layer0, layer1;
	int      tmpId(-1), id0(-1), id1(-1);
	int      nDoublets(0);
	double   rdrift;

	for (int j=0; j<ndhits; ++j){
	  hit     = doublet->fHit[j];
	  straw   = &hit->straw();
	  layer0  = straw->id().getLayer();
	  id0     = straw->index().asInt();
	 
	  if (nDoublets  == 1)      goto NEXT_DOUBLET;
	  
	  for (int k=j+1; k<ndhits; ++k){
	    hit     = doublet->fHit[k];
	    straw   = &hit->straw();
	    layer1  = straw->id().getLayer();
	    id1     = straw->index().asInt();

	    // 2015-03-22 P.Murat	    if (layer1 == layer0)  continue;
	    if (nDoublets  == 1  ) continue;
	                                                       //try to search for a doublet
	    markDoublet(doublet, j, k);                        //request of both: doublet found and oppposite
	    if ( doublet->fStrawAmbig[j] * doublet->fStrawAmbig[k] < 0) {        // ambiguity signs for the two straw hits
	      nDoublets = 1;

	      //now adjust the ambiguity sign of the other strawhits
	      for (int h=0; h<ndhits; ++h){
		hit        = doublet->fHit[h];
		straw      = &hit->straw();
		tmpLayerId = straw->id().getLayer();
		tmpId      = straw->index().asInt();
		
		if ( (h == j) || (h == k)) continue;

		if (tmpLayerId == layer0){
		  if (tmpId != id0){
		    doublet->fStrawAmbig[h] = -doublet->fStrawAmbig[j];
		  }
		}else if (tmpLayerId == layer1){
		  if (tmpId != id1){
		    doublet->fStrawAmbig[h] = -doublet->fStrawAmbig[k];
		  }
		}
		
		rdrift = hit->driftRadius();

		if ( fabs(rdrift) < fDeltaDriftDoublet){
		//set the hit ambiguity
		  hit->setAmbig(doublet->fStrawAmbig[h]);
		  hit->setAmbigUpdate(false);
		}
   
	      }//end of the loop over the hits in the multiplet
	    }
	    
	  }

	NEXT_DOUBLET:;
	}
      }
    }
//-----------------------------------------------------------------------------
// drift signs for doublet hits are defined, set hit drift signs
//-----------------------------------------------------------------------------
    if (_debug > 0) {
      printf("----------------------------------------------------------\n");
      printf("[KalFitHack::markMultiplets]   shId    sec   panel  ambig \n");
      printf("----------------------------------------------------------\n");
    }
  }

//--------------------------------------------------------------------------------
// set hit drift sign
//-----------------------------------------------------------------------------
  void KalFitHack::findAndMarkMultiplets(KalRep* krep, std::vector<TrkStrawHit*> *hits) {

    //    DoubletCollection DCol;
    mu2e::TrkStrawHit *hit;
//-----------------------------------------------------------------------------
// set external errors
// reduce external errors for hits from doublets during first iterations
//-----------------------------------------------------------------------------
    int nhits = hits->size();
    for (int i=0; i<nhits; ++i){
      hit   = hits->at(i);
      hit->setExtErr(_herr[fAnnealingStep]);
    }

    //    if (fUseDoublets == 0) return;
    if ( (fAnnealingStep >= fILoopUseDoublets) || (fUseDoublets == 0) ) return;
//-----------------------------------------------------------------------------
// create list of doublets 
//-----------------------------------------------------------------------------
    findDoublets (krep, hits, &fListOfDoublets);

    if (_debug >0) {
      printf("[KalFitHack::findAndMarkMultiplets] found %lu multiplets\n",
	     fListOfDoublets.size());
    }
//-----------------------------------------------------------------------------
// resolve drift signs for hits in doublets. Choose the combination of drift 
// signs for which the 2-hit segment slope is the closest to that of the track 
//-----------------------------------------------------------------------------
    if (fListOfDoublets.size() > 0) markMultiplets(&fListOfDoublets);
  }

//-----------------------------------------------------------------------------
// 2015-02-25 P.Murat: new resolver
// assume that coordinates are rotated into the coordinate system where Y axis 
// is pointed along the wire by a rotation along the Z axis
// sign ordering convention:    ++, +-, --, -+ is defined by KalFitHack::fSign
//-----------------------------------------------------------------------------
  void KalFitHack::findLines(Hep3Vector* Pos, double* R, double* Slopes) {
//-----------------------------------------------------------------------------
    double            nx[4], ny[4], dx, dy; //, ux, uy;
    double            alpha, dr, dr2, lx, ly;
    //    int               invert[4];

    dx  = Pos[1].z()-Pos[0].z();
    dy  = Pos[1].x()-Pos[0].x();
    dr2 = dx*dx+dy*dy;
    dr  = sqrt(dr2);

    lx  = dx/dr;
    ly  = dy/dr;

    for (int i=0; i<4; i++) {
      alpha = (R[1]*fSign[i][1]-R[0]*fSign[i][0])/dr;
      nx[i] = fSign[i][0]*(-ly*alpha+lx*sqrt(1-alpha*alpha));
      ny[i] = fSign[i][0]*( lx*alpha+ly*sqrt(1-alpha*alpha));

//       invert[i] =  1;
//       if (nx[i]*dx+ny[i]*dy < 0) {
// 	   nx[i]     = -nx[i];
// 	   ny[i]     = -ny[i];
// 	   invert[i] = -1;
//       }
      
      Slopes[i] = ny[i]/nx[i];
    }
  }

//------------------------------------------------------------------------------------------
// called once per event from CalPatRec_module::produce
//-----------------------------------------------------------------------------
  void KalFitHack::makeTrack(KalFitResult& kres, CalTimePeak* TPeak, int USE_DOUBLETS) {

    kres._fit = TrkErrCode(TrkErrCode::fail);

					// test if fitable
    if (fitable(kres._tdef)) {
					// first, find t0
      TrkT0 t0;
      bool caloInitCond(false);
      if (TPeak->Cluster() != NULL) {
	caloInitCond = true;
      }

      if (_initt0) {
	if (!caloInitCond) {
	  initT0(kres._tdef, t0);
	} 
	else {
	  initCaloT0(TPeak, kres._tdef, t0);
	}
      }
      else {
	t0 = kres._tdef.t0();
      }
					// knowing t0, create the hits
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
      kres._krep = new KalRep(kres._tdef.helix(), hotlist, kres._detinter, *this, kres._tdef.particle());
      assert(kres._krep != 0);
// initialize krep t0; eventually, this should be in the constructor, FIXME!!!
      double flt0 = kres._tdef.helix().zFlight(0.0);
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

      fUseDoublets = USE_DOUBLETS;
      fListOfDoublets.clear();
//-----------------------------------------------------------------------------
// 09 - 26 - 2013 giani 
// include the calorimeter information when it is avaiable
//-----------------------------------------------------------------------------
      if ((_daveMode == 0) && caloInitCond) {
	fitTrack(kres, TPeak);
      } 
      else {
	fitTrack(kres);
      }

      if (_removefailed) kres.removeFailed();
    }
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void KalFitHack::addHits(KalFitResult&              kres   , 
			   const StrawHitCollection*  straws , 
			   std::vector<hitIndex>      indices, 
			   double                     maxchi ,
			   CalTimePeak*               tpeak   ) {

					// there must be a valid Kalman fit to add hits to
    int activity(0);
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

      for(unsigned iind=0;iind<indices.size(); ++iind){
	size_t istraw = indices[iind]._index;
	const StrawHit& strawhit(straws->at(istraw));
	const Straw& straw = tracker.getStraw(strawhit.strawIndex());
// estimate  initial flightlength
	double hflt(0);
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
	activity = 1;
// if it's outside limits, deactivate the HOT
	if (chi > maxchi || !trkhit->physicalDrift(maxchi)) {
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
      std::sort(kres._hits.begin(),kres._hits.end(),fltlencomp(kres._tdef.fitdir().fitDirection()));

// refit the track one more time with minimal external errors
//---------------------------------------------------------------------------
//2015 - 02 - 27 Gianipez added the loop for including the external errors 
//---------------------------------------------------------------------------
      if (tpeak){
//------------------------------------------------------------------------------------------
// 2015 - 03 - 09 Gainipez added the following line for forcing the fiITeration procedure
// to use findAndUseDoublets
//------------------------------------------------------------------------------------------
	for (size_t iherr=_herr.size()-2; iherr<_herr.size();++iherr) {
	  fitIteration(kres, iherr, tpeak);
	}
      }else{
	fitIteration(kres,_herr.size()-1);
      }
      kres._krep->addHistory(kres._fit,"AddHits");
    }
  }


  void KalFitHack::fitTrack(KalFitResult& kres, CalTimePeak* TPeak) {
// loop over external hit errors, ambiguity assignment, t0 tolerance
// 10-03-2013 giani changed this loop. now it loops on all the stations
// and store the last fit that converges

    for (size_t iherr=0; iherr<_herr.size();++iherr) {
      fitIteration(kres,iherr,TPeak);

      if (! kres._fit.success()) break; //commented by gianipez
     }

    if(kres._krep != 0) kres._krep->addHistory(kres._fit,"KalFitHack");
  }

//-----------------------------------------------------------------------------
// one step of the track fit
//-----------------------------------------------------------------------------
  void KalFitHack::fitIteration(KalFitResult& kres, size_t IHErr, CalTimePeak* TPeak) {
    // update external hit errors.  This isn't strictly necessary on the 1st iteration.

    //    TrkStrawHit* hit;
    //    int          nhits;
    double       oldt0 = kres._krep->t0()._t0;
    unsigned     niter(0);
    bool         changed(true);
    bool         fit_success;


    if (_debug >0) {
      printf("------------------------------------------------------------------------------------------\n");
      printf("[KalFitHack::fitIteration] BEGIN IHErr:%i \n", int(IHErr));
      printf("------------------------------------------------------------------------------------------\n");
    }

    fAnnealingStep = IHErr;

//-----------------------------------------------------------------------------------
// 2015 -02 -17 G. Pezzullo: loop over the hits and assign a smaller external error 
// for the doublets
//-----------------------------------------------------------------------------------
    fAnnealingStep = IHErr;
    //    if ((fAnnealingStep < fILoopUseDoublets) && (fUseDoublets==1)) {
//--------------------------------------------------------------------------------
// 2015-02-19 G. Pezzu: re-search multiplets using updated fit results
//-----------------------------------------------------------------------------
    findAndMarkMultiplets(kres._krep, &kres._hits);
    if (_debug>0){
      printHits(kres,"fitIteration_001");
    }
      //   }

    kres._nt0iter = 0;
    kres._fit     = TrkErrCode::succeed;

    while (kres._fit.success() && changed && niter < maxIterations()) {
      changed = false;
      _ambigresolver[IHErr]->resolveTrk(kres);
//--------------------------------------------------------------------------------
//2015-02-17 G. Pezzu: fix the ambiguity of the doublets!
//2015-02-19 G. Pezzu: re-search hit multiplets using updated fit results
//--------------------------------------------------------------------------------

	findAndMarkMultiplets(kres._krep, &kres._hits);
	if (_debug > 0) {
	  printHits(kres,"fitIteration_002");
	}

//--------------------------------------------------------------------------------
// perform the track fit
//-----------------------------------------------------------------------------
      kres._krep->resetFit();
      kres.fit();

      fit_success = kres._fit.success();
      if (! fit_success) break;
//-----------------------------------------------------------------------------
// if the fit succeeded, update track T0, and recalculate T0's of the individual 
// hits
//-----------------------------------------------------------------------------
      if (_updatet0 ) {
	// 2014-12-11: G.Pezzullo and P.Murat - temporary *FIXME*
	if (TPeak != NULL)  updateCalT0(kres,TPeak);
	else                updateT0(kres);

	changed |= fabs(kres._krep->t0()._t0-oldt0) > _t0tol[IHErr];
	oldt0    = kres._krep->t0()._t0;
      }
//-----------------------------------------------------------------------------
// drop outliers. weedHits() calls KalRep::fit(), so the fit success code may
// may be redefined there
//-----------------------------------------------------------------------------
      if(_weedhits){
	kres._nweediter = 0;
	changed        |= weedHits(kres);
	fit_success     = kres._fit.success();
      }
      niter++;
    }
//-----------------------------------------------------------------------------
// done iterating, define drift signs with respect to the final trajectory
// 2015-02-17 G. Pezzu: ::resolveTrk() updates drift signs of ALL hits, 
// so fix the ambiguity of the doublets after that
//-----------------------------------------------------------------------------
    fNIter += niter;

    kres._ninter = kres._krep->intersections();
  }


//-----------------------------------------------------------------------------
  bool KalFitHack::fitable(TrkDef const& tdef){
    return tdef.strawHitIndices().size() >= _minnstraws;
  }
  
//-----------------------------------------------------------------------------
  void KalFitHack::makeHits(KalFitResult& kres, TrkT0 const& t0) {
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

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void KalFitHack::printHits(KalFitResult& kres, const char* Caller) {
    const KalRep* Trk  = kres._krep;
    double d0(-1.), om(-1.), r(-1.), phi0(-1.), x0(-1.), y0(-1.), chi2N(-1.);

    if (Trk != 0) {
      d0    = Trk->helix(0.).d0();
      om    = Trk->helix(0.).omega();
      r     = fabs(1./om);
      phi0  = Trk->helix(0.).phi0();
      x0    =  -(1/om+d0)*sin(phi0);
      y0    =   (1/om+d0)*cos(phi0);
      chi2N = Trk->chisq()/(Trk->nDof());

      printf("[KalFitHack::printHits] BEGIN called from %s iherr:%i \n",Caller,fAnnealingStep);
      printf("---------------------------------------------------------------------------------");
      printf("-----------------------------------------------------\n");
      //      printf("%s",Prefix);
      printf("  TrkID       Address    N  NA      P       pT     costh    T0      T0Err   Omega");
      printf("      D0       Z0      Phi0   TanDip    Chi2    FCons\n");
      printf("---------------------------------------------------------------------------------");
      printf("-----------------------------------------------------\n");

      Hep3Vector trk_mom;
      //      Trk->printAll();
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
    int nhits = kres._hits.size();
    printf("--------------------------------------------------------------------");
    printf("----------------------------------------------------------------");
    printf("--------------------------------------------\n");
    printf(" ih  SInd U A     len         x        y        z      HitT    HitDt");
    printf(" Ch Pl  L  W     T0       Xs      Ys        Zs     resid sigres");
    printf(" Rdrift   mcdoca  totErr hitErr  t0Err penErr extErr\n");
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
      hit   =  kres._hits.at(it);
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

      double step_doca = -99.0;

      if (step) {
	const Hep3Vector* v1 = &straw->getMidPoint();
	HepPoint p1(v1->x(),v1->y(),v1->z());

	const Hep3Vector* v2 = &step->position();
	HepPoint    p2(v2->x(),v2->y(),v2->z());

	TrkLineTraj trstraw(p1,straw->getDirection()  ,0.,0.);
	TrkLineTraj trstep (p2,step->momentum().unit(),0.,0.);

	TrkPoca poca(trstep, 0., trstraw, 0.);
    
	step_doca = poca.doca();
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

      printf("%8.3f %8.3f %9.3f %7.3f %7.3f",
	     pos.x(),
	     pos.y(),
	     pos.z(),
	     res,
	     sigres
	     );

      if (hit->ambig() != 0) printf(" %6.3f",hit->ambig()*hit->driftRadius());
      else                   printf(" *%5.3f",hit->driftRadius());

      printf("  %7.3f  %6.3f %6.3f %6.3f %6.3f %6.3f\n",		 
	     step_doca, 
	     hit->totalErr(),
	     hit->hitErr(),
	     hit->t0Err(),
	     hit->penaltyErr(),
	     hit->extErr()
	     );
    }
  }
 

  void
  KalFitHack::makeMaterials(KalFitResult& kres) {
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
  KalFitHack::weedHits(KalFitResult& kres) {
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
  KalFitHack::unweedHits(KalFitResult& kres, double maxchi) {
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
  KalFitHack::bField() const {
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
  KalFitHack::trkVolume(trkDirection trkdir) const {
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
    t0._t0err = .5;
  }


  void
  KalFitHack::initT0(TrkDef const& tdef, TrkT0& t0) {
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
  void KalFitHack::updateCalT0(KalFitResult& kres, CalTimePeak* TPeak) {
//    2014-11-24 gianipez and Pasha removed time offset between caloriemter and tracker

    TrkT0 t0;
    double mom, vflt, path, t0flt, flt0(0.0);
    bool converged = TrkHelixUtils::findZFltlen(kres._krep->traj(),0.0,flt0);
    
    //get helix from kalrep
    HelixTraj trkHel(kres._krep->helix(flt0).params(),kres._krep->helix(flt0).covariance());
    
					// get flight distance of z=0
    t0flt = trkHel.zFlight(0.0);
    
    if (converged) {
//-----------------------------------------------------------------------------
// estimate the momentum at that point using the helix parameters.  
// This is assumed constant for this crude estimate
// compute the particle velocity
//-----------------------------------------------------------------------------
      mom  = TrkMomCalculator::vecMom(trkHel,bField(),t0flt).mag();
      vflt = kres._tdef.particle().beta(mom)*CLHEP::c_light;
//-----------------------------------------------------------------------------
// path length of the particle from the middle of the Tracker to the  calorimeter
// set dummy error value
//-----------------------------------------------------------------------------
      path      = TPeak->ClusterZ()/trkHel.sinDip();
      t0._t0    = TPeak->ClusterT0() + _dtoffset - path/vflt;
      t0._t0err = .5;
      
      kres._krep->setT0(t0,flt0);
      updateHitTimes(kres);
    }
  }
  
  bool KalFitHack::updateT0(KalFitResult& kres) {
    using namespace boost::accumulators;
    bool retval(false);
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    KalRep* krep = kres._krep;
// need to have a valid fit
    if(krep->fitValid()) {
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
  KalFitHack::updateHitTimes(KalFitResult& kres) {
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
    if (_debug > 1) {
      printf("[KalFitHack::updateHitTimes] moving forward\n");
      printHits(kres,"updateTimes_001");
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

    if (_debug > 1) {
      printf("[KalFitHack::updateHitTimes] moving backwards\n");
      printHits(kres,"updateTimes_002");
    }

  }

  void
  KalFitHack::findBoundingHits(std::vector<TrkStrawHit*>& hits,double flt0,
    std::vector<TrkStrawHit*>::reverse_iterator& ilow,
    std::vector<TrkStrawHit*>::iterator& ihigh) {
    ilow = hits.rbegin();
    ihigh = hits.begin();
    while(ilow != hits.rend() && (*ilow)->fltLen() > flt0 )++ilow;
    while(ihigh != hits.end() && (*ihigh)->fltLen() < flt0 )++ihigh;
  }

}

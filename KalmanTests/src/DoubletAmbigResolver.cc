///////////////////////////////////////////////////////////////////////////////
// class to resolve hit ambiguities one hit at a time, assuming a reasonable track
// fit as input
///////////////////////////////////////////////////////////////////////////////
#include "KalmanTests/inc/DoubletAmbigResolver.hh"
#include "KalmanTests/inc/KalFitResult.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalSite.hh"
#include "KalmanTrack/KalHit.hh"
#include "TrkBase/TrkPoca.hh"
#include <vector>
#include <algorithm>
#include <functional>

///using namespace CLHEP;

namespace mu2e {
  typedef std::vector<TrkStrawHit*>::iterator TSHI;


//-----------------------------------------------------------------------------
  DoubletAmbigResolver::DoubletAmbigResolver(fhicl::ParameterSet const& pset,
					     int                        Iherr) :
    AmbigResolver(pset),
    _debugLevel  (pset.get<int>   ("debugLevel"      ,0        )),
    _mindrift    (pset.get<double>("HitMinDrift"     ,0.2      )),
    _zeropenalty (pset.get<double>("ZeroDriftPenalty",0.2      )),
    _penalty     (pset.get<bool>  ("HitAmbigPenalty" ,false    )),
    _expnorm     (pset.get<double>("HitExpNorm"      ,0.03907  )),
    _lambda      (pset.get<double>("HitLambda"       ,0.1254   )),
    _offset      (pset.get<double>("HitOffset"       ,0.073    )),
    _slope       (pset.get<double>("HitSlope"        ,-0.002374)),

    _sigmaSlope       (pset.get<double>("sigmaSlope")),
    _maxDoubletChi2   (pset.get<double>("maxDoubletChi2",9.)),
    _scaleErrDoublet  (pset.get<double>("scaleErrDoublet")),
    _minDriftDoublet  (pset.get<double>("minDriftDoublet")),
    _deltaDriftDoublet(pset.get<double>("deltaDriftDoublet"))
  {
    _iherr        = Iherr;
    _decisionMode = 0;
//-----------------------------------------------------------------------------
// initialize sequence of drift signs
//-----------------------------------------------------------------------------
    double s[4][2] = { 1, 1, 1, -1, -1, -1, -1, 1} ;
    for (int i=0; i<4; i++) {
      for (int j=0; j<2; j++) {
	_sign[i][j] = s[i][j];
      }
    }
  }


//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  DoubletAmbigResolver::~DoubletAmbigResolver() {}

//-----------------------------------------------------------------------------
// first step: create list of doublets
//-----------------------------------------------------------------------------
  void DoubletAmbigResolver::findDoublets (KalFitResult& KRes) const {
    mu2e::TrkStrawHit *hit;
    const Straw       *straw;

    int               nhits, station, panel;
    int               oldStation(-1), oldPanel(-1), idlast(0);
    int               trkshsize, shId, layer, istraw;
    
    CLHEP::Hep3Vector wdir, pos, posPanel, wpos[10], tmppos;    
    CLHEP::Hep3Vector tdir, trkpos;
    HepPoint          tpos;

    std::vector<Doublet>* dcol;

    double            flen, ds, doca, rdrift, phiPanel;
    double            endTrk(0.0);//Krep->endFoundRange();

    if (_debugLevel > 1){
      printf("[KalFitHack::findDoublets]-------------------------------------------------\n");
      printf("[KalFitHack::findDoublets]  i  shId  ch  panel  il   iw   driftR       doca\n");
      printf("[KalFitHack::findDoublets]-------------------------------------------------\n");
    }
    
    dcol = &KRes._listOfDoublets;

    dcol->clear();
    nhits = KRes._hits.size();
    
    int multipletIndex(0);
    for (int i=0; i<nhits; ++i) {
//-----------------------------------------------------------------------------
// use active hits only
//-----------------------------------------------------------------------------
      hit       = KRes._hits.at(i);
      //      idoublet  = -1;
      //      if (hit->isActive() == 0)                             goto END_OF_LOOP;
      straw     = &hit ->straw();
      wdir      = straw->getDirection();
      pos       = straw->getMidPoint();
      station   = straw->id().getDevice();
      panel     = straw->id().getSector();
      shId      = straw->index().asInt();
//-----------------------------------------------------------------------------
// set an external error
//-----------------------------------------------------------------------------
      hit->setExtErr(AmbigResolver::_extErr);
//-----------------------------------------------------------------------------
// track info 
//-----------------------------------------------------------------------------
      HelixTraj trkHel(KRes._krep->helix(endTrk).params(),KRes._krep->helix(endTrk).covariance());
      flen   = trkHel.zFlight(pos.z());
      KRes._krep->traj().getInfo(flen, tpos, tdir);
//-----------------------------------------------------------------------------
// try to extrapolate helix a bit more accurately
//-----------------------------------------------------------------------------
      ds = (pos.z()-tpos.z())/tdir.z();
      KRes._krep->traj().getInfo(flen+ds, tpos, tdir);
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

      if (_debugLevel > 1) {
	layer  = straw->id().getLayer();
	istraw = straw->id().getStraw();
	printf("[KalFitHack::findDoublets] %2i  %5i %3i  %4i  %3i %3i %8.3f %8.3f\n",
	       i, shId, station, panel, layer, istraw, rdrift, doca);
      }
//-----------------------------------------------------------------------------
// do not use straw hits with small drift radii
//-----------------------------------------------------------------------------
//      if (rdrift < fMinHitDrift)                           goto END_OF_LOOP;

      if (station != oldStation) { 
//-----------------------------------------------------------------------------
// new chamber : create new doublet candidate
//-----------------------------------------------------------------------------
	dcol->push_back(Doublet(multipletIndex, station, panel, wdir, tdir, trkpos, hit));
	oldStation = station;
	oldPanel   = panel;
	//	idoublet   = idlast;
	++idlast;
	++multipletIndex;
      } 
      else {
	if (panel == oldPanel) {
//-----------------------------------------------------------------------------
// same chamber, same panel : add one more hit to the last doublet
//-----------------------------------------------------------------------------
	  dcol->at(idlast-1).addStrawHit(tdir, trkpos, hit);
	  //	  idoublet = idlast-1;
	}
	else {
//-----------------------------------------------------------------------------
// same chamber, different panel : new doublet candidate
//-----------------------------------------------------------------------------
	  dcol->push_back(Doublet(multipletIndex, station, panel, wdir, tdir, trkpos, hit)); // 
	  oldStation = station;
	  oldPanel   = panel;
	  //	  idoublet   = idlast;
	  ++idlast;
	  ++multipletIndex;
	}
      }
      // END_OF_LOOP:; 
//-----------------------------------------------------------------------------
// in the very end of the loop, for each hit store association with the doublet
//-----------------------------------------------------------------------------
    }

//-----------------------------------------------------------------------------
// list of doublets is formed, the rest of this routine - diagnostics only
//-----------------------------------------------------------------------------
    Doublet             *doublet;
    CLHEP::HepRotationZ rot;

    int      ndoublets = dcol->size();

    if (_debugLevel >0) {
      printf("[KalFitHack::findDoublets] iherr:%i: found %i multiplets\n",_iherr,ndoublets);
      printf("--------------------------------------------------------------");
      printf("------------------------------------------------------------------------\n");
      printf("  i  shId ch pnl lay str      x        y         z      sinphi");
      printf("tksphi    xtrk     ytrk      ztrk      xr      yr     zr      doca   rdr\n");
      printf("--------------------------------------------------------------");
      printf("------------------------------------------------------------------------\n");
    }

    for (int i=0; i<ndoublets; ++i){
      doublet   = &dcol->at(i);
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

	if (_debugLevel > 1) {
	  printf(" %2i %5i %2i %3i %3i %3i %9.3f %9.3f %9.3f  %6.3f  ",
		 i, shId, doublet->fStationId, doublet->fPanelId, 
		 straw->id().getLayer(),
		 straw->id().getStraw(),
		 straw->getMidPoint().x(), 
		 straw->getMidPoint().y(), 
		 straw->getMidPoint().z(),
		 wdir.y()
		 );

	  printf("%6.3f %9.3f %8.3f %9.3f %8.3f %4.1f %9.3f %6.3f %5.3f\n",
		 tdir.x()/tdir.z(),
		 trkpos.x(), trkpos.y(), trkpos.z(),
		 wpos[j].x(), wpos[j].y(), wpos[j].z(),
		 doca,
		 rdrift
		 );
	}
      }
    }
  } 

//---------------------------------------------------------------------------
// resolve drift ambiguity for a single (non-doublet)  hit
//---------------------------------------------------------------------------
  void DoubletAmbigResolver::resolveSingleHit(KalFitResult& Kres, mu2e::TrkStrawHit* Hit) const {

    double                     doca[2];
    std::vector<TrkStrawHit*>  hits;

    const Straw* straw = &Hit->straw();
    
    const CLHEP::Hep3Vector& wdir = straw->getDirection();
    const CLHEP::Hep3Vector& wmid = straw->getMidPoint();
//-----------------------------------------------------------------------------
// calculate residuals for two hit positions corresponding to two different 
// drift signs
//-----------------------------------------------------------------------------
    hits.push_back(Hit);
    const TrkDifTraj* traj = findTraj(hits,Kres._krep);

    double dmin = Hit->timeDiffDist()-Hit->timeDiffDistErr();
    double dmax = Hit->timeDiffDist()+Hit->timeDiffDistErr();
    TrkLineTraj wtraj(HepPoint(wmid.x(),wmid.y(),wmid.z()),wdir,dmin,dmax);

    TrkPoca poca(*traj,Hit->fltLen(),wtraj,Hit->hitLen());
    if (poca.status().success()) {
//-----------------------------------------------------------------------------
// doca(hit) = doca(wire)-_iamb*radius
//-----------------------------------------------------------------------------
      doca[0] = poca.doca()-Hit->driftRadius();
      doca[1] = poca.doca()+Hit->driftRadius();

      double err = AmbigResolver::_extErr;
      double x0  = sqrt(doca[0]*doca[0]+err*err);
      double x1  = sqrt(doca[1]*doca[1]+err*err);

      int    ih;
      if (x0 < x1) ih = 0;
      else         ih = 1;
      
      if (fabs(doca[ih]/Hit->hitErr()) < 5.) {
//-----------------------------------------------------------------------------
// hit is close enough to the trajectory
//-----------------------------------------------------------------------------
	if (ih == 0) Hit->setAmbig( 1);
	else         Hit->setAmbig(-1);
      }
      else {
//-----------------------------------------------------------------------------
// can't tell
//-----------------------------------------------------------------------------
	Hit->setExtErr(Hit->driftRadius());
	Hit->setAmbig(0);
      }
    }
    else {
//-----------------------------------------------------------------------------
// couldn't determine doca
//-----------------------------------------------------------------------------
      Hit->setExtErr(Hit->driftRadius());
      Hit->setAmbig(0);
    }

    Hit->setAmbigUpdate(false);
  }

//-----------------------------------------------------------------------------
// 2015-02-25 P.Murat: new resolver
// assume that coordinates are rotated into the coordinate system where Y axis 
// is pointed along the wire by a rotation along the Z axis
// sign ordering convention:    ++, +-, --, -+ is defined by KalFitHack::_sign
//-----------------------------------------------------------------------------
  void DoubletAmbigResolver::findLines(Hep3Vector* Pos, double* R, double* Slopes) const {
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
      alpha = (R[1]*_sign[i][1]-R[0]*_sign[i][0])/dr;
      nx[i] = _sign[i][0]*(-ly*alpha+lx*sqrt(1-alpha*alpha));
      ny[i] = _sign[i][0]*( lx*alpha+ly*sqrt(1-alpha*alpha));

//       invert[i] =  1;
//       if (nx[i]*dx+ny[i]*dy < 0) {
// 	   nx[i]     = -nx[i];
// 	   ny[i]     = -ny[i];
// 	   invert[i] = -1;
//       }
      
      Slopes[i] = ny[i]/nx[i];
    }
  }



//--------------------------------------------------------------------------------
// given a multiplet, resolve the ambiguity for hit: index0 and index1
//--------------------------------------------------------------------------------
  void DoubletAmbigResolver::markDoublet(KalFitResult& KRes, Doublet *doublet, int index0, int index1) const {
    mu2e::TrkStrawHit *hit  [2];
    const mu2e::Straw *straw[2];
   
    CLHEP::Hep3Vector spos[2], sposr[2], sdir[2], sdirr[2], wpos[10], posPanel;    
    CLHEP::Hep3Vector tpos[2], tposr[2], tdir[2], tdirr[2];
    
    CLHEP::Hep3Vector wdir, wdir1, wdir2;

    int               /*layer[2],*/ ibest, inext;
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
      //      layer [i] = straw[i]->id().getLayer();
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
	hpos[ih]     = spi[ih]+u[ih]*rdrift[ih]*_sign[is][ih];
	doca[is][ih] = (hpos[ih]-tpi[ih]).mag();
	sig          = sqrt(rdrift[ih]*rdrift[ih] +0.1*0.1); // 2.5; // 1.; // hit[ih]->hitRms();
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
      xdsl      = dsl/_sigmaSlope;
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
    int    os         = _sign[ibest][0]+_sign[ibest][1];
    
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
      doublet->fStrawAmbig[index[i]] = _sign[ibest][i];
      if (os == 0) {
	if (fabs(rdrift[0]+rdrift[1]) > 0.8) {
//-----------------------------------------------------------------------------
// OS doublet reliably resolved, reduce the error
//-----------------------------------------------------------------------------
	  if (rdrift[i] > _minDriftDoublet) {
//-----------------------------------------------------------------------------
// the hit drift radius is large - reduce the external error
//-----------------------------------------------------------------------------
	    hit[i]->setExtErr(AmbigResolver::_extErr/_scaleErrDoublet);
	    hit[i]->setAmbig(_sign[ibest][i]);
	  }
	  else {
//-----------------------------------------------------------------------------
// small drift radius : unless forced, keep the external error large and set 
// the ambiguity to zero to use the wire coordinate
//-----------------------------------------------------------------------------
	    if (KRes._decisionMode == 0) {
	      hit[i]->setExtErr(rdrift[i]);
	      hit[i]->setAmbig(0);
	    }
	    else {
	      hit[i]->setAmbig(_sign[ibest][i]);
	    }
	  }
	}
      }
      else {
//-----------------------------------------------------------------------------
// SS doublet
//-----------------------------------------------------------------------------
	if ((KRes._decisionMode == 0) && (fabs(rdrift[0]-rdrift[1]) < _deltaDriftDoublet)) {
	  if ((chi2min < _maxDoubletChi2) && (chi2min/chi2next < 0.1)) {
//-----------------------------------------------------------------------------
// the best chi2 is good enough to rely on it
//-----------------------------------------------------------------------------
	    hit[i]->setExtErr(AmbigResolver::_extErr/_scaleErrDoublet);
	    hit[i]->setAmbig (_sign[ibest][i]);
	  }
	  else {
//-----------------------------------------------------------------------------
// SS doublet with close radii, scale of uncertainty is defined by the radius
//-----------------------------------------------------------------------------
	    double err = fabs(rdrift[i]);
	    hit[i]->setExtErr(err);
	    hit[i]->setAmbig(0);
	  }
	}
	else {
//-----------------------------------------------------------------------------
// SS doublet, the two radii are different or we're forced to make a decision
//-----------------------------------------------------------------------------
	  if (chi2min < _maxDoubletChi2) {
//-----------------------------------------------------------------------------
// the best chi2 is good, the doublet drift signs are determined reliably
//-----------------------------------------------------------------------------
	    if (rdrift[i] > _minDriftDoublet) {
	      hit[i]->setExtErr(AmbigResolver::_extErr/_scaleErrDoublet);
	      hit[i]->setAmbig (_sign[ibest][i]);
	    }
	    else {
//-----------------------------------------------------------------------------
// small radius (rdrift[i] < _minDriftDoublet) - use the wire position
//-----------------------------------------------------------------------------
	      hit[i]->setExtErr(rdrift[i]);
	      hit[i]->setAmbig(0);
	    }
	  }
	  else {
//-----------------------------------------------------------------------------
// the best double chi2 is large - cant believe anything, need to treat hits as 
// separate ones - this is to be implemented yet
// a good example - one of the hits - on Dave's no-gaussial tail
//-----------------------------------------------------------------------------
	    hit[i]->setExtErr(rdrift[i]);
	    hit[i]->setAmbig(0);
	    //	    resolveSingleHit(KRes,hit[i]);
	  }
	}
      }
    }
    
    if (_debugLevel > 0) {
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
	       _sign[ibest][i],
	       lineSlopes[0], doca[0][i], chi2[0], 
	       lineSlopes[1], doca[1][i], chi2[1], 
	       lineSlopes[2], doca[2][i], chi2[2], 
	       lineSlopes[3], doca[3][i], chi2[3]
	       );
      }
    }
  }


//---------------------------------------------------------------------------
// loop over the doublets found and mark their ambiguities
//---------------------------------------------------------------------------
  void DoubletAmbigResolver::markMultiplets (KalFitResult& Kres) const {

    mu2e::TrkStrawHit    *hit;
    const mu2e::Straw    *straw;
    Doublet              *doublet;
    std::vector<Doublet> *dcol;
    int                  ndoublets, ndhits;

    dcol       = &Kres._listOfDoublets;
    ndoublets  = dcol->size();

    if (_debugLevel > 0) {
      printf("[KalFitHack::markMultiplets] BEGIN iherr:%i\n",_iherr);
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
      doublet = &dcol->at(i);
      ndhits  = doublet->fNstrawHits;
      
      if (ndhits == 1) {
//-----------------------------------------------------------------------------
// a single hit in the plane - keep its external error large, unless the drift sign 
// can be determined reliably 
//-----------------------------------------------------------------------------
	hit = doublet->fHit[0];

	resolveSingleHit(Kres,hit);
      }
      else if (ndhits == 2) {
//-----------------------------------------------------------------------------
// 2 hits in a panel - attempt to determine the drift signs
//-----------------------------------------------------------------------------
	markDoublet(Kres,doublet,0,1);
      }
      else {
//-----------------------------------------------------------------------------
// more than 2 hits in a panel
//-----------------------------------------------------------------------------
	int      tmpLayerId, layer0, layer1, jbest(-1), kbest(-1);
	int      tmpId(-1), id0(-1), id1(-1);
	double   rdrift, chi2_d, chi2_best (1.e12);
	Doublet  bd;

	for (int j=0; j<ndhits; ++j) {
	  for (int k=j+1; k<ndhits; ++k){
	    // 2015-03-22 P.Murat	    if (layer1 == layer0)  continue;
//-----------------------------------------------------------------------------
// P.Murat: logic here looks questionalble, but let's first figure what it is exactly 
// - use the first found OS doublet, w/o looking at the quality
// - after it is found, resolve drift signs right away for all hits in the multiplet
// - this could be dangerous, especially, in presence of the background 
//-----------------------------------------------------------------------------
	    markDoublet(Kres,doublet,j,k); 
	    chi2_d = doublet->Chi2Best();
	    if (chi2_d < chi2_best) {
	      jbest     = j;
	      kbest     = k;
	      chi2_best = chi2_d;
	      bd        = *doublet;
	    }
	  }
	}

	*doublet = bd;
	if (chi2_best < _maxDoubletChi2) {
//-----------------------------------------------------------------------------
// the "best" doublet is good enough, resolve drift signs for the rest hits
//----------------------------------------------------------------------------- 
	  hit     = doublet->fHit[jbest];
	  straw   = &hit->straw();
	  layer0  = straw->id().getLayer();
	  id0     = straw->index().asInt();

	  hit     = doublet->fHit[kbest];
	  straw   = &hit->straw();
	  layer1  = straw->id().getLayer();
	  id1     = straw->index().asInt();

	  for (int h=0; h<ndhits; ++h) {
	    hit        = doublet->fHit[h];
	    straw      = &hit->straw();
	    tmpLayerId = straw->id().getLayer();
	    tmpId      = straw->index().asInt();
		
	    if ((h == jbest) || (h == kbest)) continue;
//-----------------------------------------------------------------------------
// the assumption here is that in case of a triplet two hits in the same layer 
// can't have the same drift sign - which is not necessarily correct.
// *stick to it for the time being*
//-----------------------------------------------------------------------------
	    if (tmpLayerId == layer0) {
	      if (tmpId != id0) {
		doublet->fStrawAmbig[h] = -doublet->fStrawAmbig[jbest];
	      }
	    }
	    else if (tmpLayerId == layer1) {
	      if (tmpId != id1) {
		doublet->fStrawAmbig[h] = -doublet->fStrawAmbig[kbest];
	      }
	    }
		
	    rdrift = hit->driftRadius();

	    if ( fabs(rdrift) < _deltaDriftDoublet){
	      hit->setAmbig(doublet->fStrawAmbig[h]);
	      hit->setAmbigUpdate(false);
	    }
	  }   
	}
      }
    }
  }

//-----------------------------------------------------------------------------
// build list of doublets and resolve the drift signs
//-----------------------------------------------------------------------------
  void DoubletAmbigResolver::resolveTrk(KalFitResult& KRes) const {

    findDoublets (KRes);

    if (KRes._listOfDoublets.size() > 0) markMultiplets(KRes);

  }

//-----------------------------------------------------------------------------
  double DoubletAmbigResolver::penaltyError(double rdrift) const {
    static double sqrt2 = sqrt(2.0);

    // model is of an exponential plus a linear.
    double frac = _expnorm*exp(-rdrift/_lambda)/_lambda + _offset + _slope*rdrift;
//-----------------------------------------------------------------------------
// the penalty term for a discrete ambiguity error depends only on the 
// mis-assignment probability and the drift distance
//-----------------------------------------------------------------------------
    double perr         = sqrt2*rdrift*sqrt(frac);
    return perr;
  }
}

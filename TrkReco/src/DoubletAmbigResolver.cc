///////////////////////////////////////////////////////////////////////////////
// class to resolve hit ambiguities one hit at a time, assuming a reasonable track
// fit as input
///////////////////////////////////////////////////////////////////////////////
#include "TrkReco/inc/DoubletAmbigResolver.hh"
#include "TrkReco/inc/TrkStrawHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalSite.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include <vector>
#include <algorithm>
#include <functional>

using CLHEP::Hep3Vector;
using CLHEP::HepRotationZ;

namespace mu2e {
  typedef std::vector<TrkStrawHit*>::iterator TSHI;


//-----------------------------------------------------------------------------
  DoubletAmbigResolver::DoubletAmbigResolver(fhicl::ParameterSet const& PSet,
					     double			extErr,
					     int                        Iter,
					     int			Final) :
    AmbigResolver(extErr),
    _debugLevel  (PSet.get<int>   ("debugLevel"      ,0        )),
    _mindrift    (PSet.get<double>("HitMinDrift"     ,0.2      )),
    _zeropenalty (PSet.get<double>("ZeroDriftPenalty",0.2      )),
    _penalty     (PSet.get<bool>  ("HitAmbigPenalty" ,false    )),
    _expnorm     (PSet.get<double>("HitExpNorm"      ,0.03907  )),
    _lambda      (PSet.get<double>("HitLambda"       ,0.1254   )),
    _offset      (PSet.get<double>("HitOffset"       ,0.073    )),
    _slope       (PSet.get<double>("HitSlope"        ,-0.002374)),
//-----------------------------------------------------------------------------
// parameters below are used for decisions on hit drift direction assignment, 
// but not for doublet reconstruction. 
// as doublet reconstruction is used in several places, need to be able to 
// create the ambig resolver w/o specifying any numbers in the talk to's
//-----------------------------------------------------------------------------
    _sigmaSlope       (PSet.get<double>("sigmaSlope"       ,  0.025)),
    _maxDoubletChi2   (PSet.get<double>("maxDoubletChi2"   , 10.   )),
    _scaleErrDoublet  (PSet.get<double>("scaleErrDoublet"  ,  5.   )),
    _minDriftDoublet  (PSet.get<double>("minDriftDoublet"  ,  0.3  )),
    _deltaDriftDoublet(PSet.get<double>("deltaDriftDoublet",  0.3  )),
    _excludeBothHits  (PSet.get<int>   ("excludeBothHits"  ,  1    )),
    _iter(Iter),
    _Final(Final)
  {
//-----------------------------------------------------------------------------
// initialize sequence of drift signs: (1,1) (1,-1) (-1,-1) (-1,1)
// 0 and 2: SS doublets, 1 and 3: OS doublets
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
  void DoubletAmbigResolver::findDoublets (const KalRep* KRep, vector<Doublet>* DCol) const {
    mu2e::TrkStrawHit *hit;

    int               station, panel;
    int               oldStation(-1), oldPanel(-1), idlast(0);
    int               trkshsize, shId, layer, istraw;
    
    Hep3Vector  wdir;
    Hep3Vector  pos, posPanel, wpos[10], tmppos;    
    Hep3Vector  tdir, trkpos;
    HepPoint           tpos;

    std::vector<Doublet>* dcol;

    double            flen, ds, doca, rdrift, phiPanel;
    double            endTrk(0.0);//Krep->endFoundRange();

    if (_debugLevel > 1){
      printf("[KalFitHack::findDoublets]-------------------------------------------------\n");
      printf("[KalFitHack::findDoublets]  i  shId  ch  panel  il   iw   driftR       doca\n");
      printf("[KalFitHack::findDoublets]-------------------------------------------------\n");
    }

    dcol = DCol;
    dcol->clear();

    //    nhits = KRes._hits.size();
    int multipletIndex(0);
//-----------------------------------------------------------------------------
// use active hits only
//-----------------------------------------------------------------------------
    //    for (int i=0; i<nhits; ++i) {
    int i = 0;
    TrkStrawHitVector tshv;
    convert(KRep->hitVector(),tshv);
    for (auto ihit=tshv.begin(); ihit!=tshv.end(); ++ihit) {
      //      idoublet  = -1;
      //      if (hit->isActive() == 0)                             goto END_OF_LOOP;
      Straw const& straw     = (*ihit) ->straw();
      wdir      = straw.getDirection();
      pos       = straw.getMidPoint();
      station   = straw.id().getPlane();
      panel     = straw.id().getPanel();
      shId      = straw.index().asInt();
//-----------------------------------------------------------------------------
// track info 
//-----------------------------------------------------------------------------
      HelixTraj trkHel(KRep->helix(endTrk).params(),KRep->helix(endTrk).covariance());
      flen   = trkHel.zFlight(pos.z());
      KRep->traj().getInfo(flen, tpos, tdir);
//-----------------------------------------------------------------------------
// try to extrapolate helix a bit more accurately
//-----------------------------------------------------------------------------
      ds = (pos.z()-tpos.z())/tdir.z();
      KRep->traj().getInfo(flen+ds, tpos, tdir);
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
      rdrift  =  (*ihit)->driftRadius();

      if (_debugLevel > 1) {
	layer  = straw.id().getLayer();
	istraw = straw.id().getStraw();
	printf("[KalFitHack::findDoublets] %2i  %5i %3i  %4i  %3i %3i %8.3f %8.3f\n",
	       i, shId, station, panel, layer, istraw, rdrift, doca);
      }
//-----------------------------------------------------------------------------
// do not use straw hits with small drift radii
//-----------------------------------------------------------------------------
      if (station != oldStation) { 
//-----------------------------------------------------------------------------
// new chamber : create new doublet candidate
//-----------------------------------------------------------------------------
	dcol->push_back(Doublet(multipletIndex, station, panel, wdir, tdir, trkpos, *ihit));
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
	  dcol->at(idlast-1).addStrawHit(tdir, trkpos, *ihit);
	  //	  idoublet = idlast-1;
	}
	else {
//-----------------------------------------------------------------------------
// same chamber, different panel : new doublet candidate
//-----------------------------------------------------------------------------
	  dcol->push_back(Doublet(multipletIndex, station, panel, wdir, tdir, trkpos, *ihit)); // 
	  oldStation = station;
	  oldPanel   = panel;
	  //	  idoublet   = idlast;
	  ++idlast;
	  ++multipletIndex;
	}
      }
      i += 1;
    }

//-----------------------------------------------------------------------------
// list of doublets is formed, the rest of this routine - diagnostics only
//-----------------------------------------------------------------------------
    Doublet             *doublet;
    CLHEP::HepRotationZ rot;

    int      ndoublets = dcol->size();

    if (_debugLevel >0) {
      printf("[KalFitHack::findDoublets] iherr:%i: found %i multiplets\n",_iter,ndoublets);
      printf("--------------------------------------------------------------");
      printf("------------------------------------------------------------------------\n");
      printf("  i  shId ch pnl lay str      x        y         z      sinphi");
      printf("tksphi    xtrk     ytrk      ztrk      xr      yr     zr      doca   rdr\n");
      printf("--------------------------------------------------------------");
      printf("------------------------------------------------------------------------\n");
    }

    for (int i=0; i<ndoublets; ++i){
      doublet   = &dcol->at(i);
      trkshsize = doublet->fNStrawHits;
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
	Straw const& straw = hit->straw();
	shId  = straw.index().asInt();
//-----------------------------------------------------------------------------
// mid-wire position and the wire direction
//-----------------------------------------------------------------------------
	wpos[j] = straw.getMidPoint();
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
		 straw.id().getLayer(),
		 straw.id().getStraw(),
		 straw.getMidPoint().x(), 
		 straw.getMidPoint().y(), 
		 straw.getMidPoint().z(),
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
  void DoubletAmbigResolver::resolveSingleHit(KalRep*      Krep, 
					      mu2e::TrkStrawHit* Hit ) const {

    double                     doca[2], xbest, xnext;
    std::vector<TrkStrawHit*>  hits;

    Straw const& straw = Hit->straw();
    
    const Hep3Vector& wdir = straw.getDirection();
    const Hep3Vector& wmid = straw.getMidPoint();
//-----------------------------------------------------------------------------
// calculate residuals for two hit positions corresponding to two different 
// drift signs
//-----------------------------------------------------------------------------
    hits.push_back(Hit);
    const TrkDifTraj* traj = findTraj(hits,Krep);

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

      int    ibest, inext;
      if (x0 < x1) {
	ibest = 0;
	inext = 1;
      }
      else {
	ibest = 1;
	inext = 0;
      }
      
      xbest = doca[ibest]/Hit->hitErr();
      xnext = doca[inext]/Hit->hitErr();
//-----------------------------------------------------------------------------
// want the best solution to be consistent with the trajectory, another one - 
// to be inconsistent, and the two - significantly different
//-----------------------------------------------------------------------------
      if ((fabs(xbest) < 5.) && (fabs(xnext) > 5) && (fabs(xbest/xnext) < 0.2)) {
//-----------------------------------------------------------------------------
// hit is close enough to the trajectory
//-----------------------------------------------------------------------------
	if (ibest == 0) Hit->setAmbig( 1);
	else            Hit->setAmbig(-1);
      }
      else {
//-----------------------------------------------------------------------------
// can't tell
//-----------------------------------------------------------------------------
	if (_Final == 0) {
	  Hit->setExtErr(Hit->driftRadius());
	  Hit->setAmbig(0);
	}
	else {
//-----------------------------------------------------------------------------
// no point to keep hits with very large residuals 
//-----------------------------------------------------------------------------
	  if (ibest == 0) Hit->setAmbig( 1);
	  else            Hit->setAmbig(-1);
	}
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


//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  int DoubletAmbigResolver::calculateDoubletParameters(const KalRep* KRep      , 
						       Doublet*      HitDoublet, 
						       Data_t*       R         ) const {
    
    mu2e::TrkStrawHit *hit  [2];
    //    const mu2e::Straw *straw[2];
   
    Hep3Vector /*spos[2], sposr[2], */ sdir[2], sdirr[2], wpos[10], posPanel;    
    Hep3Vector /*tpos[2], tposr[2], */ tdir[2], tdirr[2];
    
    //    const Hep3Vector wdir;
    Hep3Vector wdir1, wdir2;

    //    int               /*layer[2],*/ ibest, inext;
    double            phiPanel;
    
    //    int               shId[2];
    double            trkslope, /*lineSlopes[4],*/ dxdz[2], /*chi2[4], doca[4][2],*/ xdr[4][2];
    double            dsl, xdsl, sig;

    double            sflt[2], tflt[2];
    HepPoint          spi[2] , tpi[2], hpos[4][2];
    Hep3Vector        sdi[2] , tdi[2], u[2];
    TrkPoca           poca[2];
    HepRotationZ      rot;

    //    int               index[2] = {Index0, Index1};

    //    wdir  = &HitDoublet->fShDir;
//-----------------------------------------------------------------------------
// by construction, both hits are in the same panel, could be in the same or in 
// different layers
//-----------------------------------------------------------------------------
    for (int i=0; i<2; i++) {
      hit      [i] = HitDoublet->fHit[R->index[i]];
      R->straw [i] = &hit[i]->straw();
      R->rdrift[i] = hit[i]->driftRadius();
      //      layer [i] = straw[i]->id().getLayer();
      //      shId  [i] = straw[i]->index().asInt();

      R->spos  [i] = R->straw[i]->getMidPoint();
      sdir  [i] = R->straw[i]->getDirection();

      phiPanel  = std::atan2(R->spos[i].y(),R->spos[i].x());
      rot.set(-phiPanel);
      
      R->sposr[i]  = rot*R->spos[i];
      sdirr[i]  = rot*sdir[i];
      
      R->tpos [i] = HitDoublet->fTrkPos[R->index[i]];
      tdir [i]    = HitDoublet->fTrkDir[R->index[i]];

      R->tposr[i] = rot*R->tpos[i];
      tdirr[i]  = rot*tdir[i];
      dxdz [i]  = tdirr[i].x()/tdirr[i].z();
    }
//-----------------------------------------------------------------------------
// choose the best combination of the drift signs - the one corresponding 
// to the slope closest to that of the track
// 1. slope chi2 term 
//    for the track dx/dz use average of the two dx/dz slopes 
//    calculated in the two layers corresponding to the doublet hits
//-----------------------------------------------------------------------------
    trkslope  = (dxdz[0]+dxdz[1])/2.;
    findLines(R->sposr,R->rdrift,R->lineSlopes);
      
    for (int is=0; is<4; is++) {
      dsl         = fabs(trkslope - R->lineSlopes[is]);
      xdsl        = dsl/_sigmaSlope;
      R->chi2[is] = xdsl*xdsl;
    }
//-----------------------------------------------------------------------------
// 2. add coordinate term, try to exclude both hits simultaneously
//-----------------------------------------------------------------------------
    const TrkSimpTraj                   *ptraj, *ptraj0;
    vector<KalSite*>::const_iterator   it1, it2;
    double                             s, slen;
    //    bool                               rc;
    HepPoint                           t1pos;
    Hep3Vector                         t1dir;
//-----------------------------------------------------------------------------
// make sure that nothing happens to the track itself
//-----------------------------------------------------------------------------
    s      = (hit[0]->fltLen()+hit[1]->fltLen())/2.;
    ptraj  = KRep->localTrajectory(s,slen);
    ptraj0 = ptraj;

    HelixTraj  smoothed_traj(*ptraj->parameters());
    //    const HelixTraj* psmoothed = &smoothed_traj;

    if (_excludeBothHits == 1) {
//------------------------------------------------------------------------------
// determine the trajectory w/o both hits , then need to calculate residuals
//-----------------------------------------------------------------------------
      std::vector<KalSite*>::const_iterator itt[2];
	
      for (int ih=0; ih<2; ih++) {
	for (std::vector<KalSite*>::const_iterator it=KRep->siteList().begin();
	     it!= KRep->siteList().end(); it++) {
	  const KalHit* kalhit = (*it)->kalHit();
	  if (kalhit && (kalhit->hit() == hit[ih])) {
	    itt[ih] = it;
	    break;
	  }
	}
      }

      KRep->smoothedTraj(itt[0], itt[1], &smoothed_traj);

      for (int ih=0; ih<2; ih++) {
	TrkLineTraj st(HepPoint(R->spos[ih].x(),R->spos[ih].y(),R->spos[ih].z()),sdir[ih],0.,0.); 
//-----------------------------------------------------------------------------
// need to check if I'm close in Z to the hit Z 
//-----------------------------------------------------------------------------
	smoothed_traj.getInfo(hit[ih]->fltLen(),t1pos,t1dir);

	TrkLineTraj tt (t1pos,t1dir,0.); // track, don't need
	
	poca[ih] = TrkPoca(st,0.,tt,0);
	
	sflt[ih] = poca[ih].flt1();
	tflt[ih] = poca[ih].flt2();
	
	st.getInfo(sflt[ih],spi[ih],sdi[ih]);
	tt.getInfo(tflt[ih],tpi[ih],tdi[ih]);
	
	u[ih]    = sdi[ih].cross(tdi[ih]).unit();  // direction towards the center
      }
    }
//-----------------------------------------------------------------------------
// hopefully, this is a helical parameterization around hits of interest,
// so extrapolation uncertainties are not very large
// track parameterization is defined in the global coordinate system, do not rotate
//-----------------------------------------------------------------------------
    for (int ih=0; ih<2; ih++) {
      TrkLineTraj st(HepPoint(R->spos[ih].x(),R->spos[ih].y(),R->spos[ih].z()),sdir[ih],0.,0.); 

      if (_excludeBothHits == 2) {
//-----------------------------------------------------------------------------
// exclude hits one by one - emulate "unbiased" residuals
//-----------------------------------------------------------------------------
	ptraj = ptraj0;
	
	std::vector<KalSite*>::const_iterator itt;
	
	for (std::vector<KalSite*>::const_iterator it=KRep->siteList().begin();
	       it!= KRep->siteList().end(); it++) {
	  const KalHit* kalhit = (*it)->kalHit();
	  if (kalhit && (kalhit->hit() == hit[ih])) {
	    itt = it;
	    break;
	  }
	}
      
	bool rc = KRep->smoothedTraj(itt, itt, &smoothed_traj);
	  
	if (rc == true) {
	  ptraj = &smoothed_traj;
	}
	else {
//-----------------------------------------------------------------------------
// use the default trajectory parameterization
//-----------------------------------------------------------------------------
	  printf(" ERROR in DoubletAmbigResolver::markDoublets: couldnt get Trajectory w/o 1hit, use the local one\n");
	}
      }

      ptraj->getInfo(hit[ih]->fltLen(),t1pos,t1dir);

      TrkLineTraj tt (t1pos,t1dir,0.); // track, don't need
      
      poca[ih] = TrkPoca(st,0.,tt,0);
	
      sflt[ih] = poca[ih].flt1();
      tflt[ih] = poca[ih].flt2();
	
      st.getInfo(sflt[ih],spi[ih],sdi[ih]);
      tt.getInfo(tflt[ih],tpi[ih],tdi[ih]);
	
      u[ih]    = sdi[ih].cross(tdi[ih]).unit();  // direction towards the center

      for (int is=0; is<4; is++) {
//-----------------------------------------------------------------------------
// chi2 contributions of this hit
//-----------------------------------------------------------------------------
	hpos[is][ih]    = spi[ih]+u[ih]*R->rdrift[ih]*_sign[is][ih];
	R->doca[is][ih] = (hpos[is][ih]-tpi[ih]).dot(u[ih]);
	sig             = sqrt(R->rdrift[ih]*R->rdrift[ih] +0.1*0.1); // 2.5; // 1.; // hit[ih]->hitRms();
	xdr[is][ih]     = R->doca[is][ih]/sig;
	R->chi2[is]    += xdr[is][ih]*xdr[is][ih];
      }
    }
//-----------------------------------------------------------------------------
// determine the best solution
//-----------------------------------------------------------------------------
    R->ibest    = -1;
    R->inext    = -1;
    R->chi2min  = 1.e12;
    R->chi2next = 1.e12;
    
    for (int is=0; is<4; is++) {
      if (R->chi2[is] < R->chi2min) {
	R->inext    = R->ibest;
	R->chi2next = R->chi2min;
	R->ibest    = is;
	R->chi2min  = R->chi2[is];
      }
      else if (R->chi2[is] < R->chi2next) {
	R->inext    = is;
	R->chi2next = R->chi2[is];
      }
    }
//-----------------------------------------------------------------------------
// set best solutions
//-----------------------------------------------------------------------------
    int    os                = _sign[R->ibest][0]+_sign[R->ibest][1];
    
    HitDoublet->fOs          = os;
    HitDoublet->fIBest       = R->ibest;
    HitDoublet->fINext       = R->inext;
    HitDoublet->fHitIndex[0] = R->index[0];
    HitDoublet->fHitIndex[1] = R->index[1];

    HitDoublet->fTrkDxDz     = trkslope;
    for (int is=0; is<4; is++) {
      HitDoublet->fDxDz[is]  = R->lineSlopes[is];
      HitDoublet->fChi2[is]  = R->chi2[is];
    }

    for (int i=0; i<2; i++) {
      HitDoublet->fStrawAmbig[R->index[i]] = _sign[R->ibest][i];
    }

    return 0;
  }
						       

//--------------------------------------------------------------------------------
// given a multiplet, resolve the ambiguity for hit: index0 and index1
//--------------------------------------------------------------------------------
  void DoubletAmbigResolver::markDoublet(KalRep* KRep, 
					 Doublet*      HitDoublet, 
					 int           Index0, 
					 int           Index1 ) const {
    mu2e::TrkStrawHit *hit  [2];
//     const mu2e::Straw *straw[2];
   
//     Hep3Vector spos[2], sposr[2], sdir[2], sdirr[2], wpos[10], posPanel;    
//     Hep3Vector tpos[2], tposr[2], tdir[2], tdirr[2];
    
    Hep3Vector wdir;
//     Hep3Vector wdir1, wdir2;

    int                  os; // , ibest , inext;
//     double            rdrift[2], phiPanel;
    
//     int               shId[2];
//     double            trkslope, lineSlopes[4], dxdz[2], chi2[4], doca[4][2], xdr[4][2];
//     double            dsl, xdsl, sig, chi2min, chi2next;

//     double            sflt[2], tflt[2];
//     HepPoint          spi[2] , tpi[2], hpos[4][2];
//     Hep3Vector        sdi[2] , tdi[2], u[2];
//     TrkPoca           poca[2];
//     HepRotationZ      rot;
					// use to pass working data
    Data_t            r;

    //    int               index[2] = {Index0, Index1};
//-----------------------------------------------------------------------------
// calculate doublet parameters
//-----------------------------------------------------------------------------
    r.index[0] = Index0;
    r.index[1] = Index1;

    calculateDoubletParameters(KRep,HitDoublet,&r);
//-----------------------------------------------------------------------------
// create an array with the straw hit indices within a multiplet
//-----------------------------------------------------------------------------
    wdir  = HitDoublet->fShDir;
// //-----------------------------------------------------------------------------
// // by construction, both hits are in the same panel, could be in the same or in 
// // different layers
// //-----------------------------------------------------------------------------
//     for (int i=0; i<2; i++) {
//       hit   [i] = HitDoublet->fHit[index[i]];
//       straw [i] = &hit[i]->straw();
//       //      layer [i] = straw[i]->id().getLayer();
//       rdrift[i] = hit[i]->driftRadius();
//       shId  [i] = straw[i]->index().asInt();

//       spos  [i] = straw[i]->getMidPoint();
//       sdir  [i] = straw[i]->getDirection();

//       phiPanel  = std::atan2(spos[i].y(),spos[i].x());
//       rot.set(-phiPanel);
      
//       sposr[i]  = rot*spos[i];
//       sdirr[i]  = rot*sdir[i];
      
//       tpos [i]  = HitDoublet->fTrkPos[index[i]];
//       tdir [i]  = HitDoublet->fTrkDir[index[i]];

//       tposr[i]  = rot*tpos[i];
//       tdirr[i]  = rot*tdir[i];
//       dxdz [i]  = tdirr[i].x()/tdirr[i].z();
//     }
// //-----------------------------------------------------------------------------
// // choose the best combination of the drift signs - the one corresponding 
// // to the slope closest to that of the track
// // 1. slope chi2 term 
// //    for the track dx/dz use average of the two dx/dz slopes 
// //    calculated in the two layers corresponding to the doublet hits
// //-----------------------------------------------------------------------------
//     trkslope  = (dxdz[0]+dxdz[1])/2.;
//     findLines(sposr,rdrift,lineSlopes);
      
//     for (int is=0; is<4; is++) {
//       dsl       = fabs(trkslope-lineSlopes[is]);
//       xdsl      = dsl/_sigmaSlope;
//       chi2[is]  = xdsl*xdsl;
//     }
// //-----------------------------------------------------------------------------
// // 2. add coordinate term, try to exclude both hits simultaneously
// //-----------------------------------------------------------------------------
//     const TrkSimpTraj                   *ptraj, *ptraj0;
//     vector<KalSite*>::const_iterator   it1, it2;
//     double                             s, slen;
//     //    bool                               rc;
//     HepPoint                           t1pos;
//     Hep3Vector                         t1dir;
    
//     KalRep*  KRep = KRes._krep;

//     s      = (hit[0]->fltLen()+hit[1]->fltLen())/2.;*/
//     ptraj  = KRep->localTrajectory(s,slen);
//     ptraj0 = ptraj;

//     HelixTraj  smoothed_traj(*ptraj->parameters());
//     //    const HelixTraj* psmoothed = &smoothed_traj;

//     if (_excludeBothHits == 1) {
// //------------------------------------------------------------------------------
// // determine the trajectory w/o both hits , then need to calculate residuals
// //-----------------------------------------------------------------------------
//       std::vector<KalSite*>::const_iterator itt[2];
	
//       for (int ih=0; ih<2; ih++) {
// 	for (std::vector<KalSite*>::const_iterator it=KRep->siteList().begin();
// 	     it!= KRep->siteList().end(); it++) {
// 	  const KalHit* kalhit = (*it)->kalHit();
// 	  if (kalhit && (kalhit->hitOnTrack() == hit[ih])) {
// 	    itt[ih] = it;
// 	    break;
// 	  }
// 	}
//       }

//       KRep->smoothedTraj(itt[0], itt[1], &smoothed_traj);

//       for (int ih=0; ih<2; ih++) {
// 	TrkLineTraj st(HepPoint(spos[ih].x(),spos[ih].y(),spos[ih].z()),sdir[ih],0.,0.); 
// //-----------------------------------------------------------------------------
// // need to check if I'm close in Z to the hit Z 
// //-----------------------------------------------------------------------------
// 	smoothed_traj.getInfo(hit[ih]->fltLen(),t1pos,t1dir);

// 	TrkLineTraj tt (t1pos,t1dir,0.); // track, don't need
	
// 	poca[ih] = TrkPoca(st,0.,tt,0);
	
// 	sflt[ih] = poca[ih].flt1();
// 	tflt[ih] = poca[ih].flt2();
	
// 	st.getInfo(sflt[ih],spi[ih],sdi[ih]);
// 	tt.getInfo(tflt[ih],tpi[ih],tdi[ih]);
	
// 	u[ih]    = sdi[ih].cross(tdi[ih]).unit();  // direction towards the center
//       }
//     }
// //-----------------------------------------------------------------------------
// // hopefully, this is a helical parameterization around hits of interest,
// // so extrapolation uncertainties are not very large
// // track parameterization is defined in the global coordinate system, do not rotate
// //-----------------------------------------------------------------------------
//     for (int ih=0; ih<2; ih++) {
//       TrkLineTraj st(HepPoint(spos[ih].x(),spos[ih].y(),spos[ih].z()),sdir[ih],0.,0.); 

//       if (_excludeBothHits == 2) {
// //-----------------------------------------------------------------------------
// // exclude hits one by one - emulate "unbiased" residuals
// //-----------------------------------------------------------------------------
// 	ptraj = ptraj0;
	
// 	std::vector<KalSite*>::const_iterator itt;
	
// 	for (std::vector<KalSite*>::const_iterator it=KRep->siteList().begin();
// 	       it!= KRep->siteList().end(); it++) {
// 	  const KalHit* kalhit = (*it)->kalHit();
// 	  if (kalhit && (kalhit->hitOnTrack() == hit[ih])) {
// 	    itt = it;
// 	    break;
// 	  }
// 	}
      
// 	bool rc = KRep->smoothedTraj(itt, itt, &smoothed_traj);
	  
// 	if (rc == true) {
// 	  ptraj = &smoothed_traj;
// 	}
// 	else {
// //-----------------------------------------------------------------------------
// // use the default trajectory parameterization
// //-----------------------------------------------------------------------------
// 	  printf(" ERROR in DoubletAmbigResolver::markDoublets: couldnt get Trajectory w/o 1hit, use the local one\n");
// 	}
//       }

//       ptraj->getInfo(hit[ih]->fltLen(),t1pos,t1dir);

//       TrkLineTraj tt (t1pos,t1dir,0.); // track, don't need
      
//       poca[ih] = TrkPoca(st,0.,tt,0);
	
//       sflt[ih] = poca[ih].flt1();
//       tflt[ih] = poca[ih].flt2();
	
//       st.getInfo(sflt[ih],spi[ih],sdi[ih]);
//       tt.getInfo(tflt[ih],tpi[ih],tdi[ih]);
	
//       u[ih]    = sdi[ih].cross(tdi[ih]).unit();  // direction towards the center

//       for (int is=0; is<4; is++) {
// //-----------------------------------------------------------------------------
// // chi2 contributions of this hit
// //-----------------------------------------------------------------------------
// 	hpos[is][ih] = spi[ih]+u[ih]*rdrift[ih]*_sign[is][ih];
// 	doca[is][ih] = (hpos[is][ih]-tpi[ih]).dot(u[ih]);
// 	sig          = sqrt(rdrift[ih]*rdrift[ih] +0.1*0.1); // 2.5; // 1.; // hit[ih]->hitRms();
// 	xdr[is][ih]  = doca[is][ih]/sig;
// 	chi2[is]    += xdr[is][ih]*xdr[is][ih];
//       }
//     }
// //-----------------------------------------------------------------------------
// // determine the best solution
// //-----------------------------------------------------------------------------
//     ibest    = -1;
//     inext    = -1;
//     chi2min  = 1.e12;
//     chi2next = 1.e12;
    
//     for (int is=0; is<4; is++) {
//       if (chi2[is] < chi2min) {
// 	inext    = ibest;
// 	chi2next = chi2min;
// 	ibest    = is;
// 	chi2min  = chi2[is];
//       }
//       else if (chi2[is] < chi2next) {
// 	inext    = is;
// 	chi2next = chi2[is];
//       }
//     }
//     HitDoublet->fOs          = os;
//     HitDoublet->fIBest       = ibest;
//     HitDoublet->fINext       = inext;
//     HitDoublet->fHitIndex[0] = index[0];
//     HitDoublet->fHitIndex[1] = index[1];

//     HitDoublet->fTrkDxDz = trkslope;
//     for (int is=0; is<4; is++) {
//       HitDoublet->fDxDz[is] = lineSlopes[is];
//       HitDoublet->fChi2[is] = chi2[is];
//     }
//-----------------------------------------------------------------------------
// set best solutions
//-----------------------------------------------------------------------------
    os    = _sign[r.ibest][0]+_sign[r.ibest][1];
    
    for (int i=0; i<2; i++) {
      hit[i] = HitDoublet->fHit[r.index[i]];
      hit[i]->setAmbigUpdate(false);
//-----------------------------------------------------------------------------
// update the straw hit info inside the doublet, however don't rush 
// to resolve the hit sign ambiguities, do it only when completely sure
// this code is executed after the standard HitAmbigResolver, so when not sure, 
// do nothing and default to HitAmbigResolver
//-----------------------------------------------------------------------------
//      HitDoublet->fStrawAmbig[index[i]] = _sign[ibest][i];
      if (os == 0) {
	if (r.chi2min < _maxDoubletChi2) {
	  if (fabs(r.rdrift[0]+r.rdrift[1]) > 0.8) {
//-----------------------------------------------------------------------------
// OS doublet reliably resolved, reduce the error
//-----------------------------------------------------------------------------
	    if (r.rdrift[i] > _minDriftDoublet) {
//-----------------------------------------------------------------------------
// the hit drift radius is large - reduce the external error
//-----------------------------------------------------------------------------
	      hit[i]->setExtErr(AmbigResolver::_extErr/_scaleErrDoublet);
	      hit[i]->setAmbig(_sign[r.ibest][i]);
	    }
	    else {
//-----------------------------------------------------------------------------
// small drift radius : unless forced, keep the external error large and set 
// the ambiguity to zero to use the wire coordinate
//-----------------------------------------------------------------------------
	      if (_Final == 0) {
		hit[i]->setExtErr(r.rdrift[i]);
		hit[i]->setAmbig(0);
	      }
	      else {
		hit[i]->setAmbig(_sign[r.ibest][i]);
	      }
	    }
	  }
	}
	else {
//-----------------------------------------------------------------------------
// the chi2 is large, consider hits individually
//-----------------------------------------------------------------------------
	  if (_Final == 0) {
	    hit[i]->setExtErr(r.rdrift[i]);
	    hit[i]->setAmbig(0);
	  }
	  else {
//-----------------------------------------------------------------------------
// making final decision
//-----------------------------------------------------------------------------
	    double xr = r.doca[r.ibest][i]/hit[i]->hitErr();
	    if (fabs(xr) > 5.) {
//-----------------------------------------------------------------------------
// hit is very far, reject it
//-----------------------------------------------------------------------------
	      hit[i]->setActivity(false);
	    }
	    else {
//-----------------------------------------------------------------------------
// make the best choice possible, external error should be zero at this point
//-----------------------------------------------------------------------------
	      hit[i]->setExtErr(AmbigResolver::_extErr);
	      hit[i]->setAmbig (_sign[r.ibest][i]);
	    }
	  }
	}
      }
      else {
//-----------------------------------------------------------------------------
// SS doublet
//-----------------------------------------------------------------------------
	if (fabs(r.rdrift[0]-r.rdrift[1]) < _deltaDriftDoublet) {
//-----------------------------------------------------------------------------
// the hardest case - close drift radii
//-----------------------------------------------------------------------------
	  if (r.chi2min < _maxDoubletChi2) {
	    if (r.chi2min/r.chi2next < 0.1) {
//-----------------------------------------------------------------------------
// however, the best chi2 is good enough to be reliable under any circumstances
//-----------------------------------------------------------------------------
	      hit[i]->setExtErr(AmbigResolver::_extErr/_scaleErrDoublet);
	      hit[i]->setAmbig (_sign[r.ibest][i]);
	    }
	    else {
//-----------------------------------------------------------------------------
// the best chi2 is good, but the next one is also close - try to postpone 
// the decision point
//-----------------------------------------------------------------------------
	      if (_Final == 0) {
		double err = fabs(r.rdrift[i]);
		hit[i]->setExtErr(err);
		hit[i]->setAmbig(0);
	      }
	      else {
//-----------------------------------------------------------------------------
// ... but finally decide
//-----------------------------------------------------------------------------
		hit[i]->setExtErr(AmbigResolver::_extErr/_scaleErrDoublet);
		hit[i]->setAmbig (_sign[r.ibest][i]);
	      }
	    }
	  }
	  else {
//-----------------------------------------------------------------------------
// SS doublet with close radii; 
// chi2 doesn't allow to identify the best solution unambiguously, 
// postpone decision for as long as possible 
//-----------------------------------------------------------------------------
	    if (_Final == 0) {
//-----------------------------------------------------------------------------
// don't have to make a decision, scale of uncertainty is defined by the radius
//-----------------------------------------------------------------------------
	      double err = fabs(r.rdrift[i]);
	      hit[i]->setExtErr(err);
	      hit[i]->setAmbig(0);
	    }
	    else {
//-----------------------------------------------------------------------------
// need to make final decision: consider 400 um a limit
//-----------------------------------------------------------------------------
	      if (fabs(r.doca[r.ibest][i]) < 0.4) {
		hit[i]->setExtErr(AmbigResolver::_extErr/_scaleErrDoublet);
		hit[i]->setAmbig (_sign[r.ibest][i]);
	      }
	      else {
		double err = fabs(r.rdrift[i]);
		hit[i]->setExtErr(err);
		hit[i]->setAmbig(0);
	      }
	    }
	  }
	}
	else {
//-----------------------------------------------------------------------------
// SS doublet, the two radii are different 
//-----------------------------------------------------------------------------
	  if (r.chi2min < _maxDoubletChi2) {
//-----------------------------------------------------------------------------
// the best chi2 is good, the doublet drift signs are determined reliably
//-----------------------------------------------------------------------------
	    if (r.chi2min/r.chi2next < 0.2) {
	      if (r.rdrift[i] > _minDriftDoublet) {
		hit[i]->setExtErr(AmbigResolver::_extErr/_scaleErrDoublet);
		hit[i]->setAmbig (_sign[r.ibest][i]);
	      }
	      else {
//-----------------------------------------------------------------------------
// small radius (rdrift[i] < _minDriftDoublet) - use the wire position
// 2015-04-15 P.Murat: for well-resolved doublets it may be possible to decide 
//                     in all cases - need to check
//------------------------------------------------------------*/-----------------
		hit[i]->setExtErr(AmbigResolver::_extErr);
		hit[i]->setAmbig (_sign[r.ibest][i]);
	      }
	    }
	    else {
//-----------------------------------------------------------------------------
// the best and he next chi2's are close, postpone the decision
//-----------------------------------------------------------------------------
	      if (_Final == 0) {
		hit[i]->setExtErr(r.rdrift[i]);
		hit[i]->setAmbig(0);
	      }
	      else {
//-----------------------------------------------------------------------------
// final decision
//-----------------------------------------------------------------------------
		hit[i]->setExtErr(AmbigResolver::_extErr/_scaleErrDoublet);
		hit[i]->setAmbig (_sign[r.ibest][i]);
	      }
	    }
	  }
	  else {
//-----------------------------------------------------------------------------
// the best double chi2 is large - cant believe anything, need to treat hits as 
// separate ones - this is to be implemented yet
// a good example - one of the hits - on Dave's no-gaussial tail
//-----------------------------------------------------------------------------
	    if (_Final == 0) {
	      hit[i]->setExtErr(r.rdrift[i]);
	      hit[i]->setAmbig(0);
	    //	    resolveSingleHit(KRes,hit[i]);
	    }
	    else {
//-----------------------------------------------------------------------------
// _Final = 1: make final decision - forget about 'ibest' - doublet is not good,
//            look at all residuals for this hit
//-----------------------------------------------------------------------------
	      int    best_dd = -1;
	      double max_res = 1.e6;

	      for (int dd=0; dd<4; dd++) {
		if (fabs(r.doca[dd][i]) < max_res) {
		  max_res = fabs(r.doca[dd][i]);
		  best_dd = dd;
		}
	      }

	      double herr = hit[i]->hitErr();
	      if (max_res/herr < 4) {
		hit[i]->setExtErr(AmbigResolver::_extErr/_scaleErrDoublet);
		hit[i]->setAmbig (_sign[best_dd][i]);
	      }
	      else {
		double err = fabs(r.rdrift[i]);
		hit[i]->setExtErr(err);
		hit[i]->setAmbig(0);
					// hit residual is large - try to disable?
		hit[i]->setActivity(false);
	      }
	    }
	  }
	}
      }
    }
    
    if (_debugLevel > 0) {
      for (int i=0; i<2; i++) {
	printf(" %2i %5i %2i %3i %2i %2i %8.3f %8.3f %9.3f %6.3f",
	       i, r.straw[i]->index().asInt(), 
	       HitDoublet->fStationId, HitDoublet->fPanelId, 
	       r.straw[i]->id().getLayer(),
	       r.straw[i]->id().getStraw(),
	       r.spos[i].x(), r.spos[i].y(), r.spos[i].z(),
	       wdir.y()
	       );
	printf(" %6.3f %8.3f %8.3f %9.3f %8.3f %9.3f %6.3f %6.3f",
	       HitDoublet->fTrkDxDz,
	       r.tpos[i].x(),r.tpos[i].y(),r.tpos[i].z(),
	       r.sposr[i].x(),r.tposr[i].x(),r.doca[r.ibest][i],r.rdrift[i]
	       ); 
	printf(" %2i %6.3f %6.3f %8.2e %6.3f %6.3f %8.2e %6.3f %6.3f %8.2e %6.3f %6.3f %8.2e\n",
	       _sign[r.ibest][i],
	       r.lineSlopes[0], r.doca[0][i], r.chi2[0], 
	       r.lineSlopes[1], r.doca[1][i], r.chi2[1], 
	       r.lineSlopes[2], r.doca[2][i], r.chi2[2], 
	       r.lineSlopes[3], r.doca[3][i], r.chi2[3]
	       );
      }
    }
  }


//---------------------------------------------------------------------------
// loop over the doublets found and mark their ambiguities
//---------------------------------------------------------------------------
  void DoubletAmbigResolver::markMultiplets (KalRep* Krep, std::vector<Doublet>* dcol ) const {

    mu2e::TrkStrawHit    *hit, *hitj, *hitk;
    mu2e::Straw    straw;
    Doublet              *doublet;
    int                  ndoublets, ndhits;

    ndoublets  = dcol->size();

    if (_debugLevel > 0) {
      printf("[KalFitHack::markMultiplets] BEGIN iherr:%i\n",_iter);
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
      ndhits  = doublet->fNStrawHits;
      
      if (ndhits == 1) {
//-----------------------------------------------------------------------------
// a single hit in the plane - keep its external error large, unless the drift sign 
// can be determined reliably 
//-----------------------------------------------------------------------------
	hit = doublet->fHit[0];

	resolveSingleHit(Krep,hit);
      }
      else if (ndhits == 2) {
//-----------------------------------------------------------------------------
// 2 hits in a panel - attempt to determine the drift signs
//-----------------------------------------------------------------------------
	markDoublet(Krep,doublet,0,1);
      }
      else {
//-----------------------------------------------------------------------------
// more than 2 hits in a panel
//-----------------------------------------------------------------------------
	int      tmpLayerId, layer0, layer1, jbest(-1), kbest(-1), imax;
	int      tmpId(-1), id0(-1), id1(-1), nd;
	double   rdrift, chi2_d, chi2_best (1.e12);
	Doublet  bd, *d;
	vector<Doublet> list;

	for (int j=0; j<ndhits; ++j) {
	  for (int k=j+1; k<ndhits; ++k){
//-----------------------------------------------------------------------------
// P.Murat: 
// learn logic: logic here looks questionalble, but let's first figure what it is exactly 
// - use the first found OS doublet, w/o looking at the quality
// - after it is found, resolve drift signs right away for all hits in the multiplet
// - this could be dangerous, especially, in presence of the background 
// change logic: 
//-----------------------------------------------------------------------------
	    markDoublet(Krep,doublet,j,k); 

	    list.push_back(*doublet);

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
	  hitj    = doublet->fHit[jbest];
	  straw   = hitj->straw();
	  layer0  = straw.id().getLayer();
	  id0     = straw.index().asInt();

	  hitk    = doublet->fHit[kbest];
	  straw   = hitk->straw();
	  layer1  = straw.id().getLayer();
	  id1     = straw.index().asInt();
//-----------------------------------------------------------------------------
// out of the two hits making the best doublet, select the one with the larger 
// radius
//-----------------------------------------------------------------------------
	  if (hitj->driftRadius() > hitk->driftRadius()) imax = jbest;
	  else                                           imax = kbest;

	  for (int ih=0; ih<ndhits; ++ih) {
	    hit        = doublet->fHit[ih];
	    Straw const& straw      = hit->straw();
	    tmpLayerId = straw.id().getLayer();
	    tmpId      = straw.index().asInt();

	    if (ih == jbest) {
	      hit->setAmbig(bd.fStrawAmbig[0]);
	      hit->setAmbigUpdate(false);
	    }
	    else if (ih == kbest) {
	      hit->setAmbig(bd.fStrawAmbig[1]);
	      hit->setAmbigUpdate(false);
	    }
	    else {
//-----------------------------------------------------------------------------
// the assumption here is that in case of a triplet two hits in the same layer 
// can't have the same drift sign - which is not necessarily correct.
// *stick to it for the time being*
//-----------------------------------------------------------------------------
	      if (tmpLayerId == layer0) {
		if (tmpId != id0) {
		  doublet->fStrawAmbig[ih] = -doublet->fStrawAmbig[jbest];
		}
	      }
	      else if (tmpLayerId == layer1) {
		if (tmpId != id1) {
		  doublet->fStrawAmbig[ih] = -doublet->fStrawAmbig[kbest];
		}
	      }
	      else {
//-----------------------------------------------------------------------------
// 2 best hits in the same layer, this hit is in the different one 
//-----------------------------------------------------------------------------
		nd = list.size();
		for (int i=0; i<nd; i++) {
		  d = &list[i];
		  if (d->fHitIndex[0] == ih) {
		    if (d->fHitIndex[1] == imax) {
		      doublet->fStrawAmbig[ih] = d->fStrawAmbig[0];
		    }
		  }
		  else if (d->fHitIndex[1] == ih) {
		    if (d->fHitIndex[0] == imax) {
		      doublet->fStrawAmbig[ih] = d->fStrawAmbig[1];
		    }
		  }
		}
	      }
//-----------------------------------------------------------------------------
// if the hit radius is greated than some minimal value, set the drift sign
//-----------------------------------------------------------------------------
	      rdrift = hit->driftRadius();

	      if ( fabs(rdrift) > _minDriftDoublet){
		hit->setAmbig(doublet->fStrawAmbig[ih]);
		hit->setAmbigUpdate(false);
	      }
	      else {
		if (_Final == 0) {
		  hit->setAmbig(0);
		  hit->setExtErr(rdrift);
		  hit->setAmbigUpdate(false);
		}
		else {
//-----------------------------------------------------------------------------
// final decision
//-----------------------------------------------------------------------------
		  hit->setAmbig(doublet->fStrawAmbig[ih]);
		  hit->setAmbigUpdate(false);
		}
	      }
	    }
	  }   
	}
      }
    }
  }

//-----------------------------------------------------------------------------
// build list of doublets and resolve the drift signs
//-----------------------------------------------------------------------------
  void DoubletAmbigResolver::resolveTrk(KalRep* KRep) const {

					// initialize external hit errors
    initHitErrors(KRep);
    vector<Doublet> listOfDoublets;
    findDoublets (KRep,&listOfDoublets);

    if (listOfDoublets.size() > 0) markMultiplets(KRep,&listOfDoublets);
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

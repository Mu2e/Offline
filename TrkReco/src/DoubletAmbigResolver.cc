///////////////////////////////////////////////////////////////////////////////
// class to resolve hit ambiguities one hit at a time,
// assuming a reasonable track fit as input
// DoubletAmbigResolver is instantiated from Stntuple/mod/InitTrackBlock.cc
// thus all parameters need to have defaults in the source
///////////////////////////////////////////////////////////////////////////////
/*
#include "Offline/TrkReco/inc/DoubletAmbigResolver.hh"
#include "Offline/BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalSite.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include <vector>
#include <algorithm>
#include <functional>

using CLHEP::Hep3Vector;
using CLHEP::HepRotationZ;

namespace mu2e {
  typedef std::vector<TrkStrawHit*>::iterator TSHI;


  //-----------------------------------------------------------------------------
  DoubletAmbigResolver::DoubletAmbigResolver(fhicl::ParameterSet const& PSet,
      double     extErr,
      int                        Iter,
      int      Final) :
    AmbigResolver(extErr),
    _debugLevel  (PSet.get<int>   ("debugLevel"      )),
    _mindrift    (PSet.get<double>("HitMinDrift"     )),
    _zeropenalty (PSet.get<double>("ZeroDriftPenalty")),
    _penalty     (PSet.get<bool>  ("HitAmbigPenalty" )),
    //-----------------------------------------------------------------------------
    // not sure what the next four numbers are
    //-----------------------------------------------------------------------------
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
    _sigmaSlope       (PSet.get<double>("sigmaSlope"       )),
    _maxDoubletChi2   (PSet.get<double>("maxDoubletChi2"   )),
    _scaleErrDoublet  (PSet.get<double>("scaleErrDoublet"  )),
    _minDriftDoublet  (PSet.get<double>("minDriftDoublet"  )),
    _deltaDriftDoublet(PSet.get<double>("deltaDriftDoublet")),
    _excludeBothHits  (PSet.get<int>   ("excludeBothHits"  )),             // default:1
    _minChi2Ratio     (PSet.get<double>("minChi2Ratio"     )),
    _tempScale        (PSet.get<double>("tempScale"        )),
    _penaltyScale     (PSet.get<double>("penaltyScale"     )),
    _useMeanResidual  (PSet.get<double>("useMeanResidual"  )),
    _maxMeanResidual  (PSet.get<double>("maxMeanResidual"  )),
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
      printf("[DoubletAmbigResolver::findDoublets]-------------------------------------------------\n");
      printf("[DoubletAmbigResolver::findDoublets]  i  shId  ch  panel  il   iw   driftR       doca\n");
      printf("[DoubletAmbigResolver::findDoublets]-------------------------------------------------\n");
    }

    dcol = DCol;
    dcol->clear();

    int multipletIndex(0);
    //-----------------------------------------------------------------------------
    // use all hits, not only active ones
    //-----------------------------------------------------------------------------
    int i = 0;
    TrkStrawHitVector tshv;
    convert(KRep->hitVector(),tshv);
    for (auto ihit=tshv.begin(); ihit!=tshv.end(); ++ihit) {
      //      if (hit->isActive() == 0)                             goto END_OF_LOOP;
      Straw const& straw     = (*ihit) ->straw();
      wdir      = straw.getDirection();
      pos       = straw.getMidPoint();
      station   = straw.id().getPlane();
      panel     = straw.id().getPanel();
      shId      = straw.id().asUint16();
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
        printf("[DoubletAmbigResolver::findDoublets] %2i  %5i %3i  %4i  %3i %3i %8.3f %8.3f\n",
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
        ++idlast;
        ++multipletIndex;
      }
      else {
        if (panel == oldPanel) {
          //-----------------------------------------------------------------------------
          // same chamber, same panel : add one more hit to the last doublet
          //-----------------------------------------------------------------------------
          dcol->at(idlast-1).addStrawHit(tdir, trkpos, *ihit);
        }
        else {
          //-----------------------------------------------------------------------------
          // same chamber, different panel : new doublet candidate
          //-----------------------------------------------------------------------------
          dcol->push_back(Doublet(multipletIndex, station, panel, wdir, tdir, trkpos, *ihit)); //
          oldStation = station;
          oldPanel   = panel;
          ++idlast;
          ++multipletIndex;
        }
      }
      i += 1;
    }
    //-----------------------------------------------------------------------------
    // list of doublets is formed, determine the doublet parameters but not
    // assine the hit drift directions
    //-----------------------------------------------------------------------------
    int ndoublets = dcol->size();
    Data_t r;

    for (int i=0; i<ndoublets; i++) {
      Doublet* d = &dcol->at(i);
      int nhits = d->fNStrawHits;
      if (nhits == 1) {
        //-----------------------------------------------------------------------------
        // single hit
        //-----------------------------------------------------------------------------

      }
      else if (nhits == 2) {
        //-----------------------------------------------------------------------------
        // regular case: a 2-hit doublet
        //-----------------------------------------------------------------------------
        r.index[0] = 0;
        r.index[1] = 1;
        calculateDoubletParameters(KRep,d,&r);
      }
      else {
        //-----------------------------------------------------------------------------
        // a multiplet - more than 2 hits in a panel
        //-----------------------------------------------------------------------------
        //  vector<Doublet> list;
        Doublet         bd;
        double          chi2_d, chi2_best(1.e12);

        for (int j=0; j<d->fNStrawHits; ++j) {
          for (int k=j+1; k<d->fNStrawHits; ++k) {
            //--------------------------------------------------------------------------
            // find the best combination
            //-----------------------------------------------------------------------------
            r.index[0] = j;
            r.index[1] = k;
            calculateDoubletParameters(KRep,d,&r);

            chi2_d = d->chi2Best();
            if (chi2_d < chi2_best) {
              chi2_best = chi2_d;
              bd        = *d;
            }
          }
        }
        *d =bd;
      }
    }
    //-----------------------------------------------------------------------------
    // list of doublets is formed, the rest of this routine - diagnostics only
    //-----------------------------------------------------------------------------
    Doublet             *doublet;
    CLHEP::HepRotationZ rot;

    if (_debugLevel > 1) {

      printf("[DoubletAmbigResolver::findDoublets] _iter:%i: found %i multiplets\n",_iter,ndoublets);
      printf("--------------------------------------------------------------");
      printf("------------------------------------------------------------------------\n");
      printf("  i  shId pl pn il  is      x        y         z      sinphi  ");
      printf("tksphi    xtrk     ytrk      ztrk      xr      yr     zr       doca   rdr\n");
      printf("--------------------------------------------------------------");
      printf("------------------------------------------------------------------------\n");

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
          shId  = straw.id().asUint16();
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

          printf(" %2i %5i %2i %2i %2i %3i %9.3f %9.3f %9.3f  %6.3f  ",
              i, shId, doublet->fStationId, doublet->fPanelId,
              straw.id().getLayer(),
              straw.id().getStraw(),
              straw.getMidPoint().x(),
              straw.getMidPoint().y(),
              straw.getMidPoint().z(),
              wdir.y()
              );

          printf("%6.3f %9.3f %8.3f %9.3f %8.3f %4.1f %9.3f %7.3f %5.3f\n",
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
  void DoubletAmbigResolver::resolveSingleHit(KalRep*            Krep,
      mu2e::TrkStrawHit* Hit ) const {

    double                     doca[2], xbest, xnext, perr;
    std::vector<TrkStrawHit*>  hits;

    Straw const& straw = Hit->straw();

    const Hep3Vector& wdir = straw.getDirection();
    const Hep3Vector& wmid = straw.getMidPoint();

    if (_useMeanResidual) perr = _meanResidual;
    else                  perr = fabs(Hit->driftRadius()/sqrt(3));
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
      double r = Hit->driftRadius();

      doca[0] = poca.doca()-r;
      doca[1] = poca.doca()+r;

      double err = Hit->temperature()*0.0625;    // "external" error
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

      double drift_err(0.220);   // assume 220 microns - from

      xbest = doca[ibest]/drift_err;
      xnext = doca[inext]/drift_err;
      //-----------------------------------------------------------------------------
      // want the best solution to be consistent with the trajectory, another one -
      // to be inconsistent, and the two - significantly different
      //-----------------------------------------------------------------------------
      if ((fabs(xbest) < 5.) && (fabs(xnext) > 5) && (fabs(xbest/xnext) < 0.25)) {
        //-----------------------------------------------------------------------------
        // hit is close enough to the trajectory
        //-----------------------------------------------------------------------------
        if (ibest == 0) Hit->setAmbig( 1);
        else            Hit->setAmbig(-1);
      }
      else {
        //-----------------------------------------------------------------------------
        // can't tell which one is better, just set ambiguity to zero
        //-----------------------------------------------------------------------------
        if (_Final == 0) {
          //    Hit->setTemperature(Hit->driftRadius()/Hit->driftVelocity());
          Hit->setPenalty(perr);
          Hit->setAmbig(0);
        }
        else {
          //-----------------------------------------------------------------------------
          // no point to keep hits with very large residuals
          //-----------------------------------------------------------------------------
          if (ibest == 0) Hit->setAmbig( 1);
          else            Hit->setAmbig(-1);
          Hit->setPenalty(perr);
        }
      }
      if (_debugLevel > 0) {
        CLHEP::Hep3Vector hpos;
        Hit->hitPosition(hpos);
        printf(" %2i %5i %2i %2i %2i %2i %8.3f %8.3f %9.3f, %2i %2i %8.3f %8.3f\n",
            -1, Hit->comboHit().strawId().asUint16(),
            straw.id().getPlane(),
            straw.id().getPanel(),
            straw.id().getLayer(),
            straw.id().getStraw(),
            hpos.x(), hpos.y(), hpos.z(), ibest, inext, xbest, xnext);
      }
    }
    else {
      //-----------------------------------------------------------------------------
      // couldn't determine doca
      //-----------------------------------------------------------------------------
      printf(" ***** DoubletAmbigResolver::resolveSingleHit warning: cant determine POCA\n");
      //      Hit->setTemperature(Hit->driftRadius()/Hit->driftVelocity());
      Hit->setPenalty(perr);
      Hit->setAmbig(0);
    }

  }

  //-----------------------------------------------------------------------------
  // 2015-02-25 P.Murat: new resolver
  // assume that coordinates are rotated into the coordinate system where Y axis
  // is pointed along the wire by a rotation along the Z axis
  // sign ordering convention:    ++, +-, --, -+ is defined by DoubletAmbigResolver::_sign
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
      alpha     = (R[1]*_sign[i][1]-R[0]*_sign[i][0])/dr;
      nx[i]     = _sign[i][0]*(-ly*alpha+lx*sqrt(1-alpha*alpha));
      ny[i]     = _sign[i][0]*( lx*alpha+ly*sqrt(1-alpha*alpha));
      Slopes[i] = ny[i]/nx[i];
    }
  }


  //-----------------------------------------------------------------------------
  // assume that the two hits to be used are known,their indices are stored in R
  //-----------------------------------------------------------------------------
  int DoubletAmbigResolver::calculateDoubletParameters(const KalRep* KRep      ,
      Doublet*      HitDoublet,
      Data_t*       R         ) const {

    mu2e::TrkStrawHit *hit  [2];

    Hep3Vector sdir[2], sdirr[2], wpos[10], posPanel;
    Hep3Vector tdir[2], tdirr[2];

    Hep3Vector wdir1, wdir2;

    double            phiPanel;

    double            trkslope, dxdz[2], xdr[4][2];
    double            dsl, xdsl, sig;

    double            sflt[2], tflt[2];
    HepPoint          spi[2] , tpi[2], hpos[4][2];
    Hep3Vector        sdi[2] , tdi[2], u[2];
    TrkPoca           poca[2];
    HepRotationZ      rot;
    //-----------------------------------------------------------------------------
    // by construction, both hits are in the same panel, could be in the same or in
    // different layers
    // dx/dz[i] === dR/dZ of the track in the wire plane of the i-th hit
    //-----------------------------------------------------------------------------
    for (int i=0; i<2; i++) {
      hit      [i] = HitDoublet->fHit[R->index[i]];
      R->straw [i] = &hit[i]->straw();
      R->rdrift[i] = hit[i]->driftRadius();

      R->spos  [i] = R->straw[i]->getMidPoint();
      sdir     [i] = R->straw[i]->getDirection();

      phiPanel  = std::atan2(R->spos[i].y(),R->spos[i].x());
      rot.set(-phiPanel);

      R->sposr[i]  = rot*R->spos[i];
      sdirr   [i]  = rot*sdir[i];

      R->tpos [i] = HitDoublet->fTrkPos[R->index[i]];
      tdir    [i] = HitDoublet->fTrkDir[R->index[i]];

      R->tposr[i] = rot*R->tpos[i];
      tdirr[i]    = rot*tdir[i];
      dxdz [i]    = tdirr[i].x()/tdirr[i].z();
    }
    //-----------------------------------------------------------------------------
    // 2019-02-03 PM : hope, this is temporary: make sure that the two hits are different!
    // today it is not given
    //-----------------------------------------------------------------------------
    if (hit[0] == hit[1]) {
      printf(" ERROR detected in DoubletAmbigResolver::calculateDoubletParameters: doublet with two identical hits. BAIL out\n");
      R->ibest = -1;
      return -1;
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
      dsl              = fabs(trkslope - R->lineSlopes[is]);
      xdsl             = dsl/_sigmaSlope;
      R->chi2Slope[is] = xdsl*xdsl;
    }
    //-----------------------------------------------------------------------------
    // 2. add coordinate term, try to exclude both hits simultaneously
    //-----------------------------------------------------------------------------
    const TrkSimpTraj                   *ptraj, *ptraj0;
    vector<KalSite*>::const_iterator   it1, it2;
    double                             s, slen;
    HepPoint                           t1pos;
    Hep3Vector                         t1dir;
    //-----------------------------------------------------------------------------
    // make sure that nothing happens to the track itself
    //-----------------------------------------------------------------------------
    s      = (hit[0]->fltLen()+hit[1]->fltLen())/2.;
    ptraj  = KRep->localTrajectory(s,slen);
    ptraj0 = ptraj;

    HelixTraj  smoothed_traj(*ptraj->parameters());

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
        // error of the coordinate term is questionable - at initial iterations hits
        // with small drift radius have smaller errors - doesn't make sense
        //-----------------------------------------------------------------------------
        hpos[is][ih]     = spi[ih]+u[ih]*R->rdrift[ih]*_sign[is][ih];

        R->doca[is][ih]  = static_cast<Hep3Vector>(hpos[is][ih]-tpi[ih]).dot(u[ih]);

        // sig             = sqrt(R->rdrift[ih]*R->rdrift[ih] +0.1*0.1); // 2.5; // 1.; // before 2019-01-20
        sig              = hit[ih]->totalErr();

        xdr[is][ih]      = R->doca[is][ih]/sig;
        R->chi2Coord[is] = xdr[is][ih]*xdr[is][ih];
        R->chi2     [is] = R->chi2Slope[is]+R->chi2Coord[is];
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
      HitDoublet->fChi2Slope[is]  = R->chi2Slope[is];
      HitDoublet->fChi2Coord[is]  = R->chi2Coord[is];
      HitDoublet->fChi2     [is]  = R->chi2     [is];
    }

    for (int i=0; i<2; i++) {
      HitDoublet->fStrawAmbig[R->index[i]] = _sign[R->ibest][i];
    }

    return 0;
  }


  //-----------------------------------------------------------------------------
  // define hit drift sign in case the two best chi2's are close
  //-----------------------------------------------------------------------------
  void DoubletAmbigResolver::defineHitDriftSign(mu2e::TrkStrawHit* Hit, int I, Data_t* R) const {
    double perr;

    //-----------------------------------------------------------------------------
    // initialize the penalty error
    //-----------------------------------------------------------------------------
    if (_useMeanResidual) perr = _meanResidual;
    else                  perr = fabs(R->rdrift[I]);  // do we need 1./sqrt(3) here?

    if (_sign[R->ibest][I] == _sign[R->inext][I]) {
      //-----------------------------------------------------------------------------
      // both the 'best' and 'next' chi2s correspond to the same drift sign of hit 'i'
      // check other two chi2's. If both are an order of magnitude worse than the best,
      // assume that for this hit the drift direction is defined
      //-----------------------------------------------------------------------------
      int drift_dir_defined = 1;
      for (int j=0; j<4; j++) {
        if ((j != R->ibest) && (j != R->inext)) {
          if (R->chi2min/R->chi2[j] > 0.2) {
            drift_dir_defined = 0;
            break;
          }
        }
      }

      if (drift_dir_defined == 1) {
        Hit->setAmbig (_sign[R->ibest][I]);
      }
      else {
        //-----------------------------------------------------------------------------
        // all 4 chi2s are "of the same order", nothing to fish
        //-----------------------------------------------------------------------------
        Hit->setAmbig  (0  );
        Hit->setPenalty(perr);
      }
    }
    else {
      //-----------------------------------------------------------------------------
      // the best and the next solutions correspond to different drift directions
      // of this hit. Can't choose, assign large error
      //-----------------------------------------------------------------------------
      Hit->setAmbig  (0  );
      Hit->setPenalty(perr);
    }
  }


  //--------------------------------------------------------------------------------
  // given a multiplet, resolve ambiguities for a doublet made out of hits
  // 'Index0' and 'Index1'
  //--------------------------------------------------------------------------------
  void DoubletAmbigResolver::markDoublet(KalRep*       KRep      ,
      Doublet*      HitDoublet,
      int           Index0    ,
      int           Index1    ) const {
    mu2e::TrkStrawHit *hit  [2];
    Hep3Vector        wdir;
    int               os;
    // use to pass working data
    Data_t            r;
    //-----------------------------------------------------------------------------
    // calculate doublet parameters
    //-----------------------------------------------------------------------------
    r.index[0] = Index0;
    r.index[1] = Index1;

    int rc = calculateDoubletParameters(KRep,HitDoublet,&r);
    if (rc < 0) {
      //-----------------------------------------------------------------------------
      // an error detected, bail out
      //-----------------------------------------------------------------------------
      HitDoublet->fHit[Index0]->setAmbig(0);
      HitDoublet->fHit[Index0]->setPenalty(1.e5);
      HitDoublet->fHit[Index1]->setAmbig(0);
      HitDoublet->fHit[Index1]->setPenalty(1.e5);
      return;
    }
    //-----------------------------------------------------------------------------
    // create an array with the straw hit indices within a multiplet
    //-----------------------------------------------------------------------------
    wdir  = HitDoublet->fShDir;
    //-----------------------------------------------------------------------------
    // set best solutions
    //-----------------------------------------------------------------------------
    os    = _sign[r.ibest][0]+_sign[r.ibest][1];

    for (int i=0; i<2; i++) {
      hit[i] = HitDoublet->fHit[r.index[i]];
      //-----------------------------------------------------------------------------
      // os = 0: an opposite sign doublet
      // update the straw hit info inside the doublet, however don't rush
      // to resolve the hit sign ambiguities, do it only when completely sure
      // this code is executed after the standard HitAmbigResolver, so when not sure,
      // do nothing and default to HitAmbigResolver
      //-----------------------------------------------------------------------------
      if (os == 0) {
        if (r.chi2min < _maxDoubletChi2) {
          //    if (fabs(r.rdrift[0]+r.rdrift[1]) > 0.8) {
          //-----------------------------------------------------------------------------
          // OS doublet reliably resolved, reduce the error
          //-----------------------------------------------------------------------------
          if (r.rdrift[i] > _minDriftDoublet) {
            //-----------------------------------------------------------------------------
            // the hit drift radius is large - reduce the external error (no penalty in this case)
            //-----------------------------------------------------------------------------
            hit[i]->setAmbig(_sign[r.ibest][i]);
            hit[i]->setPenalty(0);
          }
          else {
            //-----------------------------------------------------------------------------
            // small drift radius : unless forced, keep the external error large and set
            // the ambiguity to zero to use the wire coordinate
            //-----------------------------------------------------------------------------
            if (_Final == 0) {
              hit[i]->setAmbig(0);
              if (_useMeanResidual) hit[i]->setPenalty(_meanResidual);
              else                  hit[i]->setPenalty(r.rdrift[i]);
            }
            else {
              hit[i]->setAmbig(_sign[r.ibest][i]);
            }
          }
          //    }
        }
        else {
          //-----------------------------------------------------------------------------
          // the chi2 is large, consider hits individually:
          // however, the chi2 may be large because of one of the hit times being misreconstructed,
          // while the other hit drift sign still could be reconstructed
          //-----------------------------------------------------------------------------
          if (_Final == 0) {
            defineHitDriftSign(hit[i],i,&r);
          }
          else {
            //-----------------------------------------------------------------------------
            // making final decision
            //-----------------------------------------------------------------------------
            double xr = r.doca[r.ibest][i]/hit[i]->totalErr();
            if (fabs(xr) > 5.) {
              //-----------------------------------------------------------------------------
              // hit is very far (more than by 5 "sigma"), reject it
              //-----------------------------------------------------------------------------
              hit[i]->setActivity(false);
            }
            else {
              //-----------------------------------------------------------------------------
              // make the best choice possible, external error should be zero at this point
              //-----------------------------------------------------------------------------
              hit[i]->setPenalty(0);
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
            if (r.chi2min/r.chi2next < _minChi2Ratio) { // 0.3
              //-----------------------------------------------------------------------------
              // however, the best chi2 is good enough to be reliable under any circumstances
              //-----------------------------------------------------------------------------
              hit[i]->setPenalty(0);
              hit[i]->setAmbig (_sign[r.ibest][i]);
            }
            else {
              //-----------------------------------------------------------------------------
              // the best chi2 is good, but the next one is also close - try to postpone
              // the decision point
              //-----------------------------------------------------------------------------
              if (_Final == 0) {
                defineHitDriftSign(hit[i],i,&r);
              }
              else {
                //-----------------------------------------------------------------------------
                // ... but finally decide
                //-----------------------------------------------------------------------------
                hit[i]->setPenalty(0);
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
              if (_useMeanResidual) hit[i]->setPenalty(_meanResidual);
              else                  hit[i]->setPenalty(fabs(r.rdrift[i]));
              hit[i]->setAmbig(0);
            }
            else {
              //-----------------------------------------------------------------------------
              // need to make final decision: consider 400 um a limit
              //-----------------------------------------------------------------------------
              if (fabs(r.doca[r.ibest][i]) < 0.4) {
                hit[i]->setAmbig (_sign[r.ibest][i]);
                hit[i]->setPenalty(0);
              }
              else {
                hit[i]->setAmbig(0);
                if (_useMeanResidual) hit[i]->setPenalty(_meanResidual);
                else                  hit[i]->setPenalty(fabs(r.rdrift[i])/sqrt(3.));
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
            if (r.chi2min/r.chi2next < _minChi2Ratio) {
              if (r.rdrift[i] > _minDriftDoublet) {
                hit[i]->setPenalty(0);
                hit[i]->setAmbig (_sign[r.ibest][i]);
              }
              else {
                //-----------------------------------------------------------------------------
                // small radius (rdrift[i] < _minDriftDoublet) - use the wire position
                // 2015-04-15 P.Murat: for well-resolved doublets it may be possible to decide
                //                     in all cases - need to check
                //-----------------------------------------------------------------------------
                if (_useMeanResidual) hit[i]->setPenalty(_meanResidual);
                else                  hit[i]->setPenalty(AmbigResolver::_tmpErr/2/sqrt(3));
                hit[i]->setAmbig (_sign[r.ibest][i]);
              }
            }
            else {
              //-----------------------------------------------------------------------------
              // the best and the next chi2's are close
              //-----------------------------------------------------------------------------
              if (_Final == 0) {
                //-----------------------------------------------------------------------------
                // postpone the decision for as long as possible
                //-----------------------------------------------------------------------------
                defineHitDriftSign(hit[i],i,&r);
              }
              else {
                //-----------------------------------------------------------------------------
                // final decision
                //-----------------------------------------------------------------------------
                hit[i]->setPenalty(0);
                hit[i]->setAmbig (_sign[r.ibest][i]);
              }
            }
          }
          else {
            //-----------------------------------------------------------------------------
            // the best doublet chi2 is large - cant believe anything, need to treat hits as
            // separate ones - this is to be implemented yet
            // a good example - one of the hits - on Dave's no-gaussial tail
            //-----------------------------------------------------------------------------
            if (_Final == 0) {
              defineHitDriftSign(hit[i],i,&r);
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

              double herr = hit[i]->totalErr();
              if (max_res/herr < 4) {
                hit[i]->setPenalty(0);
                hit[i]->setAmbig (_sign[best_dd][i]);
              }
              else {
                if (_useMeanResidual) hit[i]->setPenalty(_meanResidual);
                else                  hit[i]->setPenalty(fabs(r.rdrift[i]));

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
        printf(" %2i %5i %2i %2i %2i %2i %8.3f %8.3f %9.3f %6.3f",
            i, r.straw[i]->id().asUint16(),
            HitDoublet->fStationId, HitDoublet->fPanelId,
            r.straw[i]->id().getLayer(),
            r.straw[i]->id().getStraw(),
            r.spos[i].x(), r.spos[i].y(), r.spos[i].z(),
            wdir.y()
            );
        printf(" %6.3f %8.3f %8.3f %9.3f %8.3f %9.3f %8.3f %6.3f",
            HitDoublet->fTrkDxDz,
            r.tpos[i].x(),r.tpos[i].y(),r.tpos[i].z(),
            r.sposr[i].x(),r.tposr[i].x(),r.doca[r.ibest][i],r.rdrift[i]
            );
        printf(" %2i %7.3f %8.3f %8.2e %7.3f %8.3f %8.2e %7.3f %8.3f %8.2e %7.3f %8.3f %8.2e\n",
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

    mu2e::TrkStrawHit    *hit, *hitj, *hitk, *hmax;
    mu2e::Straw          straw;
    Doublet              *doublet;
    int                  ndoublets, ndhits;

    ndoublets  = dcol->size();

    if (_debugLevel > 0) {
      printf("[DoubletAmbigResolver::markMultiplets] BEGIN _iter:%i\n",_iter);
      printf("------------------------------------------------------");
      printf("----------------------------------------------------------------------------------------");
      printf("------------------------------------------------------------------------------------------\n");
      printf("  i  shId pl pn il is      x       y        z     ");
      printf(" w_ny    tdxdz   xtrk     ytrk      ztrk       xr      xtrkr    doca    rdr  am  ");
      printf("  s++    doca++   chi2++    s+-    doca+-   chi2+-    s--    doca--   chi2--    s-+    doca-+   chi2-+\n");
      printf("------------------------------------------------------");
      printf("----------------------------------------------------------------------------------------");
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
        int      jbest(-1), kbest(-1) ; // imax; // tmpLayerId, layer0, layer1,
        int      nd; // tmpId(-1), id0(-1), id1(-1)
        double   rdrift, chi2_d, chi2_best (1.e12);
        Doublet  bd, ad, *d;
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
            ad    = *doublet;
            markDoublet(Krep,&ad,j,k);

            list.push_back(ad);

            chi2_d = ad.chi2Best();
            if (chi2_d < chi2_best) {
              jbest     = j;
              kbest     = k;
              chi2_best = chi2_d;
              bd        = ad;
            }
          }
        }

        *doublet = bd;
        if (chi2_best < _maxDoubletChi2) {
          //-----------------------------------------------------------------------------
          // the "best" doublet is good enough, resolve drift signs for the rest hits
          //-----------------------------------------------------------------------------
          hitj    = doublet->fHit[jbest];
          hitk    = doublet->fHit[kbest];
          //-----------------------------------------------------------------------------
          // out of the two hits making the best doublet, select the one with the larger
          // radius
          //-----------------------------------------------------------------------------
          if (hitj->driftRadius() > hitk->driftRadius()) hmax = hitj;
          else                                           hmax = hitk;

          // loop over hits, assign drift signs
          for (int ih=0; ih<ndhits; ++ih) {
            hit = doublet->fHit[ih];

            if (ih == jbest) {
              hit->setAmbig(_sign[bd.fIBest][0]);
            }
            else if (ih == kbest) {
              hit->setAmbig(_sign[bd.fIBest][1]);
            }
            else {
              //-----------------------------------------------------------------------------
              // look for the doublet containing 'hit' and 'hmax'
              //-----------------------------------------------------------------------------
              nd = list.size();
              for (int i=0; i<nd; i++) {
                d = &list[i];
                if (hit == d->fHit[d->fHitIndex[0]]) {
                  if (hmax == d->fHit[d->fHitIndex[1]]) {
                    doublet->fStrawAmbig[ih] = _sign[d->fIBest][0];
                    break;
                  }
                }
                else if (hit == d->fHit[d->fHitIndex[1]]) {
                  if (hmax == d->fHit[d->fHitIndex[0]]) {
                    doublet->fStrawAmbig[ih] = _sign[d->fIBest][1];
                    break;
                  }
                }
              }
              //-----------------------------------------------------------------------------
              // if the hit radius is greated than some minimal value, set the drift sign
              //-----------------------------------------------------------------------------
              rdrift = hit->driftRadius();

              if ( fabs(rdrift) > _minDriftDoublet){
                hit->setAmbig(doublet->fStrawAmbig[ih]);
              }
              else {
                if (_Final == 0) {
                  hit->setAmbig(0);
                  if (_useMeanResidual) hit->setPenalty(_meanResidual);
                  else                  hit->setPenalty(fabs(rdrift));
                }
                else {
                  //-----------------------------------------------------------------------------
                  // final decision
                  //-----------------------------------------------------------------------------
                  hit->setAmbig(doublet->fStrawAmbig[ih]);
                }
              }
            }
          }
        }
      }
    }
  }

  //-----------------------------------------------------------------------------
  // in the beginning of iteration set external hit errors to a constant
  //-----------------------------------------------------------------------------
  //   void AmbigResolver::initHitErrors(KalRep* krep) const {
  // // get hits and cast to TrkStrawHits
  //     TrkStrawHitVector tshv;
  //     convert(krep->hitVector(),tshv);
  //     for (auto itsh=tshv.begin();itsh!=tshv.end(); ++itsh){
  //       (*itsh)->setTemperature(_tmpErr);
  //     }
  //   }

  //-----------------------------------------------------------------------------
  // in the beginning of iteration set external hit errors to a constant
  // also calculate mean residual over the trajectory
  //-----------------------------------------------------------------------------
  void DoubletAmbigResolver::initHitErrors(KalRep* krep) const {
    double res, sigres, sr2(0);
    int    na(0);

    // get hits and cast to TrkStrawHits
    TrkStrawHitVector tshv;
    convert(krep->hitVector(),tshv);
    // mark uninitialized;

    double penalty = _tmpErr*0.0625*_penaltyScale;

    for (auto itsh=tshv.begin();itsh!=tshv.end(); ++itsh) {
      (*itsh)->setPenalty    (penalty);
      (*itsh)->setTemperature(_tmpErr*_tempScale);
      if ((*itsh)->isActive()) {
        bool hasres = (*itsh)->resid(res, sigres, true);
        if (hasres) {
          sr2 += res*res;
          na += 1;
        }
      }
    }
    //-----------------------------------------------------------------------------
    // 10 is a random number with the meaning of "more than one hit"
    //-----------------------------------------------------------------------------
    DoubletAmbigResolver* dar = (DoubletAmbigResolver*) this;

    if (na > 10) dar->_meanResidual = min(sqrt(sr2/na),_maxMeanResidual);
    else         dar->_meanResidual = _maxMeanResidual;                   // mm
  }

  //-----------------------------------------------------------------------------
  // build list of doublets and resolve the drift signs
  //-----------------------------------------------------------------------------
  bool DoubletAmbigResolver::resolveTrk(KalRep* KRep) const {

    // initialize external hit errors
    initHitErrors(KRep);
    vector<Doublet> listOfDoublets;
    findDoublets (KRep,&listOfDoublets);

    if (listOfDoublets.size() > 0) markMultiplets(KRep, &listOfDoublets);

    return false; // assume no change FIXME!!!
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
*/

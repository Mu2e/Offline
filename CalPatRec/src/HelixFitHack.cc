//
// Object to perform helix fit to straw hits
//
// $Id: HelixFitHack.cc,v 1.6 2014/03/11 15:18:34 murat Exp $
// $Author: murat $ 
// $Date: 2014/03/11 15:18:34 $
//
//
// the following has to come before other BaBar includes
#include "BaBar/BaBar.hh"
#include "CalPatRec/inc/HelixFitHack.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "art/Framework/Services/Optional/TFileService.h"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/DetectorSystem.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
// Root
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TList.h"
#include "TArc.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TFile.h"
// C++
#include <vector>
#include <string>

#include "CalPatRec/inc/THackData.hh"

namespace mu2e 
{
  // statics
  double XYZP::_efac(1.0);
  StrawHitFlag XYZP::_useflag;
  // comparison functor for ordering points
  struct radcomp : public std::binary_function<VALERR, VALERR, bool> {
    bool operator()(VALERR const& r1, VALERR const& r2) { return r1._val < r2._val; }
  };

  // comparison functor for sorting by z
  struct zcomp : public std::binary_function<XYZP,XYZP,bool> {
    bool operator()(XYZP const& p1, XYZP const& p2) { return p1._pos.z() < p2._pos.z(); }
  };

  XYZP::XYZP(size_t index,StrawHit const& sh, StrawHitPosition const& shp,Straw const& straw) :
    _ind(index), _pos(shp.pos()), _phi(shp.pos().phi()), _flag(shp.flag()), _wdir(straw.getDirection()),
    _perr(_efac*shp.posRes(StrawHitPosition::phi)),_rerr(_efac*shp.posRes(StrawHitPosition::rho))
  {
    static const CLHEP::Hep3Vector _zdir(0.0,0.0,1.0);
     _sdir = _zdir.cross(_wdir);
  }

  void
  XYZP::rinfo(CLHEP::Hep3Vector const& center,VALERR& rad) const {
    //    static const double onethird(1.0/3.0);
    //    static const double invsqrt12(1./sqrt(12.0));
    // average the 1-sigma radii to account for non-linear errors
    double rvec = CLHEP::Hep3Vector(_pos - center).perp();
    //    rad._val = onethird*(rvec+rvec1+rvec2);
    rad._val = rvec;
    rad._err = _rerr;
    
  }

  void
  XYZP::finfo(CLHEP::Hep3Vector const& center,VALERR& phi) const {
    //    static const double onethird(1.0/3.0);
    //    static const double invsqrt12(1./sqrt(12.0));
    // average the 1-sigma radii to account for non-linear errors
    double phi0 = CLHEP::Hep3Vector(_pos - center).phi();
    //    rad._val = onethird*(rvec+rvec1+rvec2);
    phi._val = phi0;
    phi._err = _perr; 
  }

  bool XYZP::use() const { 
    return !_flag.hasAnyProperty(_useflag);
  }

  bool XYZP::stereo() const {
    static StrawHitFlag stereo(StrawHitFlag::stereo);
    return _flag.hasAllProperties(stereo);
  }

  void XYZP::setUse(bool use) {
    static StrawHitFlag other(StrawHitFlag::other);
    if(!use)
      _flag.merge(other);
    else
      _flag.clear(other);
  }

  void XYZP::setOutlier(){
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    _flag.merge(outlier);
  }

  bool XYZP::isOutlier(){
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    return _flag.hasAllProperties(outlier);
  }

  //-----------------------------------------------------------------------------
  // HelixFitHackResult
  //-----------------------------------------------------------------------------
  HelixFitHackResult& 
  HelixFitHackResult::operator =(HelixFitHackResult const& other) {
    if(this != &other){
      _hdef = other._hdef;
      _fit = other._fit;
      _center = other._center;
      _radius = other._radius;
      _dfdz = other._dfdz;
      _fz0 = other._fz0;
    }
    return *this;
  }
  
  void
  HelixFitHack::helixParams(HelixFitHackResult const& helix,CLHEP::HepVector& pvec,CLHEP::HepVector& perr) const {
    HelixDef const& mytrk = helix._hdef;
    static const double pi(M_PI);
    //    static const double twopi(2*pi);
    static const double halfpi(pi/2.0);
    // the helix fit introduces a radial bias due to an asymmetry in the detector (more phase space for
    // noise hits outside the circle than inside.  correct for it.
    double radius = helix._radius + _rbias;
    pvec = HepVector(5,0);
    // omega is the inverse transverse radius of the particle's circular motion.  Its
    // signed by the particle angular momentum about the cirle center.
    // This CANNOT be deduced geometrically, so must be supplied as an ad-hoc assumption
    double amsign = copysign(1.0,-mytrk.particle().charge()*bz());
    pvec[HelixTraj::omegaIndex] = amsign/radius;
    // phi0 is the azimuthal angle of the particle velocity vector at the point
    // of closest approach to the origin.  It's sign also depends on the angular
    // momentum.  To translate from the center, we need to reverse coordinates
    pvec[HelixTraj::phi0Index] = atan2(-amsign*helix._center.x(),amsign*helix._center.y());
    // d0 describes the distance to the origin at closest approach.
    // It is signed by the particle angular momentum WRT the origin.
    // The Helix fit radial bias is anti-correlated with d0; correct for it here.
    pvec[HelixTraj::d0Index] = amsign*(helix._center.perp() - helix._radius - 2*_rbias);
    // the dip angle is measured WRT the perpendicular.  It is signed by the particle Z momentum    
    pvec[HelixTraj::tanDipIndex] = amsign/(radius*helix._dfdz);
    // must change conventions here: fz0 is the phi at z=0, z0 is defined at the point of closest approach
    // resolve the loop ambiguity such that the POCA is closest to z=0.
    double dphi = deltaPhi(helix._fz0+amsign*halfpi,pvec[HelixTraj::phi0Index]);
    // choose z0 (which loop) so that f=0 is as close to z=0 as possible
    pvec[HelixTraj::z0Index] = dphi*pvec[HelixTraj::tanDipIndex]/pvec[HelixTraj::omegaIndex];
    // estimated covariance based on average performance.  These should be parameters, FIXME!!!
    perr = HepVector(5,0);
    perr[HelixTraj::d0Index] = 34.0;
    perr[HelixTraj::phi0Index] = 0.02;
    perr[HelixTraj::omegaIndex]  = 0.0002;
    perr[HelixTraj::tanDipIndex] = 0.05;
    perr[HelixTraj::z0Index] = 15.0;

  }

  HelixFitHack::HelixFitHack(fhicl::ParameterSet const& pset) :
    _hDist(0),
    _chi2nFindZ(0.0),
    _eventToLook(-1),
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _mindelta(pset.get<double>("minDelta",5000.0)),
    _minnhit(pset.get<unsigned>("minNHit",10)),
    _lambda0(pset.get<double>("lambda0",1.0)),
    _lstep(pset.get<double>("lstep",0.2)),
    _minlambda(pset.get<double>("minlambda",0.01)),
    _maxniter(pset.get<unsigned>("maxniter",50)),
    _nsigma(pset.get<double>("nsigma",15.0)),
    _nssigma(pset.get<double>("nssigma",3.0)),
    _minzsep(pset.get<double>("minzsep",50.0)),
    _maxzsep(pset.get<double>("maxzsep",500.0)),
    _maxdz(pset.get<double>("maxdz",35.0)),
    _maxdot(pset.get<double>("maxdot",0.9)),
    _maxDfDz(pset.get<double>("maxDfDz",0.01)),
    _rbias(pset.get<double>("radialBias",0.0)),
    _efac(pset.get<double>("ErrorFactor",1.0)),
    _rhomin(pset.get<double>("rhomin",350.0)),
    _rhomax(pset.get<double>("rhomax",780.0)),
    _mindist(pset.get<double>("mindist",50.0)),
    _maxdist(pset.get<double>("maxdist",500.0)),
    _pmin(pset.get<double>("minP",50.0)),
    _pmax(pset.get<double>("maxP",150.0)),
    _tdmin(pset.get<double>("minAbsTanDip",0.3)),
    _tdmax(pset.get<double>("maxAbsTanDip",2.0)),
    _rcmin(pset.get<double>("rcmin",200.0)),
    _rcmax(pset.get<double>("rcmax",350.0)),
    _forcep(pset.get<bool>("forceP",false)),
    _xyweights(pset.get<bool>("xyWeights",false)),
    _zweights(pset.get<bool>("zWeights",false)),
    _filter(pset.get<bool>("filter",true)),
    _plotall(pset.get<bool>("plotall",false)),
    _usetarget(pset.get<bool>("usetarget",true)),
    _bz(0.0)
  {
    XYZP::_efac = _efac;
    std::vector<std::string> bitnames;
    bitnames.push_back("Outlier");
    bitnames.push_back("OtherBackground");
    XYZP::_useflag = StrawHitFlag(bitnames);
    //2014-01-05 gianip[ez added the following liones for initializing the diag-output file
    //    _fOut = new TFile("HelifFitHackDiag.root","recreate");
    //    _hDist =  TH1F("_hDist","Distribution of the point distance in the pattern recognition",1000, 0., 1e6);
    //-----------------------------------//
  }

  HelixFitHack::~HelixFitHack()
  {
    // 2014-01-05 gianipez added a diag level using a TH1F
    //     _fOut->cd();
    if(_hDist){
      delete _hDist;
      _hDist = 0;
    }
    //     _hDist->Write();
    //     _fOut->Close();
    //------------------------------//

  }

  double
  HelixFitHack::bz() const {
    if(_bz == 0.0){
      // find the magnetic field Z component at the origin
      GeomHandle<BFieldManager> bfmgr;
      GeomHandle<DetectorSystem> det;
      // change coordinates to mu2e
      CLHEP::Hep3Vector vpoint(0.0,0.0,0.0);
      CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(vpoint);
      CLHEP::Hep3Vector field = bfmgr->getBField(vpoint_mu2e);
      _bz = field.z();
    }
    return _bz;
  }

  bool
  HelixFitHack::findHelix(HelixFitHackResult& myhel) {
    HelixDef const& mytrk = myhel._hdef;
    //  compute the allowed range in radius for this fit
    double pb = fabs((CLHEP::c_light*1e-3)/(bz()*mytrk.particle().charge()));
    _rmin = _pmin/(pb*sqrt(1.0+_tdmax*_tdmax));
    _rmax = _pmax/(pb*sqrt(1.0+_tdmin*_tdmin));
    //  particle charge, field, and direction affect pitch range
    _dfdzsign = copysign(1.0,-mytrk.particle().charge()*mytrk.fitdir().dzdt()*bz());
    if(_dfdzsign > 0.0){
      _smin = 1.0/(_rmax*_tdmax);
      _smax = 1.0/(_rmin*_tdmin);
    } else {
      _smax = -1.0/(_rmax*_tdmax);
      _smin = -1.0/(_rmin*_tdmin);
    }
    // loop over hits, and store the points
    XYZPVector  xyzp;
    fillXYZP(mytrk,xyzp);
    // 2013-10-18 P.Murat: print 

    int n = xyzp.size();

    XYZP* pt;
    if (_debug > 0) {
      printf("[HelixFitHack::findHelix] x0 = %12.5f y0 = %12.5f r = %12.5f chi2 = %12.5g\n",
	     0.,0.,0.,-1.);
    }

    for(int i=0; i<n; ++i) {
      pt = &xyzp[i];

      if (_debug > 0) {
	printf("[HelixFitHack::findHelix] %08x %2i %12.5f %12.5f %12.5f \n",
	       *((int*) &pt->_flag), pt->use(), pt->_pos.x(), pt->_pos.y(), pt->_pos.z()
	       );
      }
    }

    // call down
    bool retval = findHelix(xyzp,myhel);
    if((retval || _plotall) && _diag>1){
      // fill graphs for display if requested
      plotXY(mytrk,xyzp,myhel);
      plotZ(mytrk,xyzp,myhel);
    }
    return retval;
  }

  bool
  HelixFitHack::findHelix(XYZPVector& xyzp,HelixFitHackResult& myhel) {
    bool retval(false);
    // filter by geometry
    if (_filter) filterDist(xyzp);
    //12 - 10 - 2013 Gianipez added the following line for performing a patern recognition
    THackData* hack;
    hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
    if (hack->ClusterT0() > 0.) {
      doPatternRecognition(xyzp, myhel);
    }
    //--------------------------------------------------//

    if (initCircle_new(xyzp,myhel)) {
      if (_debug > 0) {
	printf("[HelixFitHack::findHelix] xyzp.size() = %i minhit = %i\n",(int) xyzp.size(), _minnhit );
      }
      if(xyzp.size() >= _minnhit){
	// initialize the circle parameters
	// solve for the circle parameters
	//	retval = findXY(xyzp,myhel);
	retval = findXY_new(xyzp,myhel);
	// extend those into the z direction
	if(retval){
	  //2013 - 12 - 25 gianipez added few lines for the ''theoretical'' calculation of dfdz 
	  if (hack->ClusterT0() > 0.) // {
	    // 	    double toll = 50.;
	    // 	    int size(xyzp.size()), k(-9999);
	    // 	    double chi2(0.0), chi2_min(1e10), dx, dx2, dy, dy2;
	    // 	    Hep3Vector p0;
	    // 	    double radius, phi0,tanLambda;
	    // 	    Hep3Vector p1,  p2, p3;
	    // 	    p3 = Hep3Vector(hack->ClusterX(),hack->ClusterY(), hack->ClusterZ());
	    // 	    p2 = xyzp[hack->SeedIndex()]._pos;
	    // 	    for(int i=0; i<size; ++i){
	    // 	      p1 = xyzp[i]._pos;
	    // 	      if( ( xyzp[i].isOutlier() ) ||
	    // 		  !( xyzp[i].stereo() )   || 
	    // 		  ( i==hack->SeedIndex()) ||
	    // 		  ( std::fabs(p1.z() - p2.z()) <toll ) ) goto NEXT_POINT;
	     
	    // 	      calculateHelixParameters(p0, radius, phi0, tanLambda, p1,  p2, p3, myhel);
	    // 	      dx = myhel._center.x() - p0.x();
	    // 	      dx2 = dx*dx;
	    // 	      dy = myhel._center.y() - p0.y();
	    // 	      dy2 = dy*dy;
	    // 	      chi2 = dx2 + dy2;
	    // 	      if(chi2 < chi2_min){
	    // 		chi2_min = chi2;
	    // 		k = i;
	    // 	      }
	    // 	    NEXT_POINT:;
	    // 	    }
	    // 	    if(k>=0){
	    // 	      p1 = xyzp[k]._pos;
	    // 	      calculateHelixParameters(p0, radius, phi0, tanLambda, p1,  p2, p3, myhel);
	      
	    // 	      myhel._fz0  = phi0;
	    // 	      myhel._dfdz = tanLambda/radius;
	    // 	      printf("[HelixFitHack::Optimize_pitch]         POINT USED         \n");
	    // 	      printf("[HelixFitHack::Optimize_pitch] p1 = %5.1f %5.1f %5.1f \n", p1.x(), p1.y(), p1.z());
	    // 	      printf("[HelixFitHack::Optimize_pitch] p2 = %5.1f %5.1f %5.1f \n", p2.x(), p2.y(), p2.z());
	    // 	      printf("[HelixFitHack::Optimize_pitch] p3 = %5.1f %5.1f %5.1f \n", p3.x(), p3.y(), p3.z());
	    // 	      printf("[HelixFitHack::Optimize_pitch] phi0 = %5.3f dfdz = %5.5f\n", phi0, myhel._dfdz);
	    // 	    }
	    	    
	    // 	  }
	    //--------------------------------------------------------------------------------//
	    retval = findZ(xyzp,myhel);

	  if(retval) 
	    retval = findXY_new(xyzp,myhel);
	  
	  // set the success
	  if(retval){
	    myhel._fit = TrkErrCode(TrkErrCode::succeed);
	  } else
	    myhel._fit = TrkErrCode(TrkErrCode::fail,4); // phi-z reconstruction failure
	} else
	  myhel._fit = TrkErrCode(TrkErrCode::fail,3); // xy reconstruction failure
      } else
	myhel._fit = TrkErrCode(TrkErrCode::fail,2); // initialization failure
    } 
    else {
      myhel._fit = TrkErrCode(TrkErrCode::fail,1); // insufficient hits
    }
    return retval;
  }
 
  bool
  HelixFitHack::findXY(XYZPVector& xyzp, HelixFitHackResult& myhel) {
    double rmed, age;
    Hep3Vector center = myhel._center;
    //    findCenterAGE(xyzp,center,rmed,age,false);
    // then, refine that using weights
    bool changed(true);
    unsigned niter(0);
    while(niter < _maxniter && changed){
      findCenterAGE(xyzp,center,rmed,age,_xyweights);
      if(_filter)
	filterXY(xyzp,center,rmed,changed);
      else
	changed = false;
      niter++;
    }
    myhel._center = center;
    myhel._radius = rmed;
    return true;
  }

  //-----------------------------------------------------------------------------
  // this routine does the initial cleanup
  //-----------------------------------------------------------------------------
  bool  HelixFitHack::findXY_new(XYZPVector& xyzp, HelixFitHackResult& Hel) {
    double chi2_min, chi2, x, y;
    int    iworst, np;
    bool   success(0);

    chi2_min = Hel._sxy.chi2DofCircle();

    //    if (chi2_min < 16.) return true;
    if (chi2_min < 16.) {

      if (_debug > 0) {

	printf("[HelixFitHack::findXY_new] don't need to restart, fit is already good!\n");
	printf("[HelixFitHack::findXY_new] x0 = %12.5f y0 = %12.5f r = %12.5f chi2 = %12.5g\n",
	       Hel._sxy.x0(),Hel._sxy.y0(),Hel._sxy.radius(),Hel._sxy.chi2DofCircle());

	np = xyzp.size();
	for (int i=0; i<np; i++) {
	  printf("[HelixFitHack::findXY_new] %08x %2i %12.5f %12.5f %12.5f \n",
		 *((int*) &xyzp[i]._flag), xyzp[i].use(), 
		 xyzp[i]._pos.x(), xyzp[i]._pos.y(), xyzp[i]._pos.z()
		 );
	}
      }
      return true;
    }
    //-----------------------------------------------------------------------------
    // chi2 per degree of freedom is large, perform the cleanup
    // 1. find the worst point
    //-----------------------------------------------------------------------------
    ::LsqSums4 sxy;
    np       = xyzp.size();
    double weight(0.0);
  NEXT_ITERATION:;
    chi2_min = 1e06;
    //chi2_min = Hel._sxy.chi2DofCircle();
    iworst   = -1;
    for (int i=0; i<np; i++) {
      if (xyzp[i].isOutlier()) goto NEXT_POINT;

      sxy.init(Hel._sxy);

      x = xyzp[i]._pos.x();
      y = xyzp[i]._pos.y();
      weight = 0.5;
      if(xyzp[i].stereo()) weight = 1.0;
      sxy.removePoint(x,y, weight);

      chi2 = sxy.chi2DofCircle();

      if (chi2 < chi2_min) {
	iworst   = i;
	chi2_min = chi2;
      }
    NEXT_POINT:;
    }

    if(iworst>=0){
      xyzp[iworst].setOutlier();
      x = xyzp[iworst]._pos.x();
      y = xyzp[iworst]._pos.y();
      weight = 0.5;
      if(xyzp[iworst].stereo()) weight = 1.0;
      Hel._sxy.removePoint(x,y, weight);
    }
    // after removal, chi2DofCircle() = chi2_min
    //    if (chi2_min >= 16.) {
    //    printf("[HelixFitHack::findXY_new] chi2_min = %5.3f\n",chi2_min);
    //    printf("[HelixFitHack::findXY_new] qn = %5.3f\n", Hel._sxy.qn());
    if (chi2_min >= 100.) {
      //-----------------------------------------------------------------------------
      // still bad chi2, repeat the cleanup cycle
      //-----------------------------------------------------------------------------
      if (Hel._sxy.qn() > 10.) {
	goto NEXT_ITERATION;
      }
    }
    else {
      success = true;
      // print circle parameters and tehn - all points with flags:

      if (_debug > 0) {
	printf("[HelixFitHack::findXY_new] x0 = %12.5f y0 = %12.5f r = %12.5f chi2 = %12.5g\n",
	       Hel._sxy.x0(),Hel._sxy.y0(),Hel._sxy.radius(),Hel._sxy.chi2DofCircle());

	for (int i=0; i<np; i++) {
	  printf("[HelixFitHack::findXY_new] %08x %2i %12.5f %12.5f %12.5f \n",
		 *((int*) &xyzp[i]._flag), xyzp[i].use(), 
		 xyzp[i]._pos.x(), xyzp[i]._pos.y(), xyzp[i]._pos.z()
		 );
	}
      }
    }

    Hel._center.set(Hel._sxy.x0(),Hel._sxy.y0(),0.);
    Hel._radius = Hel._sxy.radius();
    if (_debug > 0) {
      printf("[HelixFitHack::findXY_new] retval = %d\n",success ? 1:0);
    }
    return success;
  }


  void
  HelixFitHack::plotXY(HelixDef const& mytrk, XYZPVector const& xyzp,
		       HelixFitHackResult const& myhel) const {
    unsigned igraph = 10*mytrk.eventId()+mytrk.trackId();
    art::ServiceHandle<art::TFileService> tfs;
    char gname[100];
    snprintf(gname,100,"gshxy%i",igraph);
    char bname[100];
    snprintf(bname,100,"bshxy%i",igraph);
    char sname[100];
    snprintf(sname,100,"sshxy%i",igraph);
    char bsname[100];
    snprintf(bsname,100,"bsshxy%i",igraph);
    char title[100];
    snprintf(title,100,"StrawHit XY evt %i trk %i;mm;mm",mytrk.eventId(),mytrk.trackId());
    TH2F* g = tfs->make<TH2F>(gname,title,100,-500,500,100,-500,500);
    TH2F* b = tfs->make<TH2F>(bname,title,100,-500,500,100,-500,500);
    TH2F* s = tfs->make<TH2F>(sname,title,100,-500,500,100,-500,500);
    TH2F* bs = tfs->make<TH2F>(bsname,title,100,-500,500,100,-500,500);
    g->SetMarkerStyle(8);
    g->SetMarkerColor(kGreen);
    b->SetMarkerStyle(4);
    b->SetMarkerColor(kYellow);
    s->SetMarkerStyle(22);
    s->SetMarkerColor(kGreen);
    bs->SetMarkerStyle(26);
    bs->SetMarkerColor(kYellow);
    for(unsigned ihp=0;ihp<xyzp.size();++ihp){
      if(xyzp[ihp].stereo()){
	if(xyzp[ihp].use())
	  s->Fill(xyzp[ihp]._pos.x()-myhel._center.x(),xyzp[ihp]._pos.y()-myhel._center.y());
	else
	  bs->Fill(xyzp[ihp]._pos.x()-myhel._center.x(),xyzp[ihp]._pos.y()-myhel._center.y());
      } else {
	if(xyzp[ihp].use())
	  g->Fill(xyzp[ihp]._pos.x()-myhel._center.x(),xyzp[ihp]._pos.y()-myhel._center.y());
	else
	  b->Fill(xyzp[ihp]._pos.x()-myhel._center.x(),xyzp[ihp]._pos.y()-myhel._center.y());
      }
    }

    TArc* fitarc = new TArc(0.0,0.0,myhel._radius);
    fitarc->SetLineColor(kRed);
    fitarc->SetLineWidth(2);
    fitarc->SetFillStyle(0);
    // draw the detector boundaries
    static double innerrad(380.0);
    static double outerrad(680.0);
    TArc* indet = new TArc(-myhel._center.x(),-myhel._center.y(),innerrad);
    TArc* outdet = new TArc(-myhel._center.x(),-myhel._center.y(),outerrad);
    indet->SetLineColor(kBlue);
    indet->SetFillStyle(0);
    outdet->SetLineColor(kBlue);
    outdet->SetFillStyle(0);
    // add these to the plot
    TList* flist = g->GetListOfFunctions();
    flist->Add(fitarc);
    flist->Add(indet);
    flist->Add(outdet);
  }

  bool
  HelixFitHack::findCenterAGE(XYZPVector const& xyzp,Hep3Vector& center, double& rmed, double& age,bool useweights) {
    // this algorithm follows the method described in J. Math Imagin Vis Dec. 2010 "Robust Fitting of Circle Arcs"
    // initialize step
    double lambda = _lambda0;
    // find median and AGE for the initial center
    findAGE(xyzp,center,rmed,age,useweights);
    // loop while step is large
    unsigned niter(0);
    Hep3Vector descent(1.0,0.0,0.0);
    while(lambda*descent.mag() > _minlambda && niter < _maxniter){
      // fill the sums for computing the descent vector
      SUMS sums;
      fillSums(xyzp,center,rmed,sums,useweights);
      // descent vector cases: if the inner vs outer difference is significant (compared to the median), damp using the median sums,
      // otherwise not.  These expressions take care of the undiferentiable condition on the boundary. 
      double dx(sums._sco-sums._sci);
      double dy(sums._sso-sums._ssi);
      if(fabs(dx) < sums._scc)
        dx += (sums._sco < sums._sci) ? -sums._scc : sums._scc;
      if(fabs(dy) < sums._ssc)
        dy += (sums._sso < sums._ssi) ? -sums._ssc : sums._ssc;
      descent = Hep3Vector(dx,dy,0.0);
      // compute error function, decreasing lambda until this is better than the previous
      double agenew;
      Hep3Vector cnew = center + lambda*descent;
      findAGE(xyzp,cnew,rmed,agenew,useweights);
      // if we've improved, increase the step size and iterate
      if(agenew < age){
        lambda *= (1.0+_lstep);
      } else {
	// if we haven't improved, keep reducing the step till we do
        unsigned miter(0);
        while(agenew > age && miter < _maxniter && lambda*descent.mag() > _minlambda){
          lambda *= (1.0-_lstep);
          cnew = center + lambda*descent;
          findAGE(xyzp,cnew,rmed,agenew,useweights);
          ++miter;
        }
	// if this fails, reverse the descent drection and try again
        if(agenew > age){
	  descent *= -1.0;
	  lambda *= (1.0 +_lstep);
          cnew = center + lambda*descent;
          findAGE(xyzp,cnew,rmed,agenew,useweights);
        }
      }
      // prepare for next iteration
      if(agenew < age){
        center = cnew;
        age = agenew;
      } else {
	static const double minage(0.1);
	if(_debug > 0 && agenew-age>minage)
	  std::cout << "iteration did not improve AGE!!! lambda = " 
		    << lambda  << " age = " << age << " agenew = " << agenew << std::endl;
	break;
      }
      ++niter;
    }
    // check for convergence
    if(_debug > 0 && niter > _maxniter ){
      std::cout << "AGE didn't converge!!! " << std::endl;
      //      return false;
    }
    return true;
  }
  
  bool
  HelixFitHack::findZ(XYZPVector& xyzp,HelixFitHackResult& myhel) {
    // sort points by z
    std::sort(xyzp.begin(),xyzp.end(),zcomp());

    bool success(false);

    THackData* hack;
    int offset(0);
    hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
    if (hack->ClusterT0() > 0.) {
      offset = 1;
      if (_debug > 0) printf("HelixFitHack::FindZ offset = %d\n",offset);
    }

    int N = xyzp.size();
    double z_vec[N+offset], phi_vec[N+offset], weight_vec[N+offset];
    double phi(0.0), z(0.0), weight(0.0);

    CLHEP::Hep3Vector center = myhel._center;
    CLHEP::Hep3Vector pos;

    //calculate the mean pitch assuming conversion electron 
    //hypothesis (so assuming a pitch angle of ~ 0.65)
    double pitch = twopi*std::pow(myhel._radius,2.)*myhel._dfdz;//*tan(0.65);

    if (_debug > 0) printf("[HelixFitHack::findZ] pitch = %5.3f\n", pitch);
 
    //    printf(">>HelixFit.cc:\n");
    //     printf("z-pitch = %5.3f\n",pitch);
    double meanDphi(0.0);

    int realSize(0);
    if (_debug > 0) printf("[HelixFitHack::findZ] xyzp.size() = %u\n", N);
    for(int ixyzp=0; ixyzp < N+offset; ++ixyzp){
      if(ixyzp < N){
	if (xyzp[ixyzp].isOutlier()) goto NEXT_POINT;
	weight = 0.5;
	if(xyzp[ixyzp].stereo()) weight = 1.0;
	pos = xyzp[ixyzp]._pos;
	z   = pos.z();
	z_vec[realSize]      = z;
	weight_vec[realSize] = weight;
      }else {
	weight = 10.0;
	pos = Hep3Vector(hack->ClusterX(), hack->ClusterY(), hack->ClusterZ());
	z   = pos.z();
	z_vec[realSize]      = z;
	weight_vec[realSize] = weight;
      }
      
      phi = CLHEP::Hep3Vector(pos - center).phi();
      //      phi_fashion = phi*180.*2./twopi;
      
      //  if( ixyzp==0){
      // 	printf(">>>HelixFit.cc : findZ \n");
      // 	printf("|      z      |     phi     |    phi_fh    |\n");
      //       }
      
      if(phi<0.0) phi +=twopi;
      
      phi_vec[realSize]= phi;
     
      ++realSize;
    NEXT_POINT:;
      //       printf("| %5.3f  |   %5.3f  |   %5.3f   |\n", z, phi, phi_fashion);
      //       printf("|--------------------------------|\n");
    }
    
    if (_debug > 0) printf("[HelixFitHack::findZ] realsize = %d\n", realSize);
    
    //gianipez procedure for alligning the phi vector
    ::LsqSums4 srphi;
    int i;
    int iworst, jworst, ibest;
    double chi2,chi2min, dfdz;
    double zn, phin, wn;
    int count=0;
    double deltaPhi;
    //    bool changed(false);

    for(int i=realSize-1; i>0; --i){
      count=0;
      ibest = -999;
      
      z   = z_vec[i];
      phi = phi_vec[i];
      weight = weight_vec[i];
      if (hack->ClusterT0() > 0.) {
	z   = hack->ClusterZ();
	phi = phi_vec[realSize-1];
	weight = weight_vec[realSize-1];
      }
      if(i==realSize-1) {
	srphi.addPoint(z,phi,weight);
	dfdz = myhel._dfdz;
      }
      else
	dfdz = srphi.dfdz();

     
      deltaPhi = phi - phi_vec[i-1];
      meanDphi = (z - z_vec[i-1])*dfdz;

      while( (std::fabs(deltaPhi - meanDphi) > twopi ||
	      std::fabs(deltaPhi + twopi - meanDphi) < std::fabs(deltaPhi - meanDphi) )&&
	     count < 3){
	phi_vec[i-1] -= twopi;
	deltaPhi = phi - phi_vec[i-1];
	++count;
      }

      srphi.addPoint(z_vec[i-1], phi_vec[i-1],weight_vec[i-1]);

      if( i < realSize-1 ){
	chi2min = srphi.chi2rphiDofCircle();
	phin = phi_vec[i-1];
	zn   = z_vec[i-1];
	wn   = weight_vec[i-1];
	for(int j=0; j<2; ++j){
	  srphi.removePoint(zn, phin, wn);
	  if(j==0) phin += twopi;
	  if(j==1) phin = phi_vec[i-1] - twopi;
	  srphi.addPoint(zn, phin, wn);
	  chi2 = srphi.chi2rphiDofCircle();
	  if(chi2 < chi2min){
	    ibest = j;
	    chi2min = chi2;
	  }
	}
	srphi.removePoint(zn, phin, wn);
	if(ibest == 0) phi_vec[i-1] +=twopi;
	if(ibest == 1) phi_vec[i-1] -=twopi;
	phin = phi_vec[i-1];
	srphi.addPoint(zn, phin, wn);
      }
    }

    //--------------------------------------------------//
    for(int i=0; i<realSize; ++i){
      z = z_vec[i];
      phi = phi_vec[i];
      if(i!=realSize-1) 
	myhel._srphi.addPoint(z,phi);
      else
	myhel._srphi.addPoint(z,phi,10.0);
    
      if (_debug > 0) printf("[HelixFitHack::findZ] %5.3f     %5.3f   \n", z, phi);
    }
//-----------------------------------------------------------------------------
// perform a cleanup
//-----------------------------------------------------------------------------
  NEXT_ITERATION:;
    i=0;
    iworst = -1;
    jworst = -1;
    chi2min = 1e10;
    //    chi2min = myhel._srphi.chi2rphiDofCircle();
    for(int ixyzp=0; ixyzp < N; ++ixyzp){
      if (xyzp[ixyzp].isOutlier()) goto NEXT_P;
      srphi.init(myhel._srphi);
      pos = xyzp[ixyzp]._pos;
      z   = z_vec[i];
      phi = phi_vec[i];
      weight = 0.5;
      if(xyzp[ixyzp].stereo()) weight = 1.0;
      srphi.removePoint(z, phi, weight);
      chi2 = srphi.chi2rphiDofCircle();
      //printf("[HelixFitHack::findZ] chi2 = %5.3e chi2min = %5.3e\n", chi2, chi2min);
      if (chi2 < chi2min) {
	iworst   = ixyzp;
	chi2min = chi2;
	jworst = i;
      }
      ++i;
    NEXT_P:;
    }
    
    if(iworst>=0 && myhel._srphi.qn() > 3.){
      xyzp[iworst].setOutlier();
      z   = z_vec[jworst];
      phi = phi_vec[jworst];
      --realSize;
      weight = 0.5;
      if(xyzp[iworst].stereo()) weight = 1.0;
      myhel._srphi.removePoint(z, phi, weight);
      if (_debug > 0) printf("[HelixFitHack::findZ_removed] %5.3f     %5.3f   \n", z, phi);
    }
    if(myhel._srphi.qn()<=3) chi2min = myhel._srphi.chi2rphiDofCircle();
    
    // printf("[HelixFitHack::findZ_clean] chi2_min = %5.3f\n",chi2min);
    if (chi2min >= 1.) {
      
      if (myhel._srphi.qn() > 10.) {
	goto NEXT_ITERATION;
      }
    }else
      success = true;
    //----------------------------------------------------------------------//
  
    myhel._fz0  = myhel._srphi.phi0();
    myhel._dfdz = myhel._srphi.dfdz();

    
    if (_debug > 0) {
      printf("[HelixFitHack::findZ] phi_0 = %5.3f dfdz = %5.5f chi2N = %5.3f\n", 
	     myhel._fz0,myhel._dfdz, myhel._srphi.chi2rphiDofCircle() );
    }
    
    _chi2nFindZ = myhel._srphi.chi2rphiDofCircle();
    if(_chi2nFindZ < 0.0)  _eventToLook = myhel._hdef.eventId();
    // if( realSize> (int) _minnhit) {
    //       printf("[HelixFitHack::findZ] retval = 1\n" );
    //       return true;
    //     }else {
    //       printf("[HelixFitHack::findZ] retval = 0\n" );
    //       return false;
    //     }
    if (_debug > 0) printf("[HelixFitHack::findZ] retval = %d\n",success ? 1:0);
    return success;
    
  }

  void
  HelixFitHack::plotZ(HelixDef const& mytrk, XYZPVector const& xyzp, HelixFitHackResult const& myhel) const {
    unsigned igraph = 10*mytrk.eventId()+mytrk.trackId();
    art::ServiceHandle<art::TFileService> tfs;
    char gname[100];
    snprintf(gname,100,"gshphiz%i",igraph);
    char bname[100];
    snprintf(bname,100,"bshphiz%i",igraph);
    char sname[100];
    snprintf(sname,100,"sshphiz%i",igraph);
    char bsname[100];
    snprintf(bsname,100,"bsshphiz%i",igraph);
    char title[100];
    snprintf(title,100,"StrawHit #phi Z evt %i trk %i;mm;rad",mytrk.eventId(),mytrk.trackId());
    TH2F* g = tfs->make<TH2F>(gname,title,100,-1500,1500,100,-12.5,12.5);
    TH2F* b = tfs->make<TH2F>(bname,title,100,-1500,1500,100,-12.5,12.5);
    TH2F* s = tfs->make<TH2F>(sname,title,100,-1500,1500,100,-12.5,12.5);
    TH2F* bs = tfs->make<TH2F>(bsname,title,100,-1500,1500,100,-12.5,12.5);
    g->SetMarkerStyle(8);
    g->SetMarkerColor(kGreen);
    b->SetMarkerStyle(4);
    b->SetMarkerColor(kYellow);
    s->SetMarkerStyle(22);
    s->SetMarkerColor(kGreen);
    bs->SetMarkerStyle(26);
    bs->SetMarkerColor(kYellow);
    for(unsigned ih=0;ih<xyzp.size();++ih){
      if(xyzp[ih].stereo()){
	if(xyzp[ih].use())
	  s->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	else
	  bs->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
      } else {
	if(xyzp[ih].use())
	  g->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	else
	  b->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
      }
    }
    TF1* line = new TF1("line","[0]+[1]*x",-1500,1500);
    line->SetParameter(0,myhel._fz0);
    line->SetParameter(1,myhel._dfdz);
    line->SetLineColor(kRed);
    TList* flist = g->GetListOfFunctions();
    flist->Add(line);
  }

  bool
  HelixFitHack::initCircle(XYZPVector const& xyzp,HelixFitHackResult& myhel) {
    bool retval(false);
    static const double mind2 = _mindist*_mindist;
    using namespace boost::accumulators;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accx, accy, accr;
    // form all triples, and compute the circle center for unaligned hits.  I can aford to be choosy
    unsigned ntriple(0);
    unsigned nxyzp = xyzp.size();
    std::vector<CLHEP::Hep3Vector> pos;
    pos.reserve(xyzp.size());
    for(size_t ixyzp=0; ixyzp < nxyzp; ++ixyzp) {
      if(xyzp[ixyzp].use()){
	pos.push_back(xyzp[ixyzp]._pos);
      }
    }
    // I should randomly pick entries if there are too many, FIXME
    // add the target (if requested)
    static const CLHEP::Hep3Vector tpos(0.0,0.0,0.0);
    if(_usetarget)
      pos.push_back(tpos);

    size_t np = pos.size();
    for(size_t ip=0;ip<np;++ip){
      // pre-compute some values
      double ri2 = pow(pos[ip].x(),2) + pow(pos[ip].y(),2);
      for(size_t jp=ip+1;jp<np; ++jp){
	if(pos[ip].perpPart().diff2(pos[jp].perpPart()) > mind2){
	  double rj2 = pow(pos[jp].x(),2) + pow(pos[jp].y(),2);
	  for(size_t kp=jp+1;kp<np; ++kp){
	    if(pos[ip].perpPart().diff2(pos[kp].perpPart()) > mind2 &&
	       pos[jp].perpPart().diff2(pos[kp].perpPart()) > mind2){
	      // this effectively measures the slope difference
	      double delta = (pos[kp].x() - pos[jp].x())*(pos[jp].y() - pos[ip].y()) - 
		(pos[jp].x() - pos[ip].x())*(pos[kp].y() - pos[jp].y());
	      if(fabs(delta) > _mindelta){
		double rk2 = pow(pos[kp].x(),2) + pow(pos[kp].y(),2);
		// find circle center for this triple
		double cx = 0.5* (
				  (pos[kp].y() - pos[jp].y())*ri2 + 
				  (pos[ip].y() - pos[kp].y())*rj2 + 
				  (pos[jp].y() - pos[ip].y())*rk2 ) / delta;
		double cy = -0.5* (
				   (pos[kp].x() - pos[jp].x())*ri2 + 
				   (pos[ip].x() - pos[kp].x())*rj2 + 
				   (pos[jp].x() - pos[ip].x())*rk2 ) / delta;
		double rho = sqrt(pow(pos[ip].x()-cx,2)+pow(pos[ip].y()-cy,2));
		if(rho > _rcmin && rho<_rcmax){
		  // accumulate 
		  ++ntriple;
		  accx(cx);
		  accy(cy);
		  accr(rho);
		}
	      }
	    }
	  }
	}
      }
    }
    if(ntriple > _minnhit){
      double centx = extract_result<tag::median>(accx);
      double centy = extract_result<tag::median>(accy);
      double rho = extract_result<tag::median>(accr);
      myhel._center = CLHEP::Hep3Vector(centx,centy,0.0);
      myhel._radius = rho;
      if(_forcep){
	myhel._radius = std::max(std::min(myhel._radius,_rmax),_rmin);
	retval = true;
      } else {
	retval = myhel._radius >= _rmin && myhel._radius <= _rmax;
      }
    }
    return retval;
  }

  //-----------------------------------------------------------------------------
  // substitute for initCircle - know that there is a cluster
  //-----------------------------------------------------------------------------
  bool HelixFitHack::initCircle_new(XYZPVector const& xyzp, HelixFitHackResult& Hel) {
    bool   retval(false);
    //    static const double mind2 = _mindist*_mindist;
    double   x, y, x_0, y_0, r;
    XYZP*    p;
    // first point - center, with low weight

    //    Hel._sxy.addPoint(0.,0.,0.5);
    Hel._sxy.addPoint(0.,0.,0.05);

    //second point: the EMC cluster position!
    THackData* hack;
    
    hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
    if (hack->ClusterT0() > 0.) {
      Hel._sxy.addPoint(hack->ClusterX(), hack->ClusterY(), 10. );
      if (_debug > 0) {
	printf("[HelixFitHack::initCircle_new] x_EMC = %12.5f y_EMC = %12.5f\n",
	       hack->ClusterX(), hack->ClusterY());
      }
    }

    int   n = xyzp.size();
    //2013 - 12 - 25 gianipez added a paramter for weghting the difference between stereo and straw hits
    double weight(0.);//use the weight parameter for distinguisching between stereohits and strawhits

    // don't use the outliers
    for (int i=0; i<n; i++) {
      p = (XYZP*) &xyzp[i];
      if (! p->isOutlier()) {
	x        = p->_pos.x();
	y        = p->_pos.y();
	weight = 0.5;
	if(p->stereo()) weight = 1.0;
	Hel._sxy.addPoint(x,y, weight);
      }
    }

    x_0  = Hel._sxy.x0();
    y_0  = Hel._sxy.y0();
    r    = Hel._sxy.radius();

    if (_debug > 0) {
      printf("[HelixFitHack::initCircle_new] x0 = %12.5f y0 = %12.5f r = %12.5f chi2 = %12.5g\n",
	     x_0,y_0,r,Hel._sxy.chi2DofCircle());

      for (int i=0; i<n; i++) {
	p = (XYZP*) &xyzp[i];
	printf("[HelixFitHack::initCircle_new] %08x %2i %12.5f %12.5f %12.5f \n",
	       *((int*) &p->_flag), p->use(), p->_pos.x(), p->_pos.y(), p->_pos.z()
	       );
      }
    }

    Hel._center.set(x_0,y_0,0.0);
    Hel._radius = r;
    
    retval = (r >= _rmin) && (r <= _rmax);
    if (_debug > 0) printf("[HelixFitHack::initCircle_new] retval = %d\n",retval ? 1:0);
    return retval;
  }

  void HelixFitHack::fillXYZP(HelixDef const& mytrk, XYZPVector& xyzp) {
    //12-09-2013 gianipez modified this procedure to avoid the doubling of the 
    // same stereohitposition 
    //-------------------------------------------------------------------------
    const Tracker& tracker = getTrackerOrThrow();

    
    const std::vector<hitIndex> shIndices = mytrk.strawHitIndices();
    int size = shIndices.size();
    int stIndeces[10000];
    //initialize the vector
    for(int i=0; i<size;++i){
      stIndeces[i] = -9999;
    }
    //--------------------------------------------------------------------------------
    
    if(mytrk.strawHitPositionCollection() != 0){
      for(int i=0; i<size; ++i){
	StrawHit const& sh = mytrk.strawHitCollection()->at(shIndices[i]._index);
	Straw const& straw= tracker.getStraw(sh.strawIndex());
	StrawHitPosition const& shp = mytrk.strawHitPositionCollection()->at(shIndices[i]._index);
	
	XYZP pos(shIndices[i]._index,sh,shp,straw);
	if (_debug > 0) {
	  if(i == 0){
	    printf("----->HelixFitHack::fillXYXP\n");
	    printf("|   strawIndex    |  stereo index   |          pos          |\n");
	  }					
	  printf("| %10.3d    | %10.3d    | %5.3f, %.3f, %5.3f |\n", 
		 (int) shIndices[i]._index,
		 shp.stereoHitIndex(),
		 shp.pos().x(), shp.pos().y(), shp.pos().z());
	}

	//if the stereo index is <0, it means that no stereo hit was created
	if(shp.stereoHitIndex() < 0 ){
	  xyzp.push_back(pos);
	}else {
	  bool found(false);
	  int k=0, t=0;
	  //printf("search if the stereo index ia lerady present...\n");
	  while(!found && k< size){
	    //  printf("stIndeces[%d] = %5.3d ? %u\n",k,stIndeces[k], shp.stereoHitIndex());
	    
	    if( stIndeces[k] == shp.stereoHitIndex()) found = true;
	    if(stIndeces[k]>0) ++t;
	    ++k;
	  }
	  // printf("indexes: k = %d, t = %d\n", k, t);
	  if(!found) {
	    xyzp.push_back(pos);
	    stIndeces[t] = shp.stereoHitIndex();
	  }
	  
	}
      }//end loop onm the strawindexes
    } else {
      static const double twoinvsqrt12(2.0/sqrt(12.0));
      ConditionsHandle<TrackerCalibrations> tcal("ignored");
      for(int j=0; j<size; ++j){
	StrawHit const& sh = mytrk.strawHitCollection()->at(shIndices[j]._index);
	Straw const& straw= tracker.getStraw(sh.strawIndex());
	SHInfo shinfo;
	tcal->StrawHitInfo(straw,sh,shinfo);
	xyzp.push_back(XYZP(shIndices[j]._index,shinfo._pos,straw.getDirection(),shinfo._tdres,twoinvsqrt12*straw.getRadius()));
      }
    }
  }

  void
  HelixFitHack::findAGE(XYZPVector const& xyzp, Hep3Vector const& center,double& rmed, double& age,bool useweights) {
    using namespace boost::accumulators;
    // protection against empty data
    if(xyzp.size() == 0)return;
    // fill radial information for all points, given this center
    std::vector<VALERR> radii;
    unsigned nxyzp = xyzp.size();
    double wtot(0.0);
    for(unsigned ixyzp=0; ixyzp < nxyzp; ++ixyzp){
      if(xyzp[ixyzp].use()){
	// find radial information for this point
	VALERR rad;
	xyzp[ixyzp].rinfo(center,rad);
	radii.push_back(rad);
	// compute the normalization too
	wtot += useweights ? 1.0/rad._err : 1.0;
      }
    }
    // find the weighted median radius
    accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > accr;
    for(unsigned irad=0;irad<radii.size();++irad){
      double wt = useweights ? 1.0/radii[irad]._err : 1.0;
      accr(radii[irad]._val, weight = wt); 
    }
    rmed = extract_result<tag::weighted_median>(accr);
    // if requested force radius into range
    if(_forcep)
      rmed = std::max(std::min(rmed,_rmax),_rmin);
    // now compute the AGE
    age = 0.0;
    for(unsigned irad=0;irad<radii.size();++irad){
      double wt = useweights ? 1.0/radii[irad]._err : 1.0;
      age += wt*fabs(radii[irad]._val-rmed);
    }
    // normalize
    age *= radii.size()/wtot;
  }

  void
  HelixFitHack::fillSums(XYZPVector const& xyzp, Hep3Vector const& center,double rmed,SUMS& sums,bool useweights) {
    // initialize sums
    sums.clear();
    // compute the transverse sums
    double wtot(0.0);
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp].use()){
	// find radial information for this point
	VALERR rad;
	xyzp[ixyzp].rinfo(center,rad);
	double wt = useweights ? 1.0/rad._err : 1.0;
	wtot += wt;
	// now x,y projections
	double pcos = (xyzp[ixyzp]._pos.x()-center.x())/rad._val;
	double psin = (xyzp[ixyzp]._pos.y()-center.y())/rad._val;
	// 3 conditions: either the radius is inside the median, outside the median, or 'on' the median.  We define 'on'
	// in terms of the error
	if(fabs(rad._val -rmed) < rad._err ){
	  sums._scc += wt*fabs(pcos);
	  sums._ssc += wt*fabs(psin);
	  ++sums._nc;
	} else if (rad._val > rmed) {
	  sums._sco += wt*pcos;
	  sums._sso += wt*psin;
	  ++sums._no;
	} else {
	  sums._sci += wt*pcos;
	  sums._ssi += wt*psin;        
	  ++sums._ni;
	}
      }  
    }
    // normalize to unit weight
    unsigned nused = sums._nc + sums._no + sums._ni;
    sums._scc *= nused/wtot;
    sums._ssc *= nused/wtot;
    sums._sco *= nused/wtot;
    sums._sso *= nused/wtot;
    sums._sci *= nused/wtot;
    sums._ssi *= nused/wtot;
  }

  void
  HelixFitHack::filterDist(XYZPVector& xyzp) {
    using namespace boost::accumulators;
    static const double pi(M_PI);
    static const double twopi(2*pi);
    // first, resolve phi.   Use the average X and Y to define the initial
    // phi value, to avoid looping issues
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accx;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accy;
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp].use()){
	accx(xyzp[ixyzp]._pos.x());
	accy(xyzp[ixyzp]._pos.y());
      }
    }
    double mx = extract_result<tag::median>(accx);
    double my = extract_result<tag::median>(accy);
    double mphi = atan2(my,mx);

    /*
      for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp].use()){
      double dphi = xyzp[ixyzp]._phi - mphi;
      if(fabs(dphi) > pi){
      if(dphi > 0)
      xyzp[ixyzp]._phi -= twopi;
      else
      xyzp[ixyzp]._phi += twopi;
      }
      }
      }
      // now cut
      CLHEP::Hep3Vector mh(mx,my,0.0);
      for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      double dist = sqrt(xyzp[ixyzp]._pos.perpPart().diff2(mh));
      if(dist > _maxdist)xyzp[ixyzp].setOutlier();
      }
    */
    // 2013-10-17 G.Pezzullo, P.Murat : try to change the cleanup logic
    THackData* hack;

    hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");

    if (hack->ClusterT0() > 0.) {
      mphi = atan2(hack->fData[1],hack->fData[0]);
    }
    
    if (_debug > 0) {
      printf("[HelixFitHack::filterDist] x0 = %12.5f y0 = %12.5f r = %12.5f chi2 = %12.5g\n",
	     -1., -1., -1., -1.);
    }

    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp].use()){
	double dphi = xyzp[ixyzp]._phi - mphi;

	if (dphi >  pi) dphi -= twopi;
	if (dphi < -pi) dphi += twopi;

	if(fabs(dphi) > pi/2){
	  xyzp[ixyzp].setOutlier();
	}
      }
      if (_debug > 0) {
	printf("[HelixFit::filterDist]: %08x %3i %12.5f %12.5f %12.5f %12.5f \n",
	       *((int*) &xyzp[ixyzp]._flag),
	       ixyzp,
	       xyzp[ixyzp]._pos.x(), xyzp[ixyzp]._pos.y(),xyzp[ixyzp]._pos.z(),
	       xyzp[ixyzp]._phi);
      }
    }
  }

  void
  HelixFitHack::filterXY(XYZPVector& xyzp, Hep3Vector const& center,double rmed,bool& changed) {
    static StrawHitFlag other(StrawHitFlag::other);
    changed = false;
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      VALERR rad;
      xyzp[ixyzp].rinfo(center,rad);
      bool olduse = xyzp[ixyzp].use();
      // update the flag each iteration, until we've converged
      if(fabs(rad._val -rmed) > _nsigma*rad._err)
	xyzp[ixyzp]._flag.merge(other);
      else
	xyzp[ixyzp]._flag.clear(other);
      changed |= olduse != xyzp[ixyzp].use();
    }
  }

  double HelixFitHack::deltaPhi(double phi1, double phi2){
    static const double pi(M_PI);
    static const double twopi(2*pi);
    double dphi = fmod(phi2-phi1,twopi);
    if(dphi>pi)dphi -= twopi;
    if(dphi<-pi)dphi += twopi;
    return dphi;
  }

  void
  HelixFitHack::doPatternRecognition(XYZPVector& xyzp, HelixFitHackResult& mytrk ){
    int np = xyzp.size();
    int seedIndex(-1);//, targetIndex(-1);
    double chi2, chi2_min;
    int countGoodPoints(0), countGoodPoints_max(0);
    chi2_min = 1e10;
    // int nTargets(16);

    // for(int j=0; j<nTargets; ++j){
    for (int i=0; i<np; i++) {
      if (xyzp[i].isOutlier()) goto NEXT_POINT;
      findTrack(xyzp, i, chi2, countGoodPoints, mytrk);
      if(chi2 < 1e10) {
	if (_debug > 0) {
	  printf("[HelixFitHack:doPatternRecognition] index = %d x = %5.3f y = %5.3f z = %5.3f chi2 = %5.3e\n", i,
		 xyzp[i]._pos.x(), 
		 xyzp[i]._pos.y(),
		 xyzp[i]._pos.z(),
		 chi2);
	  printf("//------------------------------------------------------------------//\n");
	}
      }
      if(chi2 < 1e10) {
	//if(chi2 < chi2_min) {
	if((countGoodPoints>countGoodPoints_max) ||
	   ((countGoodPoints == countGoodPoints_max) && (chi2 < chi2_min))){
	  countGoodPoints_max = countGoodPoints;
	  seedIndex   = i;  
	  //targetIndex = j;
	  chi2_min    = chi2;
	}
	
      }
      
    NEXT_POINT:;
    }
    //}

    if (_debug > 0) {
      printf("[HelixFitHack:doPatternRecognition] candidate is index = %d number of good points = %d\n", 
	     seedIndex, countGoodPoints_max);
    }
    
    if(seedIndex>=0){
      THackData* hack;
      hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
      hack->fData[2] = seedIndex;
      bool cleanPattern(true);
      findTrack(xyzp, seedIndex, chi2, countGoodPoints, mytrk, cleanPattern);
    }
    
  }

  
  void
  HelixFitHack::findTrack(XYZPVector& xyzp, int seedIndex, 
			  double&chi2, int &countGoodPoints_end,
			  HelixFitHackResult& mytrk,
			  bool cleanPattern){
    THackData* hack;
    countGoodPoints_end = 0;
    hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
    //p0 is the position of the center of the helix
    //p1 is the center of the capture target
    //p2 is the point in the position seedIndex on the vector xyzp
    //p3 is the postion of the EMC cluster

    Hep3Vector p0, p1, p2, p3; 
    p1 = Hep3Vector(0., 0., 0.);//5971. - 10200.
    
    // double z0 = 5471. + double(targetIndex)*targetStep - trackerOffSet;
    //  p1.setZ(z0);
    p2 = xyzp[seedIndex]._pos;
    p3 = Hep3Vector(hack->ClusterX(),hack->ClusterY(), hack->ClusterZ());
    // p1.setZ(p2.z());
    
    double radius,phi0,tanLambda;
    //calculate now the helix paramters using p1,p2,p3
    calculateTrackParameters(p0,radius,phi0,tanLambda,
			     p1,p2,p3,
			     mytrk,
			     cleanPattern);
    //----------------------------------------------------------------------//
    double dx,dy,phi,dx2, dy2;
    int np = xyzp.size();
 
    //index for caunting the number of outlayersfound
    int k(0);
    //--------------------------------------------------//

    //index for setting the sostitute-point for the target center
    int rejectList[np];
    int  countGoodPoints(0);
    for(int i=0; i<np; ++i){
      rejectList[i] = -9999;
    }
    //------------------------------------------------------------//
    double distGoodPoint(600.);//(1e5);
    double dist(0.0);
    double tollMin(100.), tollMax(500.);
    //2014-03-10 gianipez changed the values of the following tollerance
    // for X and Y distance to avoid delta electron for corrupting the pattern-reco
    //    double tollXmin(50.), tollYmin(50.);
    double tollXmin(100.), tollYmin(100.);

    int goodPoint(-999);
    //2014-01-29 gianipez added the followign line
    int goodPoint_end(-999);
    //--------------------------------------------------//

    double weight(0.0);
    double deltaZ(0.), deltaX(0.), deltaY(0.);
    bool removeTarget(true);
    bool isStored(false);

    //parameter used for calculating ''manually'' the value of dfdz
    double z0,z1,phi_0,phi_1;
    double dfdz = tanLambda/radius;

    int nmodes=20;
    double chi2min(1e10);
    double dfdz_end, phi0_end, radius_end;
    double x0_end, y0_end;
    
    int mode;
    if(dfdz<0.){
      if (_debug > 0) {
	printf("point x = %5.3f y = %5.3f z = %5.3f returns begative dfdz\n",
	       p2.x(), p2.y(), p2.z());
      }
                                                            goto END;
    }
    //----------------------------------------------------------------------//
    chi2 = 0.0;

    //2014-01-05 gianipez added histo for diagnostic
    if(cleanPattern){
      //     art::ServiceHandle<art::TFileService> tfs;
      int igraph = 10*int(mytrk._hdef.eventId()+mytrk._hdef.trackId());
      char gname[100];
      snprintf(gname,100,"fHdist_%i",igraph);
      char title[100];
      snprintf(title,100,"evt %i trk %i",mytrk._hdef.eventId(),
 	       mytrk._hdef.trackId());
      _hDist = new TH1F(gname,title, 1000, 0., 1e6);//tfs->make<TH1F>(gname,title, 1000, 0., 1e6);
      //hack->fHdist = new TH1F(gname,title, 1000, 0., 1e6);tfs->make<TH1F>(gname,title, 1000, 0., 1e6);
    }
    //----------------------------------------//
     
    for(int j=1;j<=nmodes; ++j){
      //reset some variables used in the loop
      for(int i=0; i<np; ++i){
	rejectList[i] = -9999;
      }
      removeTarget = true;
      countGoodPoints = 0;
      goodPoint = -9999;
      chi2 = 0;
      p1 = Hep3Vector(0., 0., 0.);
      calculateTrackParameters(p0,radius,phi0,tanLambda,
			       p1,p2,p3,
			       mytrk,
			       false);
      dfdz = tanLambda/radius;

      hack->fData[3] = -9999;
      //--------------------------------------------------//

    
      //2014-03-10 Gianipez and Pasha set a limit on th dfdz value
      // still need to be optimized
      if(std::fabs(dfdz*double(j)) > _maxDfDz) break; //0.01) break;
      //now set the mode of dfdz to study
      dfdz = dfdz*double(j);

      for (int i=0; i<np; i++) {
	if (xyzp[i].isOutlier()) goto NEXT_POINT;
	weight = 0.5;
	if(xyzp[i].stereo()) weight = 1.0;
	//      phi = phi0 + (xyzp[i]._pos.x() - p0.z())*tanLambda;
	deltaZ = xyzp[i]._pos.z() - p2.z();
	phi = phi0 + (deltaZ)*dfdz;//tanLambda/radius;

	dx = p0.x() + radius*std::cos(phi);
	dx -=  xyzp[i]._pos.x();
	dx2 = dx*dx;
      
	dy = p0.y() + radius*std::sin(phi);
	dy -=  xyzp[i]._pos.y();
	dy2 = dy*dy;
	dist = dx2 + dy2;
	
	isStored = false;

	if( dist <= distGoodPoint ){
	  ++countGoodPoints;
	  for(int j=0; j<np; ++j){
	    if(rejectList[j] == i)
	      isStored = true;
	  }
	  
	  if(!isStored){
	    deltaY = std::fabs(p2.y() - xyzp[i]._pos.y());
	    deltaX = std::fabs(p2.x() - xyzp[i]._pos.x());
	    if( ( deltaZ > tollMin ) &&
		( deltaZ < tollMax ) &&
		( deltaX > tollXmin) &&
		( deltaY > tollYmin) ){
	      goodPoint = i;
	    }
	  }
	}

	if(dist < distGoodPoint){
	  chi2 += weight*(dist);
	}
      
	if(countGoodPoints == 2 &&
	   removeTarget &&
	   (goodPoint >=0) &&
	   (goodPoint != hack->fData[3]) ){
	  removeTarget = false;
	  countGoodPoints = 0;
	  i=-1;
	  
	  //	  hack->fData[3] = goodPoint;
	  p1 = xyzp[goodPoint]._pos;
	  calculateTrackParameters(p0,radius,phi0,tanLambda,
				   p1,p2,p3,
				   mytrk,
				   cleanPattern);
	  z0 = p2.z();
	  z1 = p1.z();
	  phi_0 = CLHEP::Hep3Vector(p2 - p0).phi();
	  phi_1 = CLHEP::Hep3Vector(p1 - p0).phi();
	  int hackpoint = hack->fData[3];
	  if (_debug > 0) {
	    printf("[HelixFitHack:addGoodPoint] goodPoint = %5d hackPoint = %5d\n",
		   goodPoint,
		   hackpoint);
	    printf("[HelixFitHack:addGoodPoint] old dfdz = %5.5f mode = %i \n",
		   dfdz, j);
	  }
	  calculateDfDz(phi_0, phi_1, z0, z1, dfdz);

	  if (_debug > 0) {
	    printf("[HelixFitHack:addGoodPoint] dfdz = %5.5f \n",dfdz);
	  }
	
 	  //what to do if dfdz s negative?
 	  if(dfdz<0.0) {
 	    for(int j=0; j<np; ++j){
 	      if(rejectList[j]<0){
 		rejectList[j] = goodPoint;
 		break;
 	      }
 	    }
	    if (_debug > 0) printf("[HelixFitHack:addGoodPoint] FAILED dfdz<0, resetting...\n");
 	    p1 = Hep3Vector(0., 0., 0.);
 	    calculateTrackParameters(p0,radius,phi0,tanLambda,
 				     p1,p2,p3,
 				     mytrk,
 				     false);
 	    dfdz = (tanLambda/radius)*double(j);
 	    removeTarget = true;
 	    countGoodPoints = 0;
 	  }
	  //------------------------------------------------------------//
	}
      NEXT_POINT:;
      }
      if(countGoodPoints>2){
	chi2 /= double(countGoodPoints);

	mode = j;
	if (_debug > 0) {
	  printf("[HelixFitHack:searchGoodPoint] chi2 = %5.3f nGoodPoints = %d dfdz = %5.5f mode = %i\n", 
		 chi2,
		 countGoodPoints,
		 dfdz,
		 mode);
	}
	//	if(chi2<chi2min){
	if( (countGoodPoints > countGoodPoints_end) ||
	    ((countGoodPoints==countGoodPoints_end) && (chi2 < chi2min)) ){
	  countGoodPoints_end = countGoodPoints;
	  //2014-01-29 gianipez added the following line
	  if(goodPoint>=0) goodPoint_end = goodPoint; //hack->fData[3] = goodPoint;
	  //--------------------------------------------------//
	  chi2min = chi2;
	  x0_end = p0.x();
	  y0_end = p0.y();
	  dfdz_end = dfdz;
	  phi0_end = phi0;
	  radius_end = radius;
	  if (_debug > 0) {
	    printf("[HelixFitHack:setGoodPoint] chi2min = %5.3f hackpoint = %d dfdz = %5.5f\n",   chi2,
		   goodPoint,
		   dfdz_end );
	  }
	}
      }


    }//end loop on nmodes
    
    if(countGoodPoints_end>2){
      //2014-01-29 gianipez added the following line
      hack->fData[3] = goodPoint_end;
      //--------------------------------------------------//
      chi2 = chi2min;
      if (_debug > 0) {
	printf("[HelixFitHack:EndFindTrack] chi2 = %5.3f nGoodPoints = %d dfdz = %5.5f mode = %i\n", 
	       chi2min,
	       countGoodPoints_end,
	       dfdz_end,
	       mode);
      }
      
    }else{
      chi2 = 1e10;
    }
    
    if(cleanPattern){
      if(hack->fData[3]>=0) 
	p1 = xyzp[hack->fData[3]]._pos;
      else
	p1 = Hep3Vector(0., 0., 0.);
      mytrk._center.set(x0_end, y0_end, 0.0);
      mytrk._radius = radius_end;
      //now calculate the phi coordinate at z = 0
      mytrk._fz0    = phi0_end + dfdz_end*(0.0 - p2.z()) ;
      mytrk._dfdz   = dfdz_end;

      if (_debug > 0) {
	printf("[HelixFitHack:calculateTrackParameters]              POINTS USED      \n");
	printf("[HelixFitHack:calculateTrackParameters]    %5.3f   %5.3f  %5.3f\n", p1.x(), p1.y(), p1.z());
	printf("[HelixFitHack:calculateTrackParameters]    %5.3f   %5.3f  %5.3f\n", p2.x(), p2.y(), p2.z());
	printf("[HelixFitHack:calculateTrackParameters]    %5.3f   %5.3f  %5.3f\n", p3.x(), p3.y(), p3.z());
	//     printf("[HelixFitHack:calculateTrackParameters]    %5.3f   %5.3f  %5.3f\n", x_m,y_m,0.0);
	//     printf("[HelixFitHack:calculateTrackParameters]    %5.3f   %5.3f  %5.3f\n",x_n,y_n,0.0 );
	printf("[HelixFitHack:calculateTrackParameters] x0 = %5.3f y0 = %5.3f radius = %5.3e phi_0 = %5.3f dfdz = %5.5f\n",
	       x0_end, y0_end, radius_end, mytrk._fz0, dfdz_end);
      }
    }
    
    //    if(countGoodPoints <2) chi2 = 1e10;

    if(cleanPattern){
      p0 = Hep3Vector(x0_end, y0_end, 0.0);
      for (int i=0; i<np; i++) {
	if (xyzp[i].isOutlier()) goto NEXT;

	deltaZ = xyzp[i]._pos.z() - p2.z();
	phi = phi0_end + (deltaZ)*dfdz_end;//tanLambda/radius;

	dx = p0.x() + radius_end*std::cos(phi);
	dx -=  xyzp[i]._pos.x();
	dx2 = dx*dx;
      
	dy = p0.y() + radius_end*std::sin(phi);
	dy -=  xyzp[i]._pos.y();
	dy2 = dy*dy;
	dist = dx2 + dy2;
	if (_debug > 0) {
	  printf("[HelixFitHack:doCleanUp] index = %i dist = %5.3e\n",i,dist);
	}
	_hDist->Fill(dist);
	if( (dist > distGoodPoint ) &&
	    (i != seedIndex) &&
	    (i != hack->fData[3]) ){
	  ++k;
	  xyzp[i].setOutlier();
	}
      NEXT:;
      }//end loop on strawhits for doing cleanup
    }

    
  END:;
  }
  

  void
  HelixFitHack::calculateTrackParameters(Hep3Vector& p0, double&radius,
					 double& phi0, double& tanLambda,
					 Hep3Vector p1, Hep3Vector p2,
					 Hep3Vector p3,
					 HelixFitHackResult& mytrk,
					 bool cleanPattern){
    p0.setZ(p2.z());
    
    double x_m,y_m, x_n, y_n;
    //coordianates of the mean point between p1 and p3
    x_m = (p3.x() + p1.x())/2.;
    y_m = (p3.y() + p1.y())/2.;
    //------------------------------------------------------------//

    //coordianates of the mean point between p2 and p3    
    x_n = (p3.x() + p2.x())/2.;
    y_n = (p3.y() + p2.y())/2.;
    //------------------------------------------------------------//

    //calculate now the term of the line ortoghonal to the mid point of
    //the cord which links p1 and p3
    double m = -1.*(p3.x() - p1.x())/(p3.y() - p1.y());
    double c = y_m - x_m*m;
    //the eq. is: y = x*m + c

    //calculate now the term of the line ortoghonal to the mid point of
    //the cord which links p2 and p3
    double k = -1.*(p3.x() - p2.x())/(p3.y() - p2.y());
    double t = y_n - x_n*k;
    //the eq. is: y = x*k + t
    
    //now we can calculate the x0 and y0 of p0
    double x0 = (t - c)/(m - k);//(c - t) * (k*m)/(m-k);
    p0.setX(x0);
    double y0 = m*x0 + c;//(c - t) * m / (m - k) + t;
    p0.setY(y0);
    
    //now calculate the radius,phi0, tanLambda assuming that the helix
    //crosses the point (0,0,z), which is the center of the capture
    //targets
    double deltaX = p3.x() - x0;
    double deltaY = p3.y() - y0;
    double deltaZ = p3.z() - p2.z();
    //    radius    = std::sqrt(x0*x0 + y0*y0);
    radius    = std::sqrt((x0 - p3.x())*(x0 - p3.x()) + (y0 - p3.y())*(y0 - p3.y()));

    double delta2Y = (p2.y() - y0);
    double delta2X = (p2.x() - x0);
    phi0      = std::acos( delta2X/radius);
    //now compute phi0 using y coordinates so to estabilish which sign phi0 must have
    double phi0_fromY = std::asin( delta2Y/radius);
    //    phi0     = phi0*phi0_fromY/std::fabs(phi0_fromY);
    double phi0_fromXY = std::atan2(  delta2Y, delta2X);//std::atan2(  (p2.y() - y0), (p2.x() - x0));
    //2014 - 02 - 01 gianipez resolved a bug?
    //if( phi0*phi0_fromY < 0.)
    //  phi0     = phi0*phi0_fromY/std::fabs(phi0_fromY);
    if (_debug > 0) {
      printf("[HelixFitHack:calculateTrackParameters] phi0 from X = %5.3f phi0 from Y = %5.3f phi0 from tan = %5.3f\n",
	     phi0, phi0_fromY, phi0_fromXY);
    }
    phi0 = phi0_fromXY;

    double deltaPhi_0 = std::atan2(deltaY,deltaX) - phi0;
    double deltaPhi_1 = std::atan(deltaY/deltaX) - phi0;
    if (_debug > 0) {
      printf("[HelixFitHack:calculateTrackParameters] deltaPhi_0 = %5.5f deltaPhi_1 = %5.5f\n",
	     deltaPhi_0,
	     deltaPhi_1);
    }
    
    if(deltaPhi_0 < 0.) 
      deltaPhi_0 += 2.*M_PI;

    tanLambda = (radius/deltaZ)*deltaPhi_0;//(radius/deltaZ)*std::fabs(deltaPhi_0);

    if (_debug > 0) {
      printf("[HelixFitHack:calculateTrackParameters] dfdz = %5.8f \n", 
	     tanLambda/radius);
    }
    
  }
  
  void 
  HelixFitHack::calculateDfDz(double &phi0, double &phi1, 
			      double &z0,   double &z1,
			      double &dfdz){
    double deltaPhi = phi1 - phi0;
    if(deltaPhi < 0.0) deltaPhi += 2.*M_PI;
    dfdz = deltaPhi/(z1 - z0);
    if(dfdz>0.0) 
      if (_debug > 0) {
	printf("[HeliFitHack::calculateDfDZ] phi_0 = %5.3f phi_1 = %5.3f z0 = %5.3f z1 = %5.3f dfdz = %5.5f\n", 
	       phi0, phi1, z0, z1, dfdz);
      }
  }
  
}



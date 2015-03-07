///////////////////////////////////////////////////////////////////////////////
// helix fit to straw hits
//
// $Id: HelixFitHack.cc,v 1.13 2014/06/06 21:35:08 murat Exp $
// $Author: murat $ 
// $Date: 2014/06/06 21:35:08 $
//
//  use of HackData:
//  ----------------
//  [00:13] : Giani
//  [14:15] : Pasha - parameters of fit with non-equal weights
//
// the following has to come before other BaBar includes
///////////////////////////////////////////////////////////////////////////////
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
#include "TVector2.h"
// C++
#include <vector>
#include <string>
#include <algorithm>

#include "CalPatRec/inc/THackData.hh"

namespace mu2e 
{
  // statics
  double       XYZPHack::_efac(1.0);
  StrawHitFlag XYZPHack::_useflag;

  // comparison functor for ordering points
  struct radcomp : public std::binary_function<VALERR, VALERR, bool> {
    bool operator()(VALERR const& r1, VALERR const& r2) { return r1._val < r2._val; }
  };

  // comparison functor for sorting by z
  struct zcomp : public std::binary_function<XYZPHack,XYZPHack,bool> {
    bool operator()(XYZPHack const& p1, XYZPHack const& p2) { return p1._pos.z() < p2._pos.z(); }
  };

  XYZPHack::XYZPHack(size_t index, 
		     StrawHit const& sh, 
		     StrawHitPosition const& shp, 
		     Straw const& straw,
		     StrawHitFlag const& flag) :
    _ind (index), 
    _pos (shp.pos()), 
    _phi (shp.pos().phi()), 
    _flag(flag),
    _used(0),
    _wdir(straw.getDirection()),
    _straw(&straw),
    _strawhit(&sh),
    _perr(_efac*shp.posRes(StrawHitPosition::phi)),
    _rerr(_efac*shp.posRes(StrawHitPosition::rho))
  {
    static const CLHEP::Hep3Vector _zdir(0.0,0.0,1.0);
    _sdir = _zdir.cross(_wdir);
  }

  XYZPHack::XYZPHack(size_t ind, 
		     CLHEP::Hep3Vector const& pos, 
		     CLHEP::Hep3Vector const& wdir, 
		     double werr, 
		     double serr) :
    _ind(ind),
    _pos(pos),
    _phi(_pos.phi()),
    _flag(),
    _used(0),
    _wdir(wdir),
    _sdir(wdir.y(),-wdir.x(),0.0),
    _straw(0),
    _strawhit(0),
    _perr(_efac*werr),
    _rerr(_efac*serr)
  {
  }
  
  void XYZPHack::rinfo(CLHEP::Hep3Vector const& center, VALERR& rad) const {
    //    static const double onethird(1.0/3.0);
    //    static const double invsqrt12(1./sqrt(12.0));
    // average the 1-sigma radii to account for non-linear errors
    double rvec = CLHEP::Hep3Vector(_pos - center).perp();
    //    rad._val = onethird*(rvec+rvec1+rvec2);
    rad._val = rvec;
    rad._err = _rerr;
    
  }

  void XYZPHack::finfo(CLHEP::Hep3Vector const& center,VALERR& phi) const {
    //    static const double onethird(1.0/3.0);
    //    static const double invsqrt12(1./sqrt(12.0));
    // average the 1-sigma radii to account for non-linear errors
    double phi0 = CLHEP::Hep3Vector(_pos - center).phi();
    //    rad._val = onethird*(rvec+rvec1+rvec2);
    phi._val = phi0;
    phi._err = _perr; 
  }

  bool XYZPHack::use() const { 
    return !_flag.hasAnyProperty(_useflag);
  }

  bool XYZPHack::stereo() const {
    static StrawHitFlag stereo(StrawHitFlag::stereo);
    return _flag.hasAllProperties(stereo);
  }

  void XYZPHack::setUse(bool use) {
    static StrawHitFlag other(StrawHitFlag::other);
    if(!use)
      _flag.merge(other);
    else
      _flag.clear(other);
  }

  void XYZPHack::setOutlier(){
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    _flag.merge(outlier);
  }

  bool XYZPHack::isOutlier() const {
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    return _flag.hasAllProperties(outlier);
  }

  bool XYZPHack::isCalosel() const {
    static StrawHitFlag calosel(StrawHitFlag::calosel);
    return _flag.hasAllProperties(calosel);
  }

  //-----------------------------------------------------------------------------
  // HelixFitHackResult
  //-----------------------------------------------------------------------------
  HelixFitHackResult& 
  HelixFitHackResult::operator =(HelixFitHackResult const& other) {
    if(this != &other){
      _hdef   = other._hdef;
      _fit    = other._fit;
      _center = other._center;
      _radius = other._radius;
      _dfdz   = other._dfdz;
      _fz0    = other._fz0;
    }
    return *this;
  }
  
//-----------------------------------------------------------------------------
  void HelixFitHack::helixParams(HelixFitHackResult const& helix,
				 CLHEP::HepVector&         pvec,
				 CLHEP::HepVector&         perr) const 
  {
    HelixDefHack const& mytrk = helix._hdef;
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

  //-------------------------------------------------------------------------//
  //  2014-12-26 Gianipez added the following method for asking t the helix
  // finder if the hit "index" has already been marked as used by a previous 
  // by a previous search
  int   HelixFitHack::isHitUsed(int index)
  {
    if ( (_goodPointsTrkCandidate < _minPointsTrkCandidate) ||
	 //(_chi2TrkCandidate > _maxChi2TrkCandidate) ||
	 _markCandidateHits == 0) return 0;
    
    if (index >=400) {
      printf("[HelixFitHack::isHitUsed] requested index = %i range exceeded the range allowed\n", index);
      return 1;
    }

    if (_indicesTrkCandidate[index] <= 0) {
      return 0;
    }else {
      return 1;
    }
  }

//-----------------------------------------------------------------------------
  HelixFitHack::HelixFitHack(fhicl::ParameterSet const& pset) :
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
    _mpDfDz(pset.get<double>("mostProbableDfDz")),
    _maxDfDz(pset.get<double>("maxDfDz",0.01)),
    _minDfDz(pset.get<double>("minDfDz",5e-04)),
    _distPatRec(pset.get<double>("distPatRec")),
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
    _bz(0.0),
    _x0(-9999.),
    _y0(-9999.), 
    _phi0(-9999.), 
    _radius(-9999.), 
    _dfdz(-9999.),
    _goodPointsTrkCandidate(-9999),
    _minPointsTrkCandidate(pset.get<int>("minPointsTrkCandidate")),
    _chi2TrkCandidate(1e10),
    _maxChi2TrkCandidate(pset.get<double>("maxChi2TrkCandidate")),
    _markCandidateHits(pset.get<int>("markCandidateHits")),
    _chi2xyMax(pset.get<double>("chi2xyMax")),
    _chi2zphiMax(pset.get<double>("chi2zphiMax")),
    _dfdzErr(pset.get<double>("dfdzErr")){
    XYZPHack::_efac = _efac;
    std::vector<std::string> bitnames;
    bitnames.push_back("Outlier");
    bitnames.push_back("OtherBackground");
    XYZPHack::_useflag = StrawHitFlag(bitnames);

    for (int i=0; i<400; ++i){
      _indicesTrkCandidate[i] = -9999;
      _distTrkCandidate[i]    = -9999;
      _dzTrkCandidate[i]      = -9999;
    }

    _hDist         = NULL;
    _chi2nFindZ    = 0.0,
    _eventToLook   = -1;
    
    _hDfDzRes      = NULL;//new TH1F("hDfDzRes","dfdz residuals",
    //  20, _minDfDz, _maxDfDz);
  }


//-----------------------------------------------------------------------------
  HelixFitHack::~HelixFitHack() {
    // 2014-01-05 gianipez added a diag level using a TH1F
    //     _fOut->cd();
    if(_hDist){
      delete _hDist;
      _hDist = 0;
    }
    if (_hDfDzRes){
      delete _hDfDzRes;
      _hDfDzRes = 0;
    }
  }

//-----------------------------------------------------------------------------
  double HelixFitHack::bz() const {
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

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  bool HelixFitHack::findHelix(HelixFitHackResult& myhel, const CalTimePeak* TimePeak) {

    fTimePeak = TimePeak;

    HelixDefHack const& mytrk = myhel._hdef;

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
    fillXYZP(mytrk);
				        // 2013-10-18 P.Murat: print 
    int n = _xyzp.size();

    XYZPHack* pt;
    if (_debug > 0) {
      printf("[HelixFitHack::findHelix] x0 = %12.5f y0 = %12.5f r = %12.5f chi2 = %12.5g\n",
	     0.,0.,0.,-1.);
    }

    for(int i=0; i<n; ++i) {
      pt = &_xyzp[i];

      if (_debug > 0) {
	printf("[HelixFitHack::findHelix] %08x %2i %12.5f %12.5f %12.5f \n",
	       *((int*) &pt->_flag), pt->use(), pt->_pos.x(), pt->_pos.y(), pt->_pos.z()
	       );
      }
    }
//-----------------------------------------------------------------------------
// call down
//-----------------------------------------------------------------------------
    bool retval = findHelix(_xyzp,myhel);

    if ((retval || _plotall) && _diag>1) {
					// fill graphs for display if requested
      plotXY(mytrk,_xyzp,myhel);
      plotZ (mytrk,_xyzp,myhel);
    }
    return retval;
  }
//---------------------------------------------------------------------------
// reset track paramters
//---------------------------------------------------------------------------
  void HelixFitHack::resetTrackParamters(){
//indices on the xyzp vector of: the straw hit seeding the search, the second strawhit used for recalculating the dfdz value
// fLastIndex which is an helpfool paramter for searching goodpoints during the pattern recognition
    fSeedIndex      = -9999;   
    fCandidateIndex = -9999;
    fLastIndex      = -9999;   
    fUseDefaultDfDz = 0;
//follow helix paramters
    _x0     = -9999.;
    _y0     = -9999.;
    _phi0   = -9999.;
    _radius = -9999.;
    _dfdz   = -9999.;

//quality paramters used for doing comparison between several track candidates
    _goodPointsTrkCandidate = -99999;
    _chi2TrkCandidate       = 1e10;

//vector which holds indices of the strawhit far from the predicted position
    for(int i=0; i<400; ++i){
      _indicesTrkCandidate[i] = -9999;
      _distTrkCandidate[i]    = -9999;
      _dzTrkCandidate[i]      = -9999;
    }

  }

//-----------------------------------------------------------------------------
// called internally
//-----------------------------------------------------------------------------
  bool HelixFitHack::findHelix(XYZPHackVector& xyzp, HelixFitHackResult& myhel) {
    bool retval(false);

//2014-11-09 gianipez added the following function to reset the track candidate paramters
// a new time peak is used! so the previous candidate should not be compared to the new 
// one at this level
    resetTrackParamters();

    if (_filter) filterDist(xyzp);

    doPatternRecognition(xyzp,myhel);
//---------------------------------------------------------------------------
// 2014-11-11 gianipez changed the following if() statement to test the
// possibility of spead up the pattern recognition in presence of background
//---------------------------------------------------------------------------
    
    if (_debug !=0){
      printf("[HelixFitHack::findHelix] myhel._sxy.qn() = %5.0f goodPointsTrkCandidate = %i\n", 
	     myhel._sxy.qn(), _goodPointsTrkCandidate);
    }
    
    if (myhel._sxy.qn() < _minnhit) {
    //    if (_goodPointsTrkCandidate< _minnhit) {
      myhel._fit = TrkErrCode(TrkErrCode::fail,1); // small number of hits
    }
    else if ((myhel._radius < _rmin) || (myhel._radius > _rmax)) {
      myhel._fit = TrkErrCode(TrkErrCode::fail,2); // small number of hits
    }
    else {
				// success
      retval = true;
    }
    return retval;
  }
 
//-----------------------------------------------------------------------------
  bool HelixFitHack::findXY(XYZPHackVector& xyzp, HelixFitHackResult& myhel) {
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
  bool  HelixFitHack::findXY_new(XYZPHackVector& xyzp, HelixFitHackResult& Hel,
				 int seedIndex, int *indexVec) {
    double chi2_min, chi2, x, y;
    int    iworst, np, pointsRemoved(0);
    bool   success(0);

    chi2_min = Hel._sxy.chi2DofCircle();
    if (_debug > 5) {
      printf("[HelixFitHack::findXY_new] x0 = %12.5f y0 = %12.5f r = %12.5f chi2 = %12.5g\n",
	     Hel._sxy.x0(),Hel._sxy.y0(),Hel._sxy.radius(),Hel._sxy.chi2DofCircle());
    }
    if (chi2_min < _chi2xyMax) {

      if (_debug > 5) {

	printf("[HelixFitHack::findXY_new] don't need to restart, fit is already good!\n");
	printf("[HelixFitHack::findXY_new] x0 = %12.5f y0 = %12.5f r = %12.5f chi2 = %12.5g\n",
	       Hel._sxy.x0(),Hel._sxy.y0(),Hel._sxy.radius(),Hel._sxy.chi2DofCircle());

	np = xyzp.size();
	for (int i=0; i<np; i++) {
	  printf("[HelixFitHack::findXY_new] %08x %2i %12.5f %12.5f %12.5f \n",
		 *((int*) &xyzp[i]._flag), indexVec[i], 
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
    {   ::LsqSums4 sxy;
    np       = xyzp.size();


    double weight(0.0);
  NEXT_ITERATION:;
    chi2_min = 1e06;
    iworst   = -1;
    //    for (int i=0; i<np; i++) {
    for (int i=seedIndex+1; i<np; i++) {
      if (xyzp[i].isOutlier()) goto NEXT_POINT;

      //2015-01-25 G.Pezzu avoid the use of hit rejected by the helix search
      if (indexVec[i] < 1)     goto NEXT_POINT;

      sxy.init(Hel._sxy);

      x = xyzp[i]._pos.x();
      y = xyzp[i]._pos.y();
      weight = 1.0;

      sxy.removePoint(x,y, weight);

      chi2 = sxy.chi2DofCircle();

      if (chi2 < chi2_min) {
	iworst   = i;
	chi2_min = chi2;
      }
    NEXT_POINT:;
    }

    if(iworst>=0){
      x = xyzp[iworst]._pos.x();
      y = xyzp[iworst]._pos.y();
      weight = 1.;
      Hel._sxy.removePoint(x,y, weight);
      indexVec[iworst] = 0;

      ++pointsRemoved;
    }
  
    if (chi2_min >= _chi2xyMax) {
      //-----------------------------------------------------------------------------
      // still bad chi2, repeat the cleanup cycle
      //-----------------------------------------------------------------------------
      if (Hel._sxy.qn() > 10.) {
	goto NEXT_ITERATION;
      }
    }
    else {
      success = true;


      if (_debug > 5) {
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
    }
    Hel._center.set(Hel._sxy.x0(),Hel._sxy.y0(),0.);
    Hel._radius = Hel._sxy.radius();
    if (_debug > 5) {
      printf("[HelixFitHack::findXY_new] retval = %d points removed = %i\n",success ? 1:0, pointsRemoved);
    }
    return success;
  }


//-----------------------------------------------------------------------------
// plot functions
//-----------------------------------------------------------------------------
  void HelixFitHack::plotXY(HelixDefHack const&       mytrk, 
			    XYZPHackVector const&         xyzp,
			    HelixFitHackResult const& myhel) const 
  {
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

//-----------------------------------------------------------------------------
  bool HelixFitHack::findCenterAGE(XYZPHackVector const& xyzp,
				   Hep3Vector&       center, 
				   double&           rmed, 
				   double&           age,
				   bool              useweights) 
  {
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
	if(_debug > 5 && agenew-age>minage)
	  std::cout << "iteration did not improve AGE!!! lambda = " 
		    << lambda  << " age = " << age << " agenew = " << agenew << std::endl;
	break;
      }
      ++niter;
    }
    // check for convergence
    if(_debug > 5 && niter > _maxniter ){
      std::cout << "AGE didn't converge!!! " << std::endl;
      //      return false;
    }
    return true;
  }
//----------------------------------------------------------------------------------------
// 2015-01-13  calculate the dfdz of the track using distribution of the dfdz residuals
//----------------------------------------------------------------------------------------
  void HelixFitHack::findDfDz(XYZPHackVector& xyzp, HelixFitHackResult& myhel, 
			      int seedIndex, int *indexVec) {
    
    double phi, phi_ref(-1e10), z, z_ref, dphi, dz;
    
    _hDfDzRes = new TH1F("hDfDzRes","dfdz residuals",
		      20, _minDfDz, _maxDfDz);
    
    mu2e::GeomHandle<mu2e::TTracker> ttHandle;
    const mu2e::TTracker* tracker = ttHandle.get();

    CLHEP::Hep3Vector center = myhel._center;
    CLHEP::Hep3Vector pos_ref, pos;

    if (_debug >5){
      printf("[HelixFitHack::findDfDz] x0 = %9.3f y0 = %9.3f radius = %9.3f dfdz = %9.6f straw-hits = %9.5f\n",
	     center.x(), center.y(), myhel._radius, myhel._dfdz, (myhel._sxy.qn() - 1) );// -1 to remove the EMC cluster contribute
      printInfo(myhel);
    }

    int  np, ist, nstations;

    np        = xyzp.size();
    nstations = tracker->nStations();
    
    double phiVec[30], zVec[30];
    int    nhits [30];

    for (int i=0; i<nstations; i++) {
      phiVec[i] = 0;
      zVec  [i] = 0;
      nhits [i] = 0;
    }
//-----------------------------------------------------------------------------    
// Part 1: use only contiguous parts of the trajectory
//-----------------------------------------------------------------------------
    for (int i=seedIndex; i<np; i++) {
      if ( (! xyzp[i].isOutlier()) && (indexVec[i] >0 )) {
	// didn't find an accessor returning the station number, hack
	ist = xyzp[i]._straw->id().getDevice()/2;
	pos = xyzp[i]._pos;
	phi = CLHEP::Hep3Vector(pos - center).phi();
	zVec  [ist] += pos.z();
	if (nhits == 0) phiVec[ist] = phi;
	else {
	  while (phi-phiVec[ist] >  M_PI) phi -= 2*M_PI;
	  while (phi-phiVec[ist] < -M_PI) phi += 2*M_PI;

	  phiVec[ist] = (phiVec[ist]*nhits[ist]+phi)/(nhits[ist]+1);
	}
	nhits [ist] += 1;
      }
    }

    for (int i=0; i<nstations; i++) {
      if (nhits[i] > 0) {
	zVec  [i] = zVec  [i]/nhits[i];
      }
    }

    double tollMin(100.), tollMax(800.);

    int i0(-1), first_point(1);

    for (int i=0; i<nstations; i++) {
      if (nhits[i] == 0) goto NEXT_POINT;
				        // find station corresponding to the first point
      if (first_point) {
	i0 = i;
	first_point = 0;
      }

      phi_ref = phiVec[i];
      z_ref   = zVec  [i];

      for(int j=i+1; j<nstations; ++j){
	if (nhits[j] == 0) continue;
	phi = phiVec[j];
	z   = zVec  [j];
	dz  = z - z_ref;

	if ( (phi_ref > -9999 ) && 
	     (  dz < tollMax  ) &&
	     (  dz > tollMin  )){
	  dphi    = phi - phi_ref;
	  while (dphi >  M_PI) dphi -= 2*M_PI;
	  while (dphi < -M_PI) dphi += 2*M_PI;
	  
	  _hDfDzRes->Fill(dphi/dz);
	    
	  if (_debug >5){
	    printf("[HelixFitHack::findDfDz] z_ref = %9.3f z = %9.3f phi_ref = %9.5f phi = %9.5f dz = %9.5f df/dz = %9.5f\n",
		   z_ref, z, phi_ref, phi, dz, dphi/dz);
	  }
	}
      }
	  
    NEXT_POINT:;
    }
    //    int maxBin  = _hDfDzRes->GetMaximumBin();
    _hdfdz      = _hDfDzRes->GetMean();

    int nentries = _hDfDzRes->GetEntries();

    if (_debug >5){
      
      printf("[HelixFitHack::findDfDz] nentries = %i mpvDfDz = %5.6f under: %5.0f over: %5.0f \n",
	     nentries, _hdfdz,
	     _hDfDzRes->GetBinContent(0),_hDfDzRes->GetBinContent(_hDfDzRes->GetNbinsX()+1)
	     );
      for (int i=0; i<_hDfDzRes->GetNbinsX(); i++) {
	printf(" %5.0f",_hDfDzRes->GetBinContent(i+1));
      }
      printf("\n");
    }
    delete     _hDfDzRes;
    _hDfDzRes = 0;
//-----------------------------------------------------------------------------
// Part 2: try to perform a more accurate estimate - straight line fit
//-----------------------------------------------------------------------------
    double z0, phi0, dphidz, pred;

    z0   = zVec  [i0];
    phi0 = phiVec[i0];

    dphidz = _hdfdz;

    double sx(0), sy(0), sx2(0), sxy(0), sy2(0), sn(0);
    int    iz;

    if (_debug >5) {
      printf("[HelixFitHack::findDfDz] ------------ Part 2:\n");
    }
    for (int i=i0; i<nstations; i++) {
      if (nhits[i] > 0) {
	z    = zVec[i];
	pred = phi0+(z-z0)*dphidz;
	iz   = (pred-phiVec[i])/(2*M_PI)+0.5;
	phi  = phiVec[i]+iz*(2*M_PI);
	//	sxy.AddPoint(z,phi);

	if (_debug >5) {
	  printf("[HelixFitHack::findDfDz] i=%3i z=%9.3f phiVec[i]=%9.5f iz=%2i phi=%9.5f\n",
		 i,z,phiVec[i],iz,phi);
	}
	
	sn  += 1;
	sx  += z;
	sy  += phi;
	sx2 += z*z;
	sxy += z*phi;
	sy2 += phi*phi;
      }
    }

    double xmean, ymean, x2mean, xymean, y2mean, sigxx, sigxy, sigyy;

    xmean  = sx/sn;
    ymean  = sy/sn;
    x2mean = sx2/sn;
    xymean = sxy/sn;
    y2mean = sy2/sn;

    sigxx  = x2mean-xmean*xmean;
    sigxy  = xymean-xmean*ymean;
    sigyy  = y2mean-ymean*ymean;

    _hdfdz = sigxy/sigxx;
    _hphi0 = ymean - xmean*sigxy/sigxx;
    _sdfdz = sigyy-_hdfdz*sigxy;

    if (_debug >5) {
      printf("[HelixFitHack::findDfDz] END: _hdfdz = %9.5f chi2 = %9.3f\n",_hdfdz,_sdfdz);
    }
    //    return (_hDfDzRes->GetEntries() > 0 ) ? true : false;
    if (nentries == 0) _hdfdz = _mpDfDz;
  }
    



//-----------------------------------------------------------------------------
  bool HelixFitHack::findZ(XYZPHackVector& xyzp, HelixFitHackResult& myhel, 
			   int seedIndex       , int *indexVec) {
    // sort points by z
    std::sort(xyzp.begin(),xyzp.end(),zcomp());

    bool success(false);

    int offset = 1;
    if (_debug > 5) printf("HelixFitHack::FindZ offset = %d\n",offset);

    int N = xyzp.size();

    double phi_corrected[N];
    double phi(0.0), z(0.0), weight(0.0);
    
    CLHEP::Hep3Vector center = myhel._center;
    CLHEP::Hep3Vector pos;

    // calculate the mean pitch assuming conversion electron 
    // hypothesis (so assuming a pitch angle of ~ 0.65)

    double pitch = twopi*std::pow(myhel._radius,2.)*myhel._dfdz; //*tan(0.65);

    if (_debug > 5) printf("[HelixFitHack::findZ] pitch = %5.3f\n", pitch);
 
    //    printf(">>HelixFit.cc:\n");
    //     printf("z-pitch = %5.3f\n",pitch);

    double meanDphi(0.0);
    //    int    realSize(0);


//-----------------------------------------------------------------------------
// gianipez: procedure for alligning the phi vector
//-----------------------------------------------------------------------------
    ::LsqSums4 srphi;
    int        iworst, count(0), indexWorst;
    double     chi2,chi2min, dfdz, deltaPhi;


    if (_debug > 5) printf("[HelixFitHack::findZ] xyzp.size() = %u\n", N);



//--------------------------------------------------------------------------------
// set EMC cluster info
    double zCl   = fTimePeak->ClusterZ();
    pos          = Hep3Vector(fTimePeak->ClusterX(), 
			      fTimePeak->ClusterY(), 
			      fTimePeak->ClusterZ());
    double phiCl = CLHEP::Hep3Vector(pos - center).phi();
    if (phiCl<0.0) phiCl +=twopi;
    // add the cluster phi to the LSq sum
    weight       = 10.0;
    myhel._srphi.clear();
    myhel._srphi.addPoint(zCl,phiCl,weight);

    //initilize the dfdz for the search
    dfdz = myhel._dfdz;

    count = 0;

    if (_debug > 5) printf("[HelixFitHack::findZ] %5.3f     %5.3f   \n", zCl, phiCl);

    for(int i=N-1; i>=seedIndex; --i){
      
      
      weight               = 1.;

      pos                  = xyzp[i]._pos;
      z                    = pos.z();
      
      phi = CLHEP::Hep3Vector(pos - center).phi();
      phi = TVector2::Phi_0_2pi(phi);

      //      if (phi<0.0) phi +=twopi;

      if (count>=2) dfdz = myhel._srphi.dfdz();

      deltaPhi = phiCl - phi;
      meanDphi = (zCl - z)*dfdz;

      while (deltaPhi - meanDphi >  M_PI) {
	phi += 2*M_PI;
	deltaPhi = phiCl - phi;
      }
      while (deltaPhi - meanDphi < -M_PI) {
	phi -= 2*M_PI;
	deltaPhi = phiCl - phi;
      }

      //store the corrected value of phi
      phi_corrected[i] = phi;

      if (_debug > 5) printf("[HelixFitHack::findZ] %08x %2i %6i %12.5f %12.5f \n",
			     *((int*) &xyzp[i]._flag), indexVec[i] < 0 ? 0 : indexVec[i], int(xyzp[i]._ind),  z, phi_corrected[i]);
      
      if (xyzp[i].isOutlier())                        goto NEXT_POINT;
      if (  indexVec[i] < 1  )                        goto NEXT_POINT;

      myhel._srphi.addPoint(z, phi_corrected[i], weight);

      ++count;

    NEXT_POINT:;
    }

    if (_debug > 5) {
      printf("[HelixFitHack::findZ] myhel: phi_0 = %5.3f dfdz = %5.5f chi2N = %5.3f\n", 
	     myhel._srphi.phi0(),myhel._srphi.dfdz(), myhel._srphi.chi2rphiDofCircle() );
    }
    //-----------------------------------------------------------------------------
    // perform a cleanup
    //-----------------------------------------------------------------------------
    if ( myhel._srphi.chi2rphiDofCircle() > _chi2zphiMax){
    NEXT_ITERATION:;
      iworst     = -1;
      indexWorst = -1;
      chi2min    = 1e10;
      //    chi2min = myhel._srphi.chi2rphiDofCircle();
      for(int ixyzp=seedIndex+1; ixyzp < N; ++ixyzp){
	if (xyzp[ixyzp].isOutlier())     goto NEXT_P;
	if (  indexVec[ixyzp] < 1  )     goto NEXT_P;
	srphi.init(myhel._srphi);
	pos = xyzp[ixyzp]._pos;
	z   = pos.z();//z_vec[i];
	phi = phi_corrected[ixyzp];//CLHEP::Hep3Vector(pos - center).phi();//phi_vec[i];
	weight = 1.;
	srphi.removePoint(z, phi, weight);
	chi2 = srphi.chi2rphiDofCircle();
	//printf("[HelixFitHack::findZ] chi2 = %5.3e chi2min = %5.3e\n", chi2, chi2min);
	if (chi2 < chi2min) {
	  iworst     = ixyzp;
	  indexWorst = xyzp[ixyzp]._ind;
	  chi2min    = chi2;
	}
      NEXT_P:;
      }
  
      if(iworst>=0 && myhel._srphi.qn() > 3.){
	indexVec[iworst] = 0;//.setOutlier();
	pos = xyzp[iworst]._pos;
	z   = pos.z();
	phi = phi_corrected[iworst];//CLHEP::Hep3Vector(pos - center).phi();
    
	weight = 1.;
	myhel._srphi.removePoint(z, phi, weight);
	chi2min = myhel._srphi.chi2rphiDofCircle();
	if (_debug > 5) {
	  printf("[HelixFitHack::findZ_removed] %6i %5.3f     %5.3f chi2 = %5.3f  \n", indexWorst, z, phi, chi2min);
	}
      }
      if(myhel._srphi.qn()<=3) chi2min = myhel._srphi.chi2rphiDofCircle();
  
      if (chi2min >= 1.) {
      //      if (chi2min > 0.1) {
    
	if (myhel._srphi.qn() > 10.) {
	  goto NEXT_ITERATION;
	}
      }
    }

    if ( myhel._srphi.chi2rphiDofCircle() < _chi2zphiMax){
     success = true;
   }
    //----------------------------------------------------------------------//
    
    myhel._fz0  = myhel._srphi.phi0();
    myhel._dfdz = myhel._srphi.dfdz();

    if ( myhel._dfdz < 0.) 
      success = false;
    
    if (_debug > 5) {
      printf("[HelixFitHack::findZ] myhel: phi_0 = %5.3f dfdz = %5.5f chi2N = %5.3f\n", 
	     myhel._srphi.phi0(),myhel._srphi.dfdz(), myhel._srphi.chi2rphiDofCircle() );
      printf("[HelixFitHack::findZ] srphi: phi_0 = %5.3f dfdz = %5.5f chi2N = %5.3f\n", 
	     srphi.phi0(), srphi.dfdz(), srphi.chi2rphiDofCircle() );
    }
    
    _chi2nFindZ = myhel._srphi.chi2rphiDofCircle();
    if(_chi2nFindZ < 0.0)  _eventToLook = myhel._hdef.eventId();

    if (_debug > 5) printf("[HelixFitHack::findZ] retval = %d\n",success ? 1:0);
   
    return success;
  }

//-----------------------------------------------------------------------------
  void  HelixFitHack::plotZ(HelixDefHack       const&  mytrk, 
			    XYZPHackVector     const&  xyzp , 
			    HelixFitHackResult const&  myhel) const 
{
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

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  bool HelixFitHack::initCircle(XYZPHackVector const& xyzp, HelixFitHackResult& myhel) {
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
  bool HelixFitHack::initCircle_new(XYZPHackVector const& xyzp, HelixFitHackResult& Hel) {
    bool   retval(false);
					//    static const double mind2 = _mindist*_mindist;
    double   x, y, x_0, y_0, r;
    XYZPHack*    p;
					// first point - center, with low weight
					//    Hel._sxy.addPoint(0.,0.,0.5);
    Hel._sxy.addPoint(0.,0.,0.05);
					//second point: the EMC cluster position!

    if (fTimePeak->ClusterT0() > 0.) {
      Hel._sxy.addPoint(fTimePeak->ClusterX(), fTimePeak->ClusterY(), 10. );
      if (_debug > 5) {
	printf("[HelixFitHack::initCircle_new] x_EMC = %12.5f y_EMC = %12.5f\n",
	       fTimePeak->ClusterX(), fTimePeak->ClusterY());
      }
    }

    int   n = xyzp.size();
//-----------------------------------------------------------------------------
//2013 - 12 - 25 gianipez added a paramter for weghting the difference between 
//               stereo and straw hits
//-----------------------------------------------------------------------------
    double weight(0.); // use the weight parameter to distinguish between 
                       // the stereohits and strawhits

					// don't use the outliers
    for (int i=0; i<n; i++) {
      p = (XYZPHack*) &xyzp[i];
      if (! p->isOutlier()) {
	x      = p->_pos.x();
	y      = p->_pos.y();
	weight = 1.;
	Hel._sxy.addPoint(x,y, weight);
      }
    }

    x_0  = Hel._sxy.x0();
    y_0  = Hel._sxy.y0();
    r    = Hel._sxy.radius();

    if (_debug > 5) {
      printf("[HelixFitHack::initCircle_new] x0 = %12.5f y0 = %12.5f r = %12.5f chi2 = %12.5g\n",
	     x_0,y_0,r,Hel._sxy.chi2DofCircle());

      for (int i=0; i<n; i++) {
	p = (XYZPHack*) &xyzp[i];
	printf("[HelixFitHack::initCircle_new] %08x %2i %12.5f %12.5f %12.5f \n",
	       *((int*) &p->_flag), p->use(), p->_pos.x(), p->_pos.y(), p->_pos.z()
	       );
      }
    }

    Hel._center.set(x_0,y_0,0.0);
    Hel._radius = r;
    
    retval = (r >= _rmin) && (r <= _rmax);
    if (_debug > 5) printf("[HelixFitHack::initCircle_new] retval = %d\n",retval ? 1:0);
    return retval;
  }
//-----------------------------------------------------------------------------
// 12-09-2013 gianipez modified this procedure to avoid the doubling of the 
// same stereohitposition 
//-------------------------------------------------------------------------
  void HelixFitHack::fillXYZP(HelixDefHack const& mytrk){
      
    //clear xyzp vector
    _xyzp.clear();
    
    const Tracker& tracker = getTrackerOrThrow();
    
    const std::vector<hitIndex> shIndices = mytrk.strawHitIndices();
    
    int size = shIndices.size();

    //--------------------------------------------------------------------------------
    if (mytrk.strawHitPositionCollection() != 0) {
      int loc;
      StrawHitFlag flag;
      for (int i=0; i<size; ++i){
	loc                = shIndices[i]._index;
	flag               = mytrk.strawHitFlagCollection()->at(loc);
	StrawHit const& sh = mytrk.strawHitCollection()->at(loc);
	Straw const& straw = tracker.getStraw(sh.strawIndex());
	StrawHitPosition const& shp = mytrk.strawHitPositionCollection()->at(loc);
	
	XYZPHack pos(loc,sh,shp,straw,flag);
	if (_debug > 0) {
	  if(i == 0){
	    printf("----->HelixFitHack::fillXYXP\n");
	    printf("|   strawIndex    |   straw ID   |   stereo index   |          pos          |\n");
	  }					
	  printf("| %10.3d    | %10.3d    | %10.3d    | %5.3f, %.3f, %5.3f |\n", 
		 (int) loc, sh.strawIndex().asInt(), shp.stereoHitIndex(),
		 shp.pos().x(), shp.pos().y(), shp.pos().z());
	}
	_xyzp.push_back(pos);
      }	// end loop over the strawindexes
    }
//----------------------------------------------------------------------
// 2014-11-06 gianipez added the following line for ordering the xyzp
// strawhits along their z coordinate
//----------------------------------------------------------------------
    std::sort(_xyzp.begin(), _xyzp.end(), [ ]( const XYZPHack& lhs,
					     const XYZPHack& rhs )
	      {
		return lhs._pos.z() < rhs._pos.z();
	      } );
  }

//-----------------------------------------------------------------------------
  void HelixFitHack::findAGE(XYZPHackVector const& xyzp, 
			     Hep3Vector const& center,
			     double& rmed, 
			     double& age,
			     bool useweights) 
  {
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

//-----------------------------------------------------------------------------
  void HelixFitHack::fillSums(XYZPHackVector const& xyzp, 
			      Hep3Vector const&     center,
			      double                rmed, 
			      SUMS&                 sums,
			      bool                  useweights) 
  {
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

//----------------------------------------------------------------------------
//2015-01-17 G. Pezzullo: the following procedure looks the hit with 
// z-coordinate smaller then the seeding one and calculates distance from
// prediction in order to check if they are good or outliers
//----------------------------------------------------------------------------
  void HelixFitHack::rescueHits(XYZPHackVector& xyzp, HelixFitHackResult&  mytrk){

    double      weight, radius, phi0, dfdz, x0, y0;
    dfdz        = mytrk._dfdz;
    phi0        = mytrk._fz0 + dfdz*(xyzp[fSeedIndex]._pos.z());
    x0          = mytrk._center.x();
    y0          = mytrk._center.y();
    radius      = mytrk._radius;

    double      dx,dy,phi,dx2, dy2, max_dist;
    Hep3Vector  shPos, hePos;
    
    double deltaZ(0.); // , deltaX(0.), deltaY(0.);
    double distXY(0.0);
    double dist(0.0), dist2(0.0); //help parameter for storing strawhit position residual
    int    i_last(fSeedIndex), rescuedPoints(0);
    
    TString banner="HelixFitHack::rescueHits";

    if (_debug > 0) {
      printf("[%s] x0 = %5.3f y0 = %5.3f radius = %5.3f phi0 = %5.5f dfdz = %5.6f chi2 = %5.3f \n", banner.Data(),
	     x0, y0, radius, phi0, dfdz , mytrk._sxy.chi2DofCircle());
      printf("[%s] seedIndex = %i N-points = %5.3f\n",  banner.Data(), fSeedIndex, mytrk._sxy.qn()-1);
    }
    
    for (int i=fSeedIndex-1; i>=0; --i){
      if (xyzp[i].isOutlier())                              goto NEXT_POINT;
      weight = 1.;
      shPos  = xyzp[i]._pos;
      deltaZ = shPos.z() - xyzp[i_last]._pos.z();
      phi    = phi0 + (deltaZ)*dfdz;                     // tanLambda/radius;

      hePos  = Hep3Vector(x0 + radius*std::cos(phi),
			  y0 + radius*std::sin(phi),
			  shPos.z());

      dx  = hePos.x() - shPos.x();
      dx2 = dx*dx;
      dy  = hePos.y() - shPos.y();
      dy2 = dy*dy;
      
      dist2 = dx2 + dy2;
      dist  = std::sqrt(dist2);

      if( _debug>0){
	
	printf("[%s]   measured     %10.3f  %10.3f  %10.3f  %8i \n", banner.Data(), 
	       shPos.x(), shPos.y(), shPos.z(), i);
	printf("[%s]  predicted     %10.3f  %10.3f  %10.3f  %8i \n", banner.Data(), 
	       hePos.x(), hePos.y(), hePos.z(), i);
	printf("[%s] X0 = %5.3f Y0 = %5.3f r = %5.3f dfdz = %5.5f  dist-from-prediction = %5.3f  dist-from-seedXY = %5.3f dz-from-seed = %5.3f\n", banner.Data(), 
	       x0, y0, radius, dfdz, dist, distXY, deltaZ);
      }
 
      max_dist = _distPatRec + _dfdzErr*fabs(deltaZ);
      if( dist <= max_dist ){

	//store index of last good point
	i_last = i;
	phi0   =  phi;//CLHEP::Hep3Vector(shPos - CLHEP::Hep3Vector(x0, y0, 0.)).phi();

	//add point to the helixfithack result objet
	mytrk._sxy.addPoint(shPos.x(),shPos.y(), weight);
      
	//store the index of the good point found
	_indicesTrkCandidate[i] = 1;

	//store distance along z-axis from the last point found relying on thr helix
	_dzTrkCandidate[i]      = deltaZ;

	//store distance from predition
	_distTrkCandidate[i]    = dist;

	//update helix parameters
	x0      = mytrk._sxy.x0();
	y0      = mytrk._sxy.y0();
	radius  = mytrk._sxy.radius();

	++rescuedPoints;
	
	if( _debug>0){
	  printf("[%s] rescued %08x %2i %12.5f %12.5f %12.5f \n", banner.Data(),
		 *((int*) &_xyzp[i]._flag), _indicesTrkCandidate[i], 
		 _xyzp[i]._pos.x(), _xyzp[i]._pos.y(), _xyzp[i]._pos.z()
		 );
	}
      }

    NEXT_POINT:;
    }

    //update mytrk info
    
    mytrk._center.set(x0, y0, 0.0);
    mytrk._radius = radius;

    _goodPointsTrkCandidate = mytrk._sxy.qn() - 1;//removing the EMC cluster!

    banner += "-results";
    if (_debug > 5) {
      printf("[%s] x0 = %5.3f y0 = %5.3f radius = %5.3f dfdz = %5.6f chi2 = %5.3f \n", banner.Data(),
	     x0, y0, radius, dfdz , mytrk._sxy.chi2DofCircle());
      printf("[%s] seedIndex = %i N rescued points = %i\n",  banner.Data(), fSeedIndex, rescuedPoints);
    }
  }

//--------------------------------------------------------------------------------
  void HelixFitHack::printInfo(HelixFitHackResult& mytrk){
    TString banner="HelixFitHack::printInfo";

    if (_debug > 0) {
      printf("[%s] N - points = %3.3f\n",   banner.Data(),  mytrk._sxy.qn() - 1);
 
      printf("[%s] x0 = %12.5f y0 = %12.5f r = %12.5f chi2 = %12.5f phi0 = %5.5f dfdz = %5.6f chi2 = %5.3f\n", banner.Data(),
	     mytrk._sxy.x0(),mytrk._sxy.y0(),mytrk._sxy.radius(),mytrk._sxy.chi2DofCircle(), 
	     mytrk._fz0, mytrk._dfdz , mytrk._srphi.chi2rphiDofCircle());

      int np = _xyzp.size();
      for (int i=0; i<np; i++) {
	printf("[%s] %08x %2i %12.5f %12.5f %12.5f \n", banner.Data(),
	       *((int*) &_xyzp[i]._flag), _indicesTrkCandidate[i] < 0 ? 0 : _indicesTrkCandidate[i], 
	       _xyzp[i]._pos.x(), _xyzp[i]._pos.y(), _xyzp[i]._pos.z()
	       );
      }
    }
  }

//-----------------------------------------------------------------------------
// this routine simply checks '_indicesTrkCandidate' array and for negative 
// indices sets the 'outlier' flag to the corresponding 'xyzp'
// no actual check of residuals is performed
//-----------------------------------------------------------------------------
  void HelixFitHack::filterUsingPatternRecognition(XYZPHackVector& xyzp, HelixFitHackResult& mytrk){
    int np = xyzp.size();
    if (fSeedIndex < 0) return;
    Hep3Vector pSeed = xyzp[fSeedIndex]._pos;
    //if no second strawhit has been found, use the target center in the transverse plane
    Hep3Vector pCand(0.0, 0.0, 0.0);

    if (fCandidateIndex >=0){
      pCand = xyzp[fCandidateIndex]._pos;
    }
    Hep3Vector shPos;
    double dist, dz;
    int    nActive(0);
    //    for (int i=fSeedIndex; i<np; ++i){
    for (int i=0; i<np; ++i){
      if (_debug>0){
	//	if( i==fSeedIndex+1) {
	if( i==0) {
	  dist = 0;
	  dz   = pCand.z() - pSeed.z();
	  printf("[HelixFitHack::filterUsingPatternRecognition]  filterUsingPatternRecognition() will set asOutlier the following hits using helix parameters\n");
	  printf("[HelixFitHack::filterUsingPatternRecognition] X0 = %5.3f Y0 = %5.3f r = %5.3f chi2N = %5.5f phi0 = %5.5f dfdz = %5.5f chi2N = %5.5f straw-hits = %i\n", 
		 mytrk._sxy.x0(), mytrk._sxy.y0(), _radius, mytrk._sxy.chi2DofCircle(), mytrk._srphi.phi0(), mytrk._srphi.dfdz(), mytrk._srphi.chi2rphiDofCircle(), _goodPointsTrkCandidate );// +1 for counting also the seeding hit
	  printf("[HelixFitHack::filterUsingPatternRecognition]   point  type     X         Y        Z        xyzp-index    hit index    dist       Dz    \n");
	  printf("[HelixFitHack::filterUsingPatternRecognition] ----------------------------------------------------------\n");
	  printf("[HelixFitHack::filterUsingPatternRecognition]    seeding     %8.3f  %8.3f  %10.3f  %7i \n",
		 pSeed.x(), pSeed.y(), pSeed.z(), fSeedIndex);
	  printf("[HelixFitHack::filterUsingPatternRecognition]   candidate    %8.3f  %8.3f  %10.3f  %7i \n",
		 pCand.x(), pCand.y(), pCand.z(), fCandidateIndex);
	}
      }
      shPos = xyzp[i]._pos;
      dist  = _distTrkCandidate[i];
      dz    = _dzTrkCandidate[i];
      if ( _indicesTrkCandidate[i] <= 0 ){
	xyzp[i].setOutlier();
	
	if (_debug>0){
	  printf("[HelixFitHack::filterUsingPatternRecognition]   outlier     | %8.3f | %8.3f | %10.3f | %8i | %7i | %8.3f | %8.3f |\n", 
		 shPos.x(), shPos.y(), shPos.z(), i, xyzp[i]._strawhit->strawIndex().asInt(), dist, dz);
	}
      }else {
	++nActive;
	if (_debug>0){
	  printf("[HelixFitHack::filterUsingPatternRecognition]    active     | %8.3f | %8.3f | %10.3f | %8i | %7i | %8.3f | %8.3f |\n", 
		 shPos.x(), shPos.y(), shPos.z(), i, xyzp[i]._strawhit->strawIndex().asInt(), dist, dz);
	}
      }
    } 
    if (_debug>0){
      printf("[HelixFitHack::filterUsingPatternRecognition]    ended: N active point = %i over N-hits = %i\n", nActive, np); 
    }
    _goodPointsTrkCandidate = nActive;
  }
  
//-----------------------------------------------------------------------------
  void HelixFitHack::filterDist(XYZPHackVector& xyzp) {
    using namespace boost::accumulators;
    static const double pi(M_PI);
    static const double twopi(2*pi);

    // 2013-10-17 G.Pezzullo, P.Murat : try to change the cleanup logic
    double mphi(-9999.);
    if (fTimePeak->ClusterT0() > 0.) {
      mphi = atan2(fTimePeak->ClusterY(),fTimePeak->ClusterX());
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

//-----------------------------------------------------------------------------
  void HelixFitHack::filterXY(XYZPHackVector& xyzp, Hep3Vector const& center,double rmed,bool& changed) {
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

//-----------------------------------------------------------------------------
  double HelixFitHack::deltaPhi(double phi1, double phi2){
    static const double pi(M_PI);
    static const double twopi(2*pi);
    double dphi = fmod(phi2-phi1,twopi);
    if(dphi>pi)dphi -= twopi;
    if(dphi<-pi)dphi += twopi;
    return dphi;
  }

//-----------------------------------------------------------------------------
  void HelixFitHack::doPatternRecognition(XYZPHackVector& xyzp, HelixFitHackResult& mytrk ){
    int np = xyzp.size();
    //int seedIndex(-1);//, targetIndex(-1);
    int mode(0);                   //flag for indexing event where a second strawhit has been included in the patern reco

    double chi2;
    int countGoodPoints(0);

    if (_debug != 0) printf("[HelixFitHack::doPatternRecognition]: BEGIN\n");

    for (int i=0; i<np; i++) {
      if (xyzp[i].isOutlier()) goto NEXT_POINT;
      if (xyzp[i].isCalosel()) goto NEXT_POINT;
//----------------------------------------------------------------------
// 2014-12-26 gianipez: don't start the search from an already used hit
// used in previous search
//-----------------------------------------------------------------------------
      if ( isHitUsed(i) == 1 )  goto NEXT_POINT;

      if (_debug != 0) printf("[HelixFitHack::doPatternRecognition]: calling findTrack i=%3i\n",i);
      findTrack(xyzp, i, chi2, countGoodPoints, mytrk, mode, false); 
//------------------------------------------------------------------------------
// 2015-01-22 P.Murat: what happens when the very first candidate is good enough ?
//                     where is the comparison of the found candidate with the best previous one ?
//-----------------------------------------------------------------------------
    NEXT_POINT:;
    }
//-----------------------------------------------------------------------------
// 2014-11-09 gianipez: if no track was found requiring the recalculation of dfdz
// look for a track candidate using the default value of dfdz and the target center
//-----------------------------------------------------------------------------
    if ( fUseDefaultDfDz == 0){
      for (int i=0; i<np; i++) {
	if (xyzp[i].isOutlier()) goto NEXT_P;
	if (xyzp[i].isCalosel()) goto NEXT_P;
	if (_debug >5 ) printf("[HelixFitHack::doPatternRecognition]: fUseDefaultDfDz=0, calling findTrack i=%3i\n",i);
	if ( (np -i) > _goodPointsTrkCandidate){
	  findTrack(xyzp, i, chi2, countGoodPoints, mytrk, mode, true);
	}
      NEXT_P:;
      } 
    }

    //2015-01-14 G. Pezzullo added the findDfDz procedure
    if (_debug > 5) printf("[HelixFitHack::doPatternRecognition]: ------------ calling findDfDz\n");
    if (fSeedIndex>=0) {
      findDfDz(xyzp, mytrk, fSeedIndex, _indicesTrkCandidate);
      if (_debug >5) printf("[HelixFitHack::doPatternRecognition]: findDfDz ----> phi0 = %5.5f dfdz = %5.5f \n",
			      _hphi0, _hdfdz);
    }else {
      //if seedIndex is < 0 it means that no cancidate has been found
      // cases when it occurs are usually the one where the cluster is not in the trajectory or 
      // la very low number of hits is in the time peak
      // maybe we should set a threshold on the time peak size to avoid such?
      int np = xyzp.size();
      int vIndices[np];
      for (int i=0; i<np; ++i){
	vIndices[i] = 1;
      }
      findDfDz(xyzp, mytrk, 0, vIndices);
      if (_debug > 5) {
	printf("[HelixFitHack::doPatternRecognition]: findDfDz called using seedIndex = 0 and using all hits (expect outliers!) \n");
	printf("[HelixFitHack::doPatternRecognition]: findDfDz ----> phi0 = %5.5f dfdz = %5.5f \n",
	       _hphi0, _hdfdz);
      }
    }

    int useMPVdfdz=1;
    for (int i=0; i<np; i++) {
      if (xyzp[i].isOutlier()) goto NEXT_HIT;
      if (xyzp[i].isCalosel()) goto NEXT_HIT;
	if (_debug > 5) printf("[HelixFitHack::doPatternRecognition]: useMPVdfdz=1, calling findTrack i=%3i\n",i);
	if ( (np -i) > _goodPointsTrkCandidate){
	  findTrack(xyzp, i, chi2, countGoodPoints, mytrk, mode, false, useMPVdfdz);
	}
    NEXT_HIT:;
    } 

    // 2015-01-17 G. Pezzullo added the following procedure to rescue points with z-coordinate 
    // less than the seed hit
    if (_debug != 0) printf("[HelixFitHack::doPatternRecognition]: calling rescueHits\n");

    rescueHits(xyzp, mytrk);

//-----------------------------------------------------------------------------
// finally, assume that the found helix is already close enough and refine  
// the helix parameters accounting for different weights
//-----------------------------------------------------------------------------
    if (_debug != 0)  printInfo(mytrk);
    if (refineHelixParameters(xyzp,mytrk, 0, _indicesTrkCandidate)){
      mytrk._center.setX(mytrk._sxyw.x0());
      mytrk._center.setY(mytrk._sxyw.y0());
      mytrk._radius    = mytrk._sxyw.radius();
      mytrk._sxy.init(mytrk._sxyw);
      if (_debug != 0)  printInfo(mytrk);
    }

    findZ(xyzp, mytrk, 0, _indicesTrkCandidate);

    // 2014-11-09 gianipez changed the cleanup process. now it is faster and cleaner
    if (_debug != 0) printf("[HelixFitHack::doPatternRecognition]: calling filterUsingPatternRecognition\n");

    filterUsingPatternRecognition(xyzp, mytrk);

    if (_debug != 0) printf("[HelixFitHack::doPatternRecognition]: END\n");
  }


//-----------------------------------------------------------------------------
// 1. where is the cluster - at this point it is no longer needed
// assume the 
//-----------------------------------------------------------------------------
  int HelixFitHack::refineHelixParameters(XYZPHackVector& Xyzp, HelixFitHackResult& Trk,
					  int seedIndex, int *indexVec) {

  double    wt, x0, y0, sinth2, costh, e2, x, y, dx, dy;
  double    rs( 2.5);  // mm
  double    ew(30.0);  // mm - erro along the wire   

  int np = Xyzp.size();

  double weights[np];
  int    success = 0;
  
  x0 = Trk._sxy.x0();//_center.x();
  y0 = Trk._sxy.y0();//_center.y();

  if (_debug > 5) {
    printf("[HelixFitHack::refineHelixParameters] starts x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f \n",
	   Trk._sxy.x0(), Trk._sxy.y0(), Trk._sxy.radius(), Trk._sxy.chi2DofCircle());
  }
  if (_debug > 5) {
    printf("[HelixFitHack::refineHelixParameters] i       X        Y        dx        dy         costh        sinth2         e2 \n");
  }
//-----------------------------------------------------------------------------
// add cluster with a position error of 10 mm => wt = 1/100
//-----------------------------------------------------------------------------
  Trk._sxyw.clear();
  Trk._sxyw.addPoint(fTimePeak->ClusterX(), fTimePeak->ClusterY(), 1./100.);

  for (int i=seedIndex; i<np; i++) {
    if (Xyzp[i].isOutlier())           goto NEXT_POINT;
    if ( indexVec[i] < 1)              goto NEXT_POINT;
    x  = Xyzp[i]._pos.x();
    y  = Xyzp[i]._pos.y();
    dx = x-x0;
    dy = y-y0;

    costh  = (dx*Xyzp[i]._sdir.x()+dy*Xyzp[i]._sdir.y())/sqrt(dx*dx+dy*dy);
    sinth2 = 1-costh*costh;

    e2    = ew*ew*sinth2+rs*rs*costh*costh;
    wt    = 1./e2;
    
    weights[i] = wt;

    Trk._sxyw.addPoint(Xyzp[i]._pos.x(),Xyzp[i]._pos.y(), weights[i]);

    if (_debug > 5) {
      printf("[HelixFitHack::refineHelixParameters] %3i %10.3f %10.3f %10.5f %10.5f %10.5f %10.5f %12.5e\n",i,x,y,dx,dy,costh,sinth2, e2);
    }

  NEXT_POINT: ;
  }

  //now perform some clean up if needed
  ::LsqSums4 sxy;
  int        pointsRemoved(0);
  int        iworst;
  double     wtWorst;
  double     chi2, chi2_min;
  int        chi2Updated;
  chi2 = Trk._sxyw.chi2DofCircle();
  if (chi2 <= _chi2xyMax) {
    success = 0;
    goto F_END;
  }

  x0 = Trk._sxyw.x0();
  y0 = Trk._sxyw.y0();

  if (_debug > 5) {
    printf("[HelixFitHack::refineHelixParameters] x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f \n",
	   Trk._sxyw.x0(), Trk._sxyw.y0(), Trk._sxyw.radius(), Trk._sxyw.chi2DofCircle());
  }

  NEXT_ITERATION:;
  chi2_min    = 1e6;
  iworst      = -1;
  wtWorst     = -1;
  chi2Updated = 0;
  
  for (int i=seedIndex; i<np; i++) {
    if (Xyzp[i].isOutlier())           goto NEXT_P;

    // avoid the use of hit rejected by the helix search
    if ( indexVec[i] < 1)  goto NEXT_P;

    sxy.init(Trk._sxyw);

    x  = Xyzp[i]._pos.x();
    y  = Xyzp[i]._pos.y();
//     dx = x-x0;
//     dy = y-y0;

//     costh  = (dx*Xyzp[i]._sdir.x()+dy*Xyzp[i]._sdir.y())/sqrt(dx*dx+dy*dy);
//     sinth2 = 1-costh*costh;

    //    e2    = ew*ew*sinth2+rs*rs*costh*costh;
    wt    = weights[i];//1./e2;
    sxy.removePoint(x, y, wt);
    
    chi2  = sxy.chi2DofCircle();

    if (chi2 < chi2_min){
      chi2_min    = chi2;
      chi2Updated = 1;
      iworst      = i;
      wtWorst     = wt;
    }
  NEXT_P: ;
  }

  if (iworst >= 0){
    x  = Xyzp[iworst]._pos.x();
    y  = Xyzp[iworst]._pos.y();
    //remove point from the track    
    Trk._sxyw.removePoint(x, y, wtWorst);
    if (_debug > 5) {
      printf("[HelixFitHack::refineHelixParameters]  x0 = %8.3f y0 = %8.3f radius = %8.3f  chi2 = %5.5f chi2Maxxy = %5.5f index point removed = %i\n",
	     Trk._sxyw.x0(), Trk._sxyw.y0(), Trk._sxyw.radius(), Trk._sxyw.chi2DofCircle(), _chi2xyMax, iworst);
      
    }
    x0 = Trk._sxyw.x0();
    y0 = Trk._sxyw.y0();

    //mark point as outlier
    indexVec[iworst] = 0;
    ++pointsRemoved;
  }

  if ( (chi2_min    >= _chi2xyMax) && 
       (chi2Updated == 1)              ){
    //-----------------------------------------------------------------------------
    // still bad chi2, repeat the cleanup cycle
    //-----------------------------------------------------------------------------
    if (Trk._sxyw.qn() > 10.) {
      goto NEXT_ITERATION;
    }
  }else if ( chi2_min < _chi2xyMax){
    success = 1;
  }
  
  F_END:;
  if (_debug > 5) {
    printf("[HelixFitHack::refineHelixParameters] success = %i\n", success);
    printf("[HelixFitHack::refineHelixParameters] points removed = %i chi2 = %5.5f\n", 
	   pointsRemoved, chi2_min);
  }
//-----------------------------------------------------------------------------
// update circle parameters
//-----------------------------------------------------------------------------
  Trk._cw.setX(Trk._sxyw.x0());
  Trk._cw.setY(Trk._sxyw.y0());
  Trk._rw    = Trk._sxyw.radius();
  Trk._chi2w = Trk._sxyw.chi2DofCircle();

  if (seedIndex ==0){
    THackData* hack;
    hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
    hack->fData[14] = Trk._rw ;
    hack->fData[15] = Trk._chi2w ;
  }
//-----------------------------------------------------------------------------
// and repeat fit in phi-Z
//-----------------------------------------------------------------------------
  
//  findDfDz(Xyzp,Trk,0,_indicesTrkCandidate);

  return success;
}

//-----------------------------------------------------------------------------
  void HelixFitHack::findTrack(XYZPHackVector&     xyzp               , 
			       int                 seedIndex          , 
			       double&             chi2               , 
			       int&                countGoodPoints    ,
			       HelixFitHackResult& mytrk              ,
			       int&                mode               , 
			       bool                useDefaultDfDz     ,
			       int                 useMPVdfdz          ){
//23014-11-08 gianipez added initialization of input/output paramters
    countGoodPoints = 0;
    mode            = 0;

//Legend of the points used for performing the pattern recongition
// p0 is the position of the center of the helix
// p1 is the center of the capture target
// p2 is the point in the position seedIndex on the vector xyzp
// p3 is the postion of the EMC cluster

    Hep3Vector p0, p1, p2, p3; 
    Hep3Vector shPos, hePos;

    p1 = Hep3Vector(0., 0., 0.);   //5971. - 10200.
    p2 = xyzp[seedIndex]._pos;
    p3 = Hep3Vector(fTimePeak->ClusterX(),fTimePeak->ClusterY(), fTimePeak->ClusterZ());

    double radius,phi0,tanLambda;
    
    //----------------------------------------------------------------------//
    //index for setting the sostitute-point for the target center
    double dx,dy,phi,dx2, dy2;
    int np = xyzp.size();
    int  mode0GoodPoints(0), mode1GoodPoints(0); //mode0GoodPoints is the numebr of points belonging to a trajectory
                                                 //when dfdz is not re-calculated using the function calculateDfDz()
                                                 //mode1GoodPoints is the numebr of point belonging to a trajectory 
                                                 //when dfdz is re-computed using calculateDfDz()
    int  rescuedPoints(0);

    // initialize a vector for storing indexes of starawhit candidate for dfdz 
    // recalculation
    // createa temporary array for storing indices of the good point which belong to a track candidate

    int candidateList[np];
    int markIndexList[np];

    double dzList[np];
    double distList[np];

    for(int i=0; i<np; ++i){
      candidateList[i] = -9999;
      markIndexList[i] = -9999;
      dzList[i]        = -9999;
      distList[i]      = -9999;
    }
						
    markIndexList[seedIndex] = 1;
  
    //---------------------------------------------------------------------
    //define constrains on the z coordinate of the strawhit candidate for re-calculating dfdz
    //if the candidate and the seeding strawhit are too close, the dfdz calculated could be affected 
    //by several effects which lead on the wrong estimation.
    //We are asking that the candidate straw must be at a distance along the z axes greateer than tollMin
    // and less than tollMax. These parameters still need to be optimized
    double tollMin(100.), tollMax(500.);


    //2014-03-10 gianipez changed the values of the following tollerance
    // for X and Y distance to avoid delta electron for corrupting the pattern-reco
    //     double tollXYdist(100.);//require a dist in the transverse plane > 10 cm
    

    int goodPoint(-999); // index of the strawhit candidate for dfdz and helix paramters recalculation
    //2014-01-29 gianipez added the followign line

    //few paramterss used for calculating strawhit position residuals
    double weight(0.0);
    double deltaZ(0.), deltaZfromSeed(0.), deltaX(0.), deltaY(0.);
    double distXY(0.0);
    double dist(0.0), dist2(0.0); //help parameter for storing strawhit position residual

  //   double distGoodPoint(_distPatRec);             //threshold on the squared distance of a strawhit position with 
//                                                    //respect to the predicted position


    double z0,z1,phi_0,phi_1;
//----------------------------------------------------------------------//
// 2014-11-05 gianipez set dfdz equal to the most probable value for CE //
//----------------------------------------------------------------------//
    double dfdz = _mpDfDz;//tanLambda/radius;
    //    double chi2min(1e10);
    double dfdz_end, phi0_end, /*phi0_seed,*/ radius_end;
    //    double x0_end, y0_end;
    //    int    mode_end(-1);
    

    //Two flags are needed:
    bool removeTarget(true);                 //avoid the recalculation of dfdz and helix parameters in case when
                                             //others strawhit candidates are found 

    bool isStored(false);                    //help parameter for indicating if a strahit hasbeen already used or not 
                                             //for dfdz and helix recalculation

    chi2 = 0.0; //initialize to 0 the chi-square value

    
//----------------------------------------------------------------------
// calculate helix paramters using the center of the stopping target,
// the EMC cluster which seeded the CalTimePeak and the seeding strawhit.
///The z coordinate of the target center is set to 0 because in the formula
//inside calculateTrackParameters(...) is not used its z coordinate
    p1 = Hep3Vector(0., 0., 0.);
    calculateTrackParameters(p0,radius,phi0,tanLambda,
			     p1,p2,p3,
			     mytrk,
			     false);
    //    phi0_seed = phi0;

    //      dfdz = tanLambda/radius;
    //----------------------------------------------------------------------//
    // 2014-11-05 gianipez set dfdz equal to the most probable value for CE //
    //----------------------------------------------------------------------//
    if (useMPVdfdz ==1 ){
      dfdz = _hdfdz;//_mpDfDz;
    }

    fLastIndex = -9999;
    //--------------------------------------------------//

    ::LsqSums4 sxy;
    ::LsqSums4 srphi;

    sxy.addPoint(p2.x(), p2.y());         //seeding strawhit
    sxy.addPoint(p3.x(), p3.y());         //EMC cluster position
    sxy.addPoint(0., 0., 0.1);                 //Target center in the transverse plane


    TString banner="HelixFitHack::findTrack";
    if (useMPVdfdz ==1) banner += "-useMPVdfdz";

    double dz_max;
    TString name;
    name =  banner;
    name += "-loop";
    
    int i_last = seedIndex;
    
    for (int i=seedIndex+1; i<np; i++) {
      name =  banner;
      name += "-loop";
      if (xyzp[i].isOutlier()) goto NEXT_POINT;
      weight = 1.;
//----------------------------------------------------------------------//
// 2014-12-26 Gianipez added the request that the hit has not already 
// been used by a previous search
      if (isHitUsed(i) == 1) {
	if( _debug>5){
	  //	  printf("[HelixFitHack::findTrack-loop]  XYZP-hit number = %i skipped\n", i);
	  printf("[%s]  XYZP-hit number = %i skipped\n", name.Data(), i);
	}
	goto NEXT_POINT;
      }
      shPos = xyzp[i]._pos;
    
      deltaZ = shPos.z() - xyzp[i_last]._pos.z();

      deltaZfromSeed = shPos.z() - p2.z();

      phi    = phi0 + (deltaZ)*dfdz;//tanLambda/radius;

      hePos  = Hep3Vector(p0.x() + radius*std::cos(phi),
			  p0.y() + radius*std::sin(phi),
			  shPos.z());

      dx  = hePos.x() - shPos.x();
      dx2 = dx*dx;
      dy  = hePos.y() - shPos.y();
      dy2 = dy*dy;
      
      dist2 = dx2 + dy2;
      dist  = std::sqrt(dist2);

      //calculate the distance of the strahit point to the seeding strawhit on the transverse plane
      deltaY = std::fabs(p2.y() - shPos.y());
      deltaX = std::fabs(p2.x() - shPos.x());
      distXY = std::sqrt(deltaY*deltaY + deltaX*deltaX);
   
      if( _debug>5){
	if( i==seedIndex+1) {
	  printf("[%s]  findTrack() starts with helix parameters derived from these points \n", name.Data());
	  printf("[%s]   point  type      X         Y         Z       xyzp-index \n", name.Data());
	  printf("[%s] ----------------------------------------------------------\n", name.Data());
	  printf("[%s]    seeding      %10.3f   %10.3f   %10.3f   %8i \n", name.Data(),p2.x(), p2.y(), p2.z(), 
		 seedIndex);
	  printf("[%s]   candidate     %10.3f   %10.3f   %10.3f   %8i \n", name.Data(),p1.x(), p1.y(), p1.z(), 
		 fLastIndex);
	  printf("[%s]  emc cluster    %10.3f   %10.3f   %10.3f   %8i \n", name.Data(),p3.x(), p3.y(), p3.z(), 
		 -1);
	}
	
	printf("[%s] X0 = %10.3f Y0 = %10.3f r = %10.3f dfdz = %5.5f \n", name.Data(), 
	       p0.x(), p0.y(), radius, dfdz);
	printf("[%s]   measured      %10.3f %10.3f %10.3f %8i \n", name.Data(), 
	       shPos.x(), shPos.y(), shPos.z(), i);
	printf("[%s]  predicted      %10.3f %10.3f %10.3f %8i \n", name.Data(), 
	       hePos.x(), hePos.y(), hePos.z(), i);
	printf("[%s] dist-from-prediction = %10.3f  dist-from-seedXY = %5.3f dz-from-seed = %5.3f\n", name.Data(), 
	       dist, distXY, deltaZ);
	printf("[%s] ---------------------------------------------------------\n", name.Data());

      }
      
      isStored = false;
      dz_max = _distPatRec + _dfdzErr*deltaZ;
      if ( dist <= dz_max ){
	++countGoodPoints;

	// 2014-11-12 gianipez:
	//adjust the helix-center coordinates and the radius
	// if the dfdz value has already been evaluated
	if (mode == 1 ){
	  p0.setX( sxy.x0());
	  p0.setY( sxy.y0());
	  radius  = sxy.radius();
	}

	//store index of last good point
	i_last = i; 

	//2015-01-27 G. Pezzu and P. Murat
	//phi0   = phi;
	phi0   = CLHEP::Hep3Vector(shPos - p0).phi();

	//add point to the helixfithack result objet
	sxy.addPoint(shPos.x(),shPos.y(), weight);

	//store the index of the good point found
	markIndexList[i] = 1;

	//store distance along z-axis from the last point found relying on thr helix
	dzList[i]        = deltaZ;

	//store distance from predition
	distList[i]      = dist;



	if (mode == 0){
	  ++mode0GoodPoints;	
	}else if ( (mode == 1) &&
		   (i<= fLastIndex)){
	  ++mode1GoodPoints;
	}

	for(int j=0; j<np; ++j){
	  if(candidateList[j] == i)
	    isStored = true;
	}
	  
	if(!isStored){
	  if( ( deltaZfromSeed > tollMin ) &&
	      ( deltaZfromSeed < tollMax )){	    
	    if (removeTarget) goodPoint = i;
	  }
	}
      }else {
	markIndexList[i] = 0;
	distList[i]      = 0;
	dzList[i]        = 0;
      }

      // 2014-04-23     gianipez fixed a bug
      if ( countGoodPoints >= 2 &&
	   removeTarget         &&
	   (goodPoint >=0)      &&
	   (goodPoint != fLastIndex) ) {

//recalculate helix parameters using the strawhit candidate "goodPoint"
	p1 = xyzp[goodPoint]._pos;

	p0.setX( sxy.x0());
	p0.setY( sxy.y0());
	radius  = sxy.radius();

//now calculate more accuratelly the value of dfdz using just the two strawhit positions
	z0 = p2.z();//z coordinate of the seeding point
	z1 = p1.z();
	phi_0 = CLHEP::Hep3Vector(p2 - p0).phi();
	phi_1 = CLHEP::Hep3Vector(p1 - p0).phi();

//2015-01-14 G. Pezzullo added the following condition because in case 
// we have a MPV for dfdz from the procedure findDfDZ we want just to use it
	if (useMPVdfdz == 0){
	  calculateDfDz(phi_0, phi_1, z0, z1, dfdz);
	} else if (useMPVdfdz ==1 ){
	  dfdz = _hdfdz;//_mpDfDz;
	}

	name = banner;
	name += "2strawhitsHelixDef";
	if (_debug > 5) {
	  printf("[%s] strawhit type     X        Y        Z     index\n", name.Data());
	  printf("[%s] ----------------------------------------------------\n", name.Data());
	  printf("[%s]    seeding     %5.3f  %5.3f  %5.3f   %i  \n", name.Data(),p2.x(), p2.y(), p2.z(), seedIndex);
	  printf("[%s]   candidate    %5.3f  %5.3f  %5.3f   %i  \n", name.Data(),p1.x(), p1.y(), p1.z(), goodPoint);
	  printf("[%s] x0 = %5.3f y0 = %5.3f radius = %5.3f dfdz = %5.6f chi2 = %5.3f \n", name.Data(),
		 p0.x(), p0.y(), radius, dfdz , sxy.chi2DofCircle());
	
	}

	//what to do if dfdz s negative?
	//	if( dfdz<0.0) {
	if ((dfdz > _maxDfDz) || (dfdz < _minDfDz)) {

	  for(int j=0; j<np; ++j){
	    if(candidateList[j]<0){
	      candidateList[j] = goodPoint;
	      break;
	    }
	  }
	  if (_debug > 5) printf("[%s]  dfdz = %5.5f not in range limits. Continue the search\n", 
				 name.Data(),
				 dfdz);
	  p1 = Hep3Vector(0., 0., 0.);
//----------------------------------------------------------------------//
// 2014-11-05 gianipez set dfdz equal to the most probable value for CE //
//----------------------------------------------------------------------//
	  dfdz = _mpDfDz;
				// 	    dfdz = (tanLambda/radius)*double(j);
	}else{
	  removeTarget = false;
	  mode         = 1;
	  fLastIndex   = goodPoint;
	}
	//------------------------------------------------------------//
      }
    NEXT_POINT:;
    }
    name = banner;
    name += "-results";
   

    chi2 = sxy.chi2DofCircle();

      //2014-04-23 gianipez added the following line
    //    chi2min    = sxy.chi2DofCircle();//chi2;
    //    x0_end     = p0.x();
    //    y0_end     = p0.y();

//2015-01-22 G. Pezzullo and P. Murat; update the dfdz value usingall hits
    HelixFitHackResult tmp1HelFitRes( mytrk);
    HelixFitHackResult tmp2HelFitRes( mytrk);

//     mytrk._center.set(p0.x(), p0.y(), 0.0);
//     mytrk._radius = radius;
//     mytrk._sxy    = sxy;

  

    
    //2015-01-27 G. Pezzu and P. Murat
    // initialize only the xy part, 
    // we do not need the z-phi part
    tmp1HelFitRes._sxy.init(sxy);
    tmp1HelFitRes._radius = sxy.radius();
    tmp1HelFitRes._center.set(sxy.x0(), sxy.y0(), 0.0);

    radius_end = sxy.radius();

    tmp2HelFitRes._center.set(p0.x(), p0.y(), 0.0);
    tmp2HelFitRes._radius = radius;
    
    //    if ( findXY_new(xyzp, tmp1HelFitRes, seedIndex, markIndexList)){
    if ( refineHelixParameters(xyzp, tmp1HelFitRes, seedIndex, markIndexList) >=0){
//       tmp2HelFitRes._center.set(tmp1HelFitRes._center.x(), tmp1HelFitRes._center.y(), 0.0);
//       tmp2HelFitRes._radius = tmp1HelFitRes._radius;
//       radius_end            = tmp1HelFitRes._radius;
//       sxy.init(tmp1HelFitRes._sxy);
      tmp2HelFitRes._center.set(tmp1HelFitRes._cw.x(), tmp1HelFitRes._cw.y(), 0.0);
      tmp2HelFitRes._radius = tmp1HelFitRes._rw;
      radius_end            = tmp1HelFitRes._rw;
      sxy.init(tmp1HelFitRes._sxyw);

      //update the chi2 value
      chi2 = sxy.chi2DofCircle();


      countGoodPoints = 0;
      for (int i=seedIndex; i<np; ++i){
	if ( markIndexList[i] >0 ) 
	  ++countGoodPoints;
      }
    }
    
    //     tmpHelFitRes._center.set(p0.x(), p0.y(), 0.0);
//     tmpHelFitRes._radius = radius;

    findDfDz(xyzp, tmp2HelFitRes, seedIndex, markIndexList);
    tmp2HelFitRes._dfdz = _hdfdz;

    //2015-01-23 G. Pezzu and P. Murat
    // findZ returns negative value in case it fails, so take it into account!
    // use the previous value for dfdz and phi0
    if ( findZ(xyzp, tmp2HelFitRes, seedIndex, markIndexList)){
      dfdz_end   = tmp2HelFitRes._dfdz;
      phi0_end   = tmp2HelFitRes._fz0;
      //      _hdfdz     = dfdz_end;
      srphi.init(tmp2HelFitRes._srphi);
      
      //countGoodPoints       = tmp1HelFitRes._srphi.qn();//FIX ME
      countGoodPoints = 0;
      for (int i=seedIndex; i<np; ++i){
	if ( markIndexList[i] >0 ) 
	  ++countGoodPoints;
      }
       
    }else {
      dfdz_end = _hdfdz;
      phi0_end = _hphi0;
    }
    //    phi0_end   = phi0_seed;


     if (_debug > 5) {
      printf("[%s] strawhit type     X        Y        Z     index\n", name.Data());
      printf("[%s] ----------------------------------------------------\n", name.Data());
      printf("[%s]    seeding     %5.3f  %5.3f  %5.3f   %i  \n", name.Data(),p2.x(), p2.y(), p2.z(), seedIndex);
      printf("[%s]   candidate    %5.3f  %5.3f  %5.3f   %i  \n", name.Data(),p1.x(), p1.y(), p1.z(), goodPoint);
      printf("[%s]  emc cluster   %5.3f  %5.3f  %5.3f \n", name.Data(),p3.x(), p3.y(), p3.z());
      printf("[%s] x0 = %5.3f y0 = %5.3f radius = %5.3f dfdz = %5.6f chi2 = %5.3f \n", name.Data(),
	     sxy.x0(), sxy.y0(), radius_end, dfdz_end , sxy.chi2DofCircle());
      printf("[%s] countGoodPoints = %i\n", name.Data(), countGoodPoints);
    }





     //    mode_end   = mode;
    if(mode1GoodPoints>0){
      rescuedPoints = mode1GoodPoints - mode0GoodPoints ;
    } else{
      rescuedPoints = -1;
    }

    //----------------------------------------------------------------------

    if (countGoodPoints > 2) {
      //2014-01-29 gianipez added the following line
      chi2       = sxy.chi2DofCircle();
      name = banner;
      name += "EndFindTrack";
      if (_debug > 5) {
	printf("[%s] chi2 = %5.3f nGoodPoints = %d dfdz = %5.5f mode = %i\n",
	       name.Data(),
	       sxy.chi2DofCircle(),
	       countGoodPoints,
	       dfdz_end,
	       mode);
      }
    }
    
    
    if ( (mode == 1) ||
	 useDefaultDfDz ){
      if ( (countGoodPoints > _goodPointsTrkCandidate) ||
	   ( (countGoodPoints == _goodPointsTrkCandidate) && 
	     (chi2             < _chi2TrkCandidate      ) ) ){
	//      update trackcandidate informations
	fSeedIndex      = seedIndex;
	if (mode == 1) {
	  fUseDefaultDfDz = 1;
	  fCandidateIndex = fLastIndex;
	}else if (useDefaultDfDz){
	  fCandidateIndex = -9999;
	}
//----------------------------------------------------------------------
// 2015 - 01 - 17 G. Pezzu: remove the target center from sxy 
// in order to evaluate more accuratelly the helix parameters
//----------------------------------------------------------------------
//	sxy.removePoint(0., 0., 0.1);                 //Target center in the transverse plane

	_x0     = sxy.x0();	// p0.x();
	_y0     = sxy.y0();	// p0.y();
	_phi0   = phi0_end;

	_radius = radius_end;//sxy.radius();	// radius_end;
	_dfdz   = dfdz_end;

	_goodPointsTrkCandidate = countGoodPoints;
	_chi2TrkCandidate       = chi2;
      
	
	//reset the vector holding the informations aboit:
	// -> hit belonging to the track candidate
	// -> distance in the X-Y plane from the prediction
	// -> distance from the seeding hit along the z-axes
	for(int i=0; i<400; ++i){
	  _indicesTrkCandidate[i] = -9999;
	  _distTrkCandidate[i]    = -9999;
	  _dzTrkCandidate[i]      = -9999;
	}
	
	for (int i=seedIndex; i<np; ++i){
	  _indicesTrkCandidate[i] = markIndexList[i];
	  _distTrkCandidate[i]    = distList[i];
	  _dzTrkCandidate[i]      = dzList[i];

	}

	mytrk._center.set(_x0, _y0, 0.0);
	mytrk._radius = _radius;
	//now calculate the phi coordinate at z = 0
	mytrk._fz0    = _phi0;// + _dfdz*(0.0 - p2.z()) ;
	mytrk._dfdz   = _dfdz;
	mytrk._sxy.init(sxy);
	mytrk._srphi.init(srphi);
	THackData* hack;
	hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
	int loopId(0);
	if (useDefaultDfDz == 0){
	  if (useMPVdfdz){
	    loopId = 2;
	  }else {
	    loopId = 0;
	  }
	}else {
	  loopId = 1;
	}
	hack->fData[4]  = loopId;
	hack->fData[5]  = _radius;
	hack->fData[6]  = phi0_end;
	hack->fData[7]  = dfdz_end*_radius;
	hack->fData[8]  = dfdz_end;//_dfdz;//_mpDfDz; 
	hack->fData[9]  = rescuedPoints;
	double dz = p1.z() - p2.z();
	hack->fData[10] = mode == 1 ? dz : -1.;
	hack->fData[11] = countGoodPoints;
	hack->fData[12] = sxy.chi2DofCircle();//sxy.chi2DofCircle();
	hack->fData[13] = srphi.chi2rphiDofCircle();

	int j=0;
	for (int i=seedIndex; i<np; ++i){
	  if (_indicesTrkCandidate[i] != 1) continue;
	  hack->fDist[j] = _distTrkCandidate[i];
	  hack->fDz[j]   = _dzTrkCandidate[i];
	  ++j;
	}	
      }
    }
   
  }
  

//-----------------------------------------------------------------------------
  void HelixFitHack::calculateTrackParameters(Hep3Vector&         p0, 
					      double&             radius,
					      double&             phi0, 
					      double&             tanLambda,
					      Hep3Vector          p1, 
					      Hep3Vector          p2,
					      Hep3Vector          p3,
					      HelixFitHackResult& mytrk,
					      bool                cleanPattern) 
  {
    
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
    if (_debug > 5) {
      printf("[HelixFitHack:calculateTrackParameters] phi0 from X = %5.3f phi0 from Y = %5.3f phi0 from tan = %5.3f p1.z = %10.3f p2.z = %10.3f p3.z = %10.3f\n",
	     phi0, phi0_fromY, phi0_fromXY,p1.z(),p2.z(),p3.z());
    }
    phi0 = phi0_fromXY;

    double deltaPhi_0 = std::atan2(deltaY,deltaX) - phi0;
    double deltaPhi_1 = std::atan(deltaY/deltaX) - phi0;
    if (_debug > 5) {
      printf("[HelixFitHack:calculateTrackParameters] deltaPhi_0 = %5.5f deltaPhi_1 = %5.5f\n",
	     deltaPhi_0,
	     deltaPhi_1);
    }
    
    if(deltaPhi_0 < 0.) 
      deltaPhi_0 += 2.*M_PI;

    tanLambda = (radius/deltaZ)*deltaPhi_0;//(radius/deltaZ)*std::fabs(deltaPhi_0);

    if (_debug > 5) {
      printf("[HelixFitHack:calculateTrackParameters] dfdz = %5.8f \n", 
	     tanLambda/radius);
    }
    
  }
  
  void  HelixFitHack::calculateDfDz(double &phi0, double &phi1, 
				    double &z0,   double &z1,
				    double &dfdz){
    double deltaPhi = phi1 - phi0;
    while(deltaPhi < -M_PI) deltaPhi += 2.*M_PI;
    while(deltaPhi >  M_PI) deltaPhi -= 2.*M_PI;
    dfdz = deltaPhi/(z1 - z0);
    if(dfdz>0.0) 
      if (_debug > 5) {
	printf("[HeliFitHack::calculateDfDZ] df = %5.3f dz = %5.3f dfdz = %5.5f\n",
	       deltaPhi, (z1 - z0), dfdz);
      }
  }
  
}



//
// Object to perform helix fit to straw hits
//
// $Id: HelixFit.cc,v 1.5 2013/01/26 18:17:13 brownd Exp $
// $Author: brownd $ 
// $Date: 2013/01/26 18:17:13 $
//
//
// the following has to come before other BaBar includes
#include "BaBar/BaBar.hh"
#include "TrkPatRec/inc/HelixFit.hh"
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

namespace mu2e 
{
// comparison functor for ordering points
  struct radcomp : public std::binary_function<VALERR, VALERR, bool> {
    bool operator()(VALERR const& r1, VALERR const& r2) { return r1._val < r2._val; }
  };

  // comparison functor for sorting by z
  struct zcomp : public std::binary_function<XYZP*,XYZP*,bool> {
    bool operator()(XYZP* const& p1, XYZP* const& p2) { return p1->_pos.z() < p2->_pos.z(); }
  };
 
  double StereoXYZP::_sfactor(1.0);

  StereoXYZP::StereoXYZP(HitXYZP& h1,HitXYZP& h2) :  _h1(&h1),_h2(&h2){
    double invdet = 1.0/(-_h1->_wdir.x()*_h2->_wdir.y() + _h1->_wdir.y()*_h2->_wdir.x());
    Hep3Vector dp = _h2->_pos - _h1->_pos;
    _d1 = invdet*(-_h2->_wdir.y()*dp.x() + _h2->_wdir.x()*dp.y());
    _d2 = invdet*(-_h1->_wdir.y()*dp.x() + _h1->_wdir.x()*dp.y());
    _pos = 0.5*(_h1->_pos + _d1*_h1->_wdir +_h2->_pos + _d2*_h2->_wdir);
    _phi = _pos.phi();
    _dz = fabs(_h2->_pos.z()-_h1->_pos.z());
    // sign the wire direction to point radially outwards
    CLHEP::Hep3Vector w1 = copysign(1.0,_h1->_wdir.dot(_pos))*_h1->_wdir;
    CLHEP::Hep3Vector w2 = copysign(1.0,_h2->_wdir.dot(_pos))*_h2->_wdir;
    // compute the radial and azimuthal projections, assuming a
    // worst-case pitch of 45 degrees
    double cosf = w1.dot(w2);
    static const double invsqrt12(1.0/sqrt(12));
    _rerr = _sfactor*invsqrt12*_dz*sqrt(0.5*(1.0+cosf));
    _ferr = _sfactor*invsqrt12*_dz*sqrt(0.5*(1.0-cosf));
    setUse(true,good);
    setUse(true,stereopoint);
  }

  void
  StereoXYZP::rinfo(CLHEP::Hep3Vector const& center,VALERR& rad) const {
    CLHEP::Hep3Vector rvec = (_pos - center).perpPart();
    rad._val = rvec.mag();
    double alpha = rvec.angle(_pos.perpPart());
    rad._err = _rerr*fabs(cos(alpha))+_ferr*fabs(sin(alpha));
  }

  void
  StereoXYZP::finfo(CLHEP::Hep3Vector const& center,VALERR& phi) const {
    CLHEP::Hep3Vector rvec = (_pos-center).perpPart();
    phi._val = rvec.phi();
    double alpha = rvec.angle(_pos.perpPart());
    phi._err = (_ferr*fabs(cos(alpha)+_rerr*fabs(sin(alpha))))/rvec.perp();
  }

  void
  StereoXYZP::setUse(bool use,useBit ibit) {
    XYZP::setUse(use,ibit);
    if(ibit==good){
      _h1->setUse(use,stereohit);
      _h2->setUse(use,stereohit);
    }
  }
  bool 
  StereoXYZP::use(useBit ibit) const {
   if(ibit == good)
    return XYZP::use() && !(XYZP::use(outlier));
  else
    return XYZP::use(ibit);// require this NOT be an outlier
  }

  void
  HitXYZP::rinfo(CLHEP::Hep3Vector const& center,VALERR& rad) const {
//    static const double onethird(1.0/3.0);
//    static const double invsqrt12(1./sqrt(12.0));
// average the 1-sigma radii to account for non-linear errors
    double rvec = CLHEP::Hep3Vector(_pos - center).perp();
    double rvec1 = CLHEP::Hep3Vector(_pos +_werr*_wdir - center).perp();
    double rvec2 = CLHEP::Hep3Vector(_pos -_werr*_wdir - center).perp();
//    rad._val = onethird*(rvec+rvec1+rvec2);
    rad._val = rvec;
    rad._err = std::max(std::max(fabs(rvec1-rvec),fabs(rvec2-rvec)),_serr);
  }

  void
  HitXYZP::finfo(CLHEP::Hep3Vector const& center,VALERR& phi) const {
//    static const double onethird(1.0/3.0);
//    static const double invsqrt12(1./sqrt(12.0));
// average the 1-sigma radii to account for non-linear errors
    double rvec = CLHEP::Hep3Vector(_pos - center).perp();
    double phi0 = CLHEP::Hep3Vector(_pos - center).phi();
    double phi1 = CLHEP::Hep3Vector(_pos +_werr*_wdir - center).phi();
    double phi2 = CLHEP::Hep3Vector(_pos -_werr*_wdir - center).phi();
//    rad._val = onethird*(rvec+rvec1+rvec2);
    phi._val = phi0;
    phi._err = std::max(std::max(fabs(phi1-phi0),fabs(phi2-phi0)),_serr/rvec);
  }

  bool 
  HitXYZP::use(useBit ibit) const {
   if(ibit == good)
    return XYZP::use() && !XYZP::use(outlier) && !XYZP::use(stereohit);
  else
    return XYZP::use(ibit);// require this NOT be an outlier
  }



  HelixFitResult& 
  HelixFitResult::operator =(HelixFitResult const& other) {
    if(this != &other){
      _tdef = other._tdef;
      _fit = other._fit;
      _center = other._center;
      _radius = other._radius;
      _dfdz = other._dfdz;
      _fz0 = other._fz0;
    }
    return *this;
  }
  
  XYZPVector::XYZPVector(std::vector<HitXYZP> const& hits) : _hxyzp(hits) {
    _xyzp.reserve(_hxyzp.size());
    for(size_t ih=0;ih<_hxyzp.size();++ih)
      _xyzp.push_back(&_hxyzp[ih]);
  }

  XYZPVector::XYZPVector(std::vector<HitXYZP> const& hits,std::vector<StereoXYZP> const& shits) : _hxyzp(hits), _sxyzp(shits) {
    _xyzp.reserve(_hxyzp.size()+_sxyzp.size());
    for(size_t ih=0;ih<_hxyzp.size();++ih)
      _xyzp.push_back(&_hxyzp[ih]);
    for(size_t is=0;is<_sxyzp.size();++is)
      _xyzp.push_back(&_sxyzp[is]);
  }

  void
  HelixFit::helixParams(HelixFitResult const& helix,CLHEP::HepVector& pvec,CLHEP::HepVector& perr) const {
    TrkDef const& mytrk = helix._tdef;
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

  HelixFit::HelixFit(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _mindelta(pset.get<double>("minDelta",5000.0)),
    _minnhit(pset.get<unsigned>("minNHit",10)),
    _minnstereo(pset.get<unsigned>("minnStereo",0)),
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
    _rbias(pset.get<double>("radialBias",0.0)),
    _sfac(pset.get<double>("strawSizeFactor",1.0)),
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
    _sfactor(pset.get<double>("stereoFactor",1.0)),
    _forcep(pset.get<bool>("forceP",false)),
    _xyweights(pset.get<bool>("xyWeights",true)),
    _zweights(pset.get<bool>("zWeights",false)),
    _filter(pset.get<bool>("filter",true)),
    _plotall(pset.get<bool>("plotall",false)),
    _usetarget(pset.get<bool>("usetarget",true)),
    _allstereo(pset.get<bool>("allstereo",true)),
    _bz(0.0),_sdist(0),_spull(0)
  {
    StereoXYZP::_sfactor = _sfactor;
  }

  HelixFit::~HelixFit()
  {}

  double
  HelixFit::bz() const {
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
  HelixFit::findHelix(HelixFitResult& myhel) {
    TrkDef const& mytrk = myhel._tdef;
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
    std::vector<HitXYZP> hxyzp;
    fillHitXYZP(mytrk,hxyzp);
    // find the stereo pairs for these
    std::vector<StereoXYZP> sxyzp;
    findStereoPairs(hxyzp,sxyzp);
    XYZPVector xyzp(hxyzp,sxyzp);
// call down
    bool retval = findHelix(xyzp,myhel);
    if((retval || _plotall) && _diag>1){
      // fill graphs for display if requested
      plotXY(mytrk,xyzp._hxyzp,xyzp._sxyzp,myhel);
      plotZ(mytrk,xyzp._hxyzp,xyzp._sxyzp,myhel);
    }
    return retval;
  }

  bool
  HelixFit::findHelix(XYZPVector& xyzp,HelixFitResult& myhel) {
    bool retval(false);
// filter by geometry
    if(_filter)filterDist(xyzp._xyzp);
    if(xyzp._hxyzp.size() >= _minnhit && xyzp._sxyzp.size() >= _minnstereo){
      // initialize the circle parameters
      if(initCircle(xyzp._xyzp,myhel)){
	// solve for the circle parameters
	retval = findXY(xyzp._xyzp,myhel);
	// extend those into the z direction
	if(retval){
	  retval = findZ(xyzp._xyzp,myhel);
	  // set the success
	  if(retval){
	    myhel._fit = TrkErrCode(TrkErrCode::succeed);
	  } else
	    myhel._fit = TrkErrCode(TrkErrCode::fail,4); // phi-z reconstruction failure
	} else
	  myhel._fit = TrkErrCode(TrkErrCode::fail,3); // xy reconstruction failure
      } else
	myhel._fit = TrkErrCode(TrkErrCode::fail,2); // initialization failure
    } else
      myhel._fit = TrkErrCode(TrkErrCode::fail,1); // insufficient hits
    return retval;
  }
  
  bool
  HelixFit::findXY(std::vector<XYZP*>& xyzp,HelixFitResult& myhel) {
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

  void
  HelixFit::plotXY(TrkDef const& mytrk, std::vector<HitXYZP> const& hxyzp,
    std::vector<StereoXYZP> const& sxyzp,
    HelixFitResult const& myhel) const {
    unsigned igraph = 10*mytrk.eventId()+mytrk.trackId();
    art::ServiceHandle<art::TFileService> tfs;
    if(_sdist == 0){
      _sdist = tfs->make<TH1F>("sdist","Selected Stereo distance",100,-500,500);
      _spull = tfs->make<TH1F>("spull","Selected Stereo pull",100,-10,10);
    }
    //      TGraph* graph = tfs->make<TGraph>(hxyzp.size());
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
    for(unsigned ihp=0;ihp<hxyzp.size();++ihp){
      if(hxyzp[ihp].use())
	g->Fill(hxyzp[ihp]._pos.x()-myhel._center.x(),hxyzp[ihp]._pos.y()-myhel._center.y());
      else
	b->Fill(hxyzp[ihp]._pos.x()-myhel._center.x(),hxyzp[ihp]._pos.y()-myhel._center.y());
    }
    for(size_t isp=0;isp<sxyzp.size();++isp){
      _sdist->Fill(sxyzp[isp]._d1);
      _sdist->Fill(sxyzp[isp]._d2);
      _spull->Fill(sxyzp[isp]._d1/sxyzp[isp]._h1->_werr);
      _spull->Fill(sxyzp[isp]._d2/sxyzp[isp]._h2->_werr);
      if(sxyzp[isp].use())
	s->Fill(sxyzp[isp]._pos.x()-myhel._center.x(),sxyzp[isp]._pos.y()-myhel._center.y());
      else 
	bs->Fill(sxyzp[isp]._pos.x()-myhel._center.x(),sxyzp[isp]._pos.y()-myhel._center.y());
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
  HelixFit::findCenterAGE(std::vector<XYZP*> const& xyzp,Hep3Vector& center, double& rmed, double& age,bool useweights) {
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
  HelixFit::findZ(std::vector<XYZP*>& xyzp,HelixFitResult& myhel) {
    using namespace boost::accumulators;
// sort points by z
    std::sort(xyzp.begin(),xyzp.end(),zcomp());
// find phi information
    std::vector<FZ> finfo;
    finfo.reserve(xyzp.size());
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp]->use()){
	FZ fz;
	xyzp[ixyzp]->finfo(myhel._center,fz._phi);
	fz._z = xyzp[ixyzp]->_pos.z();
	finfo.push_back(fz);
      }
    }
// look for residual delta rays
//    art::ServiceHandle<art::TFileService> tfs;
//    static unsigned igraph(0);
//    ++igraph;
//    char gname[100];
//    snprintf(gname,100,"zphi%i",igraph);
//    TH1F* g = tfs->make<TH1F>(gname,"findZ phi",50,-1.1*pi,1.1*pi);
//    for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
//      double dphi = deltaPhi(0.0,finfo[iphi]._phi._val);
//      g->Fill(dphi);
//      if((pi-dphi)/pi<0.1)g->Fill(dphi-twopi);
//      if((pi+dphi)/pi<0.1)g->Fill(dphi+twopi);
//    }
//    int imax = g->GetMaximumBin();
//    double maxn = g->GetBinContent(imax);
//    double fmax = g->GetBinCenter(imax);
//    double frac = maxn/finfo.size();
// test
//    if(frac > _dfrac){
// mark the hits 
//    }
// find the initial slope
    if(finfo.size() > _minnhit){
      // make initial estimate of dfdz using 'nearby' pairs
      accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > accf;
//      accumulator_set<double, stats<tag::median(with_p_square_quantile) > > acctest;
      for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
	for(unsigned jphi=iphi+1; jphi < finfo.size(); ++jphi){
	  double dz = finfo[jphi]._z - finfo[iphi]._z;
	  if(dz > _minzsep && dz < _maxzsep){
	    double dphi = deltaPhi(finfo[iphi]._phi._val,finfo[jphi]._phi._val);
	    double slope = dphi/dz;
//	    if(slope > _smin && slope < _smax){ 
	      double wt = _zweights ? dz/(finfo[iphi]._phi._err+finfo[jphi]._phi._err) : 1.0;
	      accf(slope,weight=wt);
//	      acctest(slope);
//	    }
	  }
	}
      }
      double dfdz = extract_result<tag::weighted_median>(accf);
      // if the sign of dfdz disagrees, abort
      if( dfdz * _dfdzsign < 0.0)
	return false;
      // if requested, restrict the range
      if(_forcep)
	dfdz = std::max(std::min(dfdz,_smax),_smin);
      else
	if(dfdz > _smax || dfdz < _smin) return false;
//      double dfdztest = extract_result<tag::weighted_median>(acctest);
// find phi at z intercept.  Bootstrap using the mode, since phi looping
// hasn't been resolved yet
      TH1F hphi("hphi","phi value",50,-1.1*pi,1.1*pi);
      for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
      	double phiex = finfo[iphi]._z*dfdz;
	double dphi = deltaPhi(phiex,finfo[iphi]._phi._val);
	hphi.Fill(dphi);
	if((pi-dphi)/pi<0.1)hphi.Fill(dphi-twopi);
	if((pi+dphi)/pi<0.1)hphi.Fill(dphi+twopi);
      }
      double fz0 = hphi.GetBinCenter(hphi.GetMaximumBin());
// refine this using the median
//      accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > acci;
//      for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
//	double wt = _zweight ? 1.0/finfo[iphi]._phi._err : 1.0;
//	double phiex = fz0 + finfo[iphi]._z*dfdz;
//	double dphi = deltaPhi(phiex,finfo[iphi]._phi._val);
//	acci(dphi, weight=wt);
//      }
//      fz0 += extract_result<tag::weighted_median>(acci);
      if(fz0>pi)fz0 -= twopi;
      if(fz0<-pi)fz0 += twopi;
  // iterate over slope and ambiguity resolution
      bool changed(true);
      unsigned niter(0);
      while(changed && niter < _maxniter){
  // resolve phi over the full z range using the initial slope and intercept
	changed = false;
	++niter;
	for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
	  double phiex = fz0 + finfo[iphi]._z*dfdz;
	  int nloop = (int)rint((phiex-finfo[iphi]._phi._val)/twopi);
	  finfo[iphi]._phi._val += nloop*twopi;
	  changed |= nloop != 0;
	}
	// make a long-range estimate of slope
	accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > accf2;
//	accumulator_set<double, stats<tag::median(with_p_square_quantile) > > acctest2;
	for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
	  for(unsigned jphi=iphi+1; jphi < finfo.size(); ++jphi){
	    double dz = finfo[jphi]._z -finfo[iphi]._z;
	    if(dz > _minzsep){
	      double dphi = finfo[jphi]._phi._val-finfo[iphi]._phi._val;
	      double dphiex = dz*dfdz;
	      double ferr = finfo[iphi]._phi._err+finfo[jphi]._phi._err;
	      if(!_filter || fabs(dphi-dphiex) < _nsigma*ferr){
		double slope = dphi/dz;
//		if(slope > _smin && slope < _smax){ 
// limit the weight so as not to count more than 1 loop
		  double wt = _zweights ? std::min(dz,_maxzsep)/ferr : 1.0;
		  accf2(slope,weight=wt);
//		  acctest2(slope);
//		}
	      }
	    }
	  }
	}
	dfdz = extract_result<tag::weighted_median>(accf2);
//	dfdztest = extract_result<tag::weighted_median>(acctest2);
	// find phi at z intercept
	accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > acci2;
	for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
	  double phiex = fz0+finfo[iphi]._z*dfdz;
	  if(!_filter || fabs(finfo[iphi]._phi._val-phiex) < _nsigma*finfo[iphi]._phi._err){
	    double wt = _zweights ? 1.0/finfo[iphi]._phi._err : 1.0;
	    acci2(finfo[iphi]._phi._val - finfo[iphi]._z*dfdz, weight=wt);
	  }
	}
	fz0 = fmod(extract_result<tag::weighted_median>(acci2),twopi);
	if(fz0>pi)fz0 -= twopi;
	if(fz0<-pi)fz0 += twopi;
      }
      myhel._dfdz = dfdz;
      myhel._fz0 = fz0;
      // fix the phi for the hit points
      for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
	double phiex = myhel._fz0 + xyzp[ixyzp]->_pos.z()*myhel._dfdz;
	FZ fz;
	xyzp[ixyzp]->finfo(myhel._center,fz._phi);
	int nloop = (int)rint((phiex - fz._phi._val)/twopi);
	xyzp[ixyzp]->_phi = fz._phi._val + nloop*twopi;
	if(_filter && fabs(xyzp[ixyzp]->_phi-phiex)> _nsigma*fz._phi._err){
	  xyzp[ixyzp]->setUse(false);
	}
      }
      return true;
    } else
      return false;
  }

  void
  HelixFit::plotZ(TrkDef const& mytrk, std::vector<HitXYZP> const& hxyzp, std::vector<StereoXYZP> const& spair, HelixFitResult const& myhel) const {
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
    for(unsigned ih=0;ih<hxyzp.size();++ih){
      //        graph->SetPoint(ih,hxyzp[ih]._pos.x()-myhel._center.x(),hxyzp[ih]._pos.y()-myhel._center.y());
      if(hxyzp[ih].use())
	g->Fill(hxyzp[ih]._pos.z(),hxyzp[ih]._phi);
      else
	b->Fill(hxyzp[ih]._pos.z(),hxyzp[ih]._phi);
    }
    for(size_t isp=0;isp<spair.size();++isp){
      //      double phi = atan2(spair[isp]._pos.y()-myhel._center.y(),spair[isp]._pos.x()-myhel._center.x());
      //      double phiex = myhel._fz0 + spair[isp]._pos.z()*myhel._dfdz;
      //      int nloop = (int)rint((phiex-phi)/twopi);
      //      phi += nloop*twopi;
      if(spair[isp].use())
	s->Fill(spair[isp]._pos.z(),spair[isp]._phi);
      else 
	bs->Fill(spair[isp]._pos.z(),spair[isp]._phi);
    }

    TF1* line = new TF1("line","[0]+[1]*x",-1500,1500);
    line->SetParameter(0,myhel._fz0);
    line->SetParameter(1,myhel._dfdz);
    line->SetLineColor(kRed);
    TList* flist = g->GetListOfFunctions();
    flist->Add(line);
  }

  bool
  HelixFit::initCircle(std::vector<XYZP*> const& xyzp,HelixFitResult& myhel) {
    bool retval(false);
    static const double mind2 = _mindist*_mindist;
    using namespace boost::accumulators;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accx, accy, accr;
    // form all triples, and compute the circle center for unaligned hits.  I can aford to be choosy
    unsigned ntriple(0);
    unsigned nxyzp = xyzp.size();
    std::vector<CLHEP::Hep3Vector> pos;
    pos.reserve(xyzp.size());
    for(size_t ixyzp=0; ixyzp < nxyzp; ++ixyzp){
      if(xyzp[ixyzp]->use()){
	pos.push_back(xyzp[ixyzp]->_pos);
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
      // use the center to estimate the radius
      //      accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accr;
      //      for(unsigned ixyzp=0; ixyzp<nxyzp; ++ixyzp){
      //	if(xyzp[ixyzp]->use()&&(xyzp[ixyzp]->use(XYZP::stereopoint)||!_stereoinit))
      //	  accr((pos(ip) - myhel._center).perpPart().mag());
      //      }
      //      myhel._radius = extract_result<tag::median>(accr);
      // restrict the range
      if(_forcep){
	myhel._radius = std::max(std::min(myhel._radius,_rmax),_rmin);
	retval = true;
      } else {
	retval = myhel._radius >= _rmin && myhel._radius <= _rmax;
      }
    }
    return retval;
  }

  void
  HelixFit::fillHitXYZP(TrkDef const& mytrk, std::vector<HitXYZP>& xyzp) {
    // calibration and tracker 
    const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    // loop over straw hits, and store their positions
    for(std::vector<hitIndex>::const_iterator istr=mytrk.strawHitIndices().begin();
	istr != mytrk.strawHitIndices().end(); ++istr){
      StrawHit const& sh = mytrk.strawHitCollection()->at(istr->_index);
      Straw const& straw= tracker.getStraw(sh.strawIndex());
      SHInfo shinfo;
      tcal->StrawHitInfo(straw,sh,shinfo);
      xyzp.push_back(HitXYZP(istr->_index,shinfo._pos,straw.getDirection(),shinfo._tdres,_sfac*straw.getRadius()));
    } 
  }

  void
  HelixFit::findStereoPairs(std::vector<HitXYZP>& xyzp, std::vector<StereoXYZP>& pairs) {
    size_t nxy = xyzp.size();
    for(size_t ixy=0;ixy<nxy;++ixy){
      if(xyzp[ixy].use()){
	for(size_t jxy=ixy+1;jxy<nxy;++jxy){
	  // don't reuse hits already part of a stereo pair, unless requested	  
	  if(xyzp[jxy].use() || _allstereo ){
	    HitXYZP& p1 = xyzp[ixy];
	    HitXYZP& p2 = xyzp[jxy];
	    // ignore nearly parallel wires or hits in different stations or hits in the same plane
	    double dz = fabs(p1._pos.z()-p2._pos.z());  
	    if(fabs(p1._wdir.dot(p2._wdir)) < _maxdot && dz > 0.0 && dz < _maxdz ) {
	      StereoXYZP sp(p1,p2);
	      // find the 2-d intersection between these lints
	      double rho = sp._pos.perp();
	      // if pair within TD limits and the physical detector, record
	      if(fabs(sp._d1) < _nssigma*p1._werr && fabs(sp._d2) < _nssigma*p2._werr &&
		  rho > _rhomin && rho < _rhomax){
		sp.setUse(true);
		pairs.push_back(sp);
		break;
	      }
	    }
	  }
	}
      } 
    }
  }

  void
  HelixFit::findAGE(std::vector<XYZP*> const& xyzp, Hep3Vector const& center,double& rmed, double& age,bool useweights) {
    using namespace boost::accumulators;
    // protection against empty data
    if(xyzp.size() == 0)return;
    // fill radial information for all points, given this center
    std::vector<VALERR> radii;
    unsigned nxyzp = xyzp.size();
    double wtot(0.0);
    for(unsigned ixyzp=0; ixyzp < nxyzp; ++ixyzp){
      if(xyzp[ixyzp]->use()){
	// find radial information for this point
	VALERR rad;
	xyzp[ixyzp]->rinfo(center,rad);
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
  HelixFit::fillSums(std::vector<XYZP*> const& xyzp, Hep3Vector const& center,double rmed,SUMS& sums,bool useweights) {
    // initialize sums
    sums.clear();
    // compute the transverse sums
    double wtot(0.0);
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp]->use()){
	// find radial information for this point
	VALERR rad;
	xyzp[ixyzp]->rinfo(center,rad);
	double wt = useweights ? 1.0/rad._err : 1.0;
	wtot += wt;
	// now x,y projections
	double pcos = (xyzp[ixyzp]->_pos.x()-center.x())/rad._val;
	double psin = (xyzp[ixyzp]->_pos.y()-center.y())/rad._val;
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
  HelixFit::filterDist(std::vector<XYZP*>& xyzp) {
    using namespace boost::accumulators;
    static const double pi(M_PI);
    static const double twopi(2*pi);
    // first, resolve phi.   Use the average X and Y to define the initial
    // phi value, to avoid looping issues
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accx;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accy;
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp]->use()){
	accx(xyzp[ixyzp]->_pos.x());
	accy(xyzp[ixyzp]->_pos.y());
      }
    }
    double mx = extract_result<tag::median>(accx);
    double my = extract_result<tag::median>(accy);
    double mphi = atan2(my,mx);
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp]->use()){
	double dphi = xyzp[ixyzp]->_phi - mphi;
	if(fabs(dphi) > pi){
	  if(dphi > 0)
	    xyzp[ixyzp]->_phi -= twopi;
	  else
	    xyzp[ixyzp]->_phi += twopi;
	}
      }
    }
    // now cut
    CLHEP::Hep3Vector mh(mx,my,0.0);
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      double dist = sqrt(xyzp[ixyzp]->_pos.perpPart().diff2(mh));
      if(dist > _maxdist){
	xyzp[ixyzp]->setUse(false);
	xyzp[ixyzp]->setUse(true,XYZP::outlier);
      }
    }
  }

  void
   HelixFit::filterXY(std::vector<XYZP*>& xyzp, Hep3Vector const& center,double rmed,bool& changed) {
     changed = false;
     for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
       VALERR rad;
       xyzp[ixyzp]->rinfo(center,rad);
       bool use = fabs(rad._val -rmed) < _nsigma*rad._err;
       bool olduse = xyzp[ixyzp]->use();
       xyzp[ixyzp]->setUse(use);
       changed |= olduse != xyzp[ixyzp]->use();
     }
   }

  double
  HelixFit::deltaPhi(double phi1, double phi2){
    static const double pi(M_PI);
    static const double twopi(2*pi);
    double dphi = fmod(phi2-phi1,twopi);
    if(dphi>pi)dphi -= twopi;
    if(dphi<-pi)dphi += twopi;
    return dphi;
  }
}


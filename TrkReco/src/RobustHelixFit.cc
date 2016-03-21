//
// Object to perform helix fit to straw hits
//
// $Id: HelixFit.cc,v 1.12 2014/07/10 14:47:26 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/07/10 14:47:26 $
//
//
// the following has to come before other BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/RobustHelixFit.hh"
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
// C++
#include <vector>
#include <string>
using namespace CLHEP;
namespace mu2e 
{
// statics
  double XYZP::_efac(1.0);
  StrawHitFlag XYZP::_useflag;
  StrawHitFlag XYZP::_dontuseflag;
  int XYZP::_debug(0);
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
    _perr(_efac*shp.posRes(StrawHitPosition::phi)),_rerr(_efac*shp.posRes(StrawHitPosition::rho)),_conversion(false)
  {
    static const CLHEP::Hep3Vector _zdir(0.0,0.0,1.0);
    _sdir = _zdir.cross(_wdir);
  }

  XYZP::XYZP(CLHEP::Hep3Vector const& pos, double size) :  _ind(-1),_pos(pos),
  _phi(0.0),
  _flag(StrawHitFlag::stereo),
  _wdir(CLHEP::Hep3Vector(1,0,0)),_sdir(CLHEP::Hep3Vector(0,1,0)),
  _perr(size),_rerr(size),_conversion(false) {}

  XYZP::XYZP(size_t ind,CLHEP::Hep3Vector const& pos, CLHEP::Hep3Vector const& wdir,
      double werr, double serr) :
    _ind(ind),_pos(pos),_phi(_pos.phi()),_wdir(wdir),_sdir(wdir.y(),-wdir.x(),0.0),
    _perr(_efac*werr),_rerr(_efac*serr),_conversion(false) {}

 void
  XYZP::rinfo(CLHEP::Hep3Vector const& center,VALERR& rad) const {
//    static const double onethird(1.0/3.0);
//    static const double invsqrt12(1./sqrt(12.0));
// average the 1-sigma radii to account for non-linear errors
    double rvec = CLHEP::Hep3Vector(_pos - center).perp();
//    rad._val = onethird*(rvec+rvec1+rvec2);
    rad._val = rvec;
    rad._err = _rerr;
    if(_debug > 1)std::cout << "rinfo : r = " << rad._val << " rerr = " << rad._err  << std::endl;
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
    if(_debug > 1)std::cout << "finfo : phi = " << phi._val << " ferr = " << phi._err << std::endl;
 }

  bool 
  XYZP::use() const {
    return (!_flag.hasAnyProperty(_dontuseflag))
      && (_flag.hasAllProperties(_useflag) || _useflag.empty());
  }

  bool 
  XYZP::stereo() const {
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

  void
  XYZP::setOutlier(){
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    _flag.merge(outlier);
  }

  HelixFitResult& 
    HelixFitResult::operator =(HelixFitResult const& other) {
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
  RobustHelixFit::helixParams(HelixFitResult const& helix,CLHEP::HepVector& pvec,CLHEP::HepVector& perr) const {
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

  RobustHelixFit::RobustHelixFit(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _mindelta(pset.get<double>("minDelta",5000.0)),
    _minnhit(pset.get<unsigned>("minNHit",10)),
    _lambda0(pset.get<double>("lambda0",1.0)),
    _lstep(pset.get<double>("lstep",0.2)),
    _minlambda(pset.get<double>("minlambda",0.01)),
    _maxniter(pset.get<unsigned>("maxniter",50)),
    _nsigma(pset.get<double>("nsigma",5.0)),
    _minzsep(pset.get<double>("minzsep",100.0)),
    _maxzsep(pset.get<double>("maxzsep",700.0)),
    _rbias(pset.get<double>("radialBias",0.0)),
    _efac(pset.get<double>("ErrorFactor",1.0)),
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
    _stereoinit(pset.get<bool>("stereoinit",false)),
    _stereofit(pset.get<bool>("stereofit",false)),
    _targetinit(pset.get<bool>("targetinit",true)),
    _targetinter(pset.get<bool>("targetintersect",true)),
    _targetradius(pset.get<double>("targetradius",75.0)),
    _trackerradius(pset.get<double>("trackerradius",700.0)),
    _bz(0.0)
  {
    XYZP::_efac = _efac;
    std::vector<std::string> bitnames;
    bitnames.push_back("Outlier");
    bitnames.push_back("OtherBackground");
    XYZP::_dontuseflag = StrawHitFlag(bitnames);
    if(_stereofit)XYZP::_useflag = StrawHitFlag(StrawHitFlag::stereo);
    XYZP::_debug = _debug;
}

  RobustHelixFit::~RobustHelixFit()
  {}

  double
  RobustHelixFit::bz() const {
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
  RobustHelixFit::findHelix(HelixFitResult& myhel,bool plothelix) {
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
// call down
    bool retval = findHelix(xyzp,myhel);
    if(_diag>1 && plothelix){
      // fill graphs for display if requested
      plotXY(mytrk,xyzp,myhel);
      plotZ(mytrk,xyzp,myhel);
    }
    return retval;
  }

  bool
  RobustHelixFit::findHelix(XYZPVector& xyzp,HelixFitResult& myhel) {
    bool retval(false);
// filter by geometry
    if(_filter)filterDist(xyzp);
    if(xyzp.size() >= _minnhit){
      // initialize the circle parameters
      if(initCircle(xyzp,myhel)){
	// solve for the circle parameters
	retval = findXY(xyzp,myhel);
	// extend those into the z direction
	if(retval){
	  retval = findZ(xyzp,myhel);
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
  RobustHelixFit::findXY(XYZPVector& xyzp,HelixFitResult& myhel) {
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
  RobustHelixFit::plotXY(HelixDef const& mytrk, XYZPVector const& xyzp,
    HelixFitResult const& myhel) const {

    // Check that there are more than 10 CE hits first
    int n_ce_hits = 0;
    for (XYZPVector::const_iterator i_hit = xyzp.begin(); i_hit != xyzp.end(); ++i_hit) {
      if ((*i_hit).conversion()) {
	++n_ce_hits;
      }
    }

    if (n_ce_hits > 10) {
      static unsigned igraph = 0;
      igraph++;
      art::ServiceHandle<art::TFileService> tfs;
      char ce_stereo_used_name[100];
      snprintf(ce_stereo_used_name,100,"ce_stereo_used_shxy%i",igraph);
      char ce_stereo_notused_name[100];
      snprintf(ce_stereo_notused_name,100,"ce_stereo_notused_shxy%i",igraph);
      char ce_notstereo_used_name[100];
      snprintf(ce_notstereo_used_name,100,"ce_notstereo_used_shxy%i",igraph);
      char ce_notstereo_notused_name[100];
      snprintf(ce_notstereo_notused_name,100,"ce_notstereo_notused_shxy%i",igraph);
      char bkg_stereo_used_name[100];
      snprintf(bkg_stereo_used_name,100,"bkg_stereo_used_shxy%i",igraph);
      char bkg_stereo_notused_name[100];
      snprintf(bkg_stereo_notused_name,100,"bkg_stereo_notused_shxy%i",igraph);
      char bkg_notstereo_used_name[100];
      snprintf(bkg_notstereo_used_name,100,"bkg_notstereo_used_shxy%i",igraph);
      char bkg_notstereo_notused_name[100];
      snprintf(bkg_notstereo_notused_name,100,"bkg_notstereo_notused_shxy%i",igraph);
      char title[100];
      snprintf(title,100,"StrawHit XY trk %i;mm;rad",igraph);
      TH2F* ce_stereo_used = tfs->make<TH2F>(ce_stereo_used_name,title,100,-500,500,100,-500,500);
      TH2F* ce_stereo_notused = tfs->make<TH2F>(ce_stereo_notused_name,title,100,-500,500,100,-500,500);
      TH2F* ce_notstereo_used = tfs->make<TH2F>(ce_notstereo_used_name,title,100,-500,500,100,-500,500);
      TH2F* ce_notstereo_notused = tfs->make<TH2F>(ce_notstereo_notused_name,title,100,-500,500,100,-500,500);
      TH2F* bkg_stereo_used = tfs->make<TH2F>(bkg_stereo_used_name,title,100,-500,500,100,-500,500);
      TH2F* bkg_stereo_notused = tfs->make<TH2F>(bkg_stereo_notused_name,title,100,-500,500,100,-500,500);
      TH2F* bkg_notstereo_used = tfs->make<TH2F>(bkg_notstereo_used_name,title,100,-500,500,100,-500,500);
      TH2F* bkg_notstereo_notused = tfs->make<TH2F>(bkg_notstereo_notused_name,title,100,-500,500,100,-500,500);

      ce_stereo_used->SetMarkerStyle(kFullTriangleUp);
      ce_stereo_used->SetMarkerColor(kRed);
      ce_stereo_notused->SetMarkerStyle(kOpenTriangleUp);
      ce_stereo_notused->SetMarkerColor(kRed);
      ce_notstereo_used->SetMarkerStyle(kFullCircle);
      ce_notstereo_used->SetMarkerColor(kRed);
      ce_notstereo_notused->SetMarkerStyle(kOpenCircle);
      ce_notstereo_notused->SetMarkerColor(kRed);
      bkg_stereo_used->SetMarkerStyle(kFullTriangleUp);
      bkg_stereo_used->SetMarkerColor(kGreen);
      bkg_stereo_notused->SetMarkerStyle(kOpenTriangleUp);
      bkg_stereo_notused->SetMarkerColor(kGreen);
      bkg_notstereo_used->SetMarkerStyle(kFullCircle);
      bkg_notstereo_used->SetMarkerColor(kGreen);
      bkg_notstereo_notused->SetMarkerStyle(kOpenCircle);
      bkg_notstereo_notused->SetMarkerColor(kGreen);

      for(unsigned ih=0;ih<xyzp.size();++ih){
	if(xyzp[ih].conversion()){
	  if (xyzp[ih].use()) {
	    if (xyzp[ih].stereo()) {
	      ce_stereo_used->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }
	    else {
	      ce_notstereo_used->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }	      
	  }
	  else {
	    if (xyzp[ih].stereo()) {
	      ce_stereo_notused->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }
	    else {
	      ce_notstereo_notused->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }
	  }
	}
	else {
	  if (xyzp[ih].use()) {
	    if (xyzp[ih].stereo()) {
	      bkg_stereo_used->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }
	    else {
	      bkg_notstereo_used->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }	      
	  }
	  else {
	    if (xyzp[ih].stereo()) {
	      bkg_stereo_notused->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }
	    else {
	      bkg_notstereo_notused->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }
	  }
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

      TArc* target = new TArc(-myhel._center.x(),-myhel._center.y(),_targetradius);
      target->SetLineColor(kBlack);
      target->SetFillStyle(0);
      // add these to the plot
      TList* flist = ce_stereo_used->GetListOfFunctions();
      flist->Add(fitarc);
      flist->Add(indet);
      flist->Add(outdet);
      flist->Add(target);

      if (mytrk.strawDigiMCCollection() != 0) {
	// Plot the MC true CE hits
	char mctruth_name[100];
	snprintf(mctruth_name,100,"mctshxy%i",igraph);
	TH2F* mct = tfs->make<TH2F>(mctruth_name,title,100,-500,500,100,-500,500);
	mct->SetMarkerStyle(5);
	mct->SetMarkerColor(kMagenta);
	
	for(std::vector<hitIndex>::const_iterator istr=mytrk.strawHitIndices().begin();
	    istr != mytrk.strawHitIndices().end(); ++istr){

	  StrawDigiMC const& mcdigi = mytrk.strawDigiMCCollection()->at(istr->_index);
	  // use TDC channel 0 to define the MC match
	  StrawDigi::TDCChannel itdc = StrawDigi::zero;
	  if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
	  art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	  art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	  int gid(-1);
	  if(spp->genParticle().isNonnull())
	    gid = spp->genParticle()->generatorId().id();
	
	  bool conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);
	  if (conversion) {
	    mct->Fill(spmcp->position().x()-myhel._center.x(),spmcp->position().y()-myhel._center.y());
	  }
	}
      }
    }
  }

  bool
  RobustHelixFit::findCenterAGE(XYZPVector const& xyzp,Hep3Vector& center, double& rmed, double& age,bool useweights) {
// this algorithm follows the method described in J. Math Imagin Vis Dec. 2010 "Robust Fitting of Circle Arcs" (Volume 40, Issue 2, pp. 147-161)
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
// if we're constraining to intersect the target, adjust the center and radius if necessary
      if(_targetinter) {
	double rperigee = center.perp()-rmed;
	if(fabs(rperigee) > _targetradius){
// adjust both center position and radius till they touch the target, holding phi constant.  This keeps the circle near the hits.  Sign matters!
	  double dr;
	  if(rperigee > 0)
	    dr = 0.5*(rperigee - _targetradius);  // change is 1/2 the difference
	  else
	    dr = 0.5*(rperigee + _targetradius);
	  rmed += dr;
	  // direction radially outwards from origin to the center
	  Hep3Vector pdir = Hep3Vector(center.x(),center.y(),0.0).unit();
	  // the center moves opposite the radius
	  center -= dr*pdir;
	}
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
  RobustHelixFit::findZ(XYZPVector& xyzp,HelixFitResult& myhel) {
    using namespace boost::accumulators;
// sort points by z
    std::sort(xyzp.begin(),xyzp.end(),zcomp());
// find phi information
    std::vector<FZ> finfo;
    finfo.reserve(xyzp.size());
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp].use() && (xyzp[ixyzp].stereo() || (!_stereoinit)) ){
	FZ fz;
	xyzp[ixyzp].finfo(myhel._center,fz._phi);
	fz._z = xyzp[ixyzp]._pos.z();
	finfo.push_back(fz);
      }
    }
// find the initial slope
    if(finfo.size() > _minnhit){
      // make initial estimate of dfdz using 'nearby' pairs
      accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > accf;
//      accumulator_set<double, stats<tag::mean > > acctest;
      for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
	for(unsigned jphi=iphi+1; jphi < finfo.size(); ++jphi){
	  double dz = finfo[jphi]._z - finfo[iphi]._z;
	  if(dz > _minzsep && dz < _maxzsep){
	    double dphi = deltaPhi(finfo[iphi]._phi._val,finfo[jphi]._phi._val);
	    double slope = dphi/dz;
	    if(slope > _smin && slope < _smax){ 
	      double wt = _zweights ? dz/(finfo[iphi]._phi._err+finfo[jphi]._phi._err) : 1.0;
	      accf(slope,weight=wt);
//	      acctest(slope);
	    }
	  }
	}
      }
      double  dfdz = extract_result<tag::weighted_median>(accf);
//      double dfdztest = extract_result<tag::mean>(acctest);
      // if the sign of dfdz disagrees, abort
      if( dfdz * _dfdzsign < 0.0)
	return false;
      // if requested, restrict the range
      if(_forcep)
	dfdz = std::max(std::min(dfdz,_smax),_smin);
      else
	if(dfdz > _smax || dfdz < _smin) return false;
// find phi at z intercept.  Use a histogram technique since phi looping
// hasn't been resolved yet
      TH1F hphi("hphi","phi value",50,-1.1*pi,1.1*pi);
      for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
      	double phiex = finfo[iphi]._z*dfdz;
	double dphi = deltaPhi(phiex,finfo[iphi]._phi._val);
	hphi.Fill(dphi);
	hphi.Fill(dphi-twopi);
	hphi.Fill(dphi+twopi);
      }
      double fz0 = hphi.GetBinCenter(hphi.GetMaximumBin());
// refine this using the median.  Not sure this helps, needs to be tested
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
// reset for full search
      if(_stereoinit){
	finfo.clear();
	for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
	  if(xyzp[ixyzp].use()){
	    FZ fz;
	    xyzp[ixyzp].finfo(myhel._center,fz._phi);
	    fz._z = xyzp[ixyzp]._pos.z();
	    finfo.push_back(fz);
	  }
	}
      }
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
//	accumulator_set<double, stats<tag::mean > > acctest2;
	for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
	  for(unsigned jphi=iphi+1; jphi < finfo.size(); ++jphi){
	    double dz = finfo[jphi]._z -finfo[iphi]._z;
	    if(dz > _minzsep){
	      double dphi = finfo[jphi]._phi._val-finfo[iphi]._phi._val;
	      double dphiex = dz*dfdz;
	      double ferr = finfo[iphi]._phi._err+finfo[jphi]._phi._err;
	      if(!_filter || fabs(dphi-dphiex) < _nsigma*ferr){
		double slope = dphi/dz;
		if(slope > _smin && slope < _smax){ 
// limit the weight so as not to count more than 1 loop
		  double wt = _zweights ? std::min(dz,_maxzsep)/ferr : 1.0;
		  accf2(slope,weight=wt);
//		  acctest2(slope);
		}
	      }
	    }
	  }
	}
	dfdz = extract_result<tag::weighted_median>(accf2);
//	dfdztest = extract_result<tag::mean>(acctest2);
	// find phi at z intercept
	accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > acci2;
//	accumulator_set<double, stats<tag::mean>> acci2test;
	for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
	  double phiex = fz0+finfo[iphi]._z*dfdz;
	  if(!_filter || fabs(finfo[iphi]._phi._val-phiex) < _nsigma*finfo[iphi]._phi._err){
	    double wt = _zweights ? 1.0/finfo[iphi]._phi._err : 1.0;
	    acci2(finfo[iphi]._phi._val - finfo[iphi]._z*dfdz, weight=wt);
//	    acci2test(finfo[iphi]._phi._val - finfo[iphi]._z*dfdz);
	  }
	}
//	double fz0test = fmod(extract_result<tag::mean>(acci2test),twopi);
	fz0 = fmod(extract_result<tag::weighted_median>(acci2),twopi);
	if(fz0>pi)fz0 -= twopi;
	if(fz0<-pi)fz0 += twopi;
      }
      myhel._dfdz = dfdz;
      myhel._fz0 = fz0;
      // fix the phi for the hit points
      for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
	double phiex = myhel._fz0 + xyzp[ixyzp]._pos.z()*myhel._dfdz;
	FZ fz;
	xyzp[ixyzp].finfo(myhel._center,fz._phi);
	int nloop = (int)rint((phiex - fz._phi._val)/twopi);
	xyzp[ixyzp]._phi = fz._phi._val + nloop*twopi;
	if(_filter && fabs(xyzp[ixyzp]._phi-phiex)> _nsigma*fz._phi._err) xyzp[ixyzp].setOutlier();
      }
      return true;
    } else
      return false;
  }

  void
  RobustHelixFit::plotZ(HelixDef const& mytrk, XYZPVector const& xyzp, HelixFitResult const& myhel) const {
    // Check that there are more than 10 CE hits first
    int n_ce_hits = 0;
    for (XYZPVector::const_iterator i_hit = xyzp.begin(); i_hit != xyzp.end(); ++i_hit) {
      if ((*i_hit).conversion()) {
	++n_ce_hits;
      }
    }

    if (n_ce_hits > 10) {
      static unsigned igraph = 0;
      igraph++;
      art::ServiceHandle<art::TFileService> tfs;

      char ce_stereo_used_name[100];
      snprintf(ce_stereo_used_name,100,"ce_stereo_used_shphiz%i",igraph);
      char ce_stereo_notused_name[100];
      snprintf(ce_stereo_notused_name,100,"ce_stereo_notused_shphiz%i",igraph);
      char ce_notstereo_used_name[100];
      snprintf(ce_notstereo_used_name,100,"ce_notstereo_used_shphiz%i",igraph);
      char ce_notstereo_notused_name[100];
      snprintf(ce_notstereo_notused_name,100,"ce_notstereo_notused_shphiz%i",igraph);
      char bkg_stereo_used_name[100];
      snprintf(bkg_stereo_used_name,100,"bkg_stereo_used_shphiz%i",igraph);
      char bkg_stereo_notused_name[100];
      snprintf(bkg_stereo_notused_name,100,"bkg_stereo_notused_shphiz%i",igraph);
      char bkg_notstereo_used_name[100];
      snprintf(bkg_notstereo_used_name,100,"bkg_notstereo_used_shphiz%i",igraph);
      char bkg_notstereo_notused_name[100];
      snprintf(bkg_notstereo_notused_name,100,"bkg_notstereo_notused_shphiz%i",igraph);
      char title[100];
      snprintf(title,100,"StrawHit #phi Z trk %i;mm;rad",igraph);
      TH2F* ce_stereo_used = tfs->make<TH2F>(ce_stereo_used_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* ce_stereo_notused = tfs->make<TH2F>(ce_stereo_notused_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* ce_notstereo_used = tfs->make<TH2F>(ce_notstereo_used_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* ce_notstereo_notused = tfs->make<TH2F>(ce_notstereo_notused_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* bkg_stereo_used = tfs->make<TH2F>(bkg_stereo_used_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* bkg_stereo_notused = tfs->make<TH2F>(bkg_stereo_notused_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* bkg_notstereo_used = tfs->make<TH2F>(bkg_notstereo_used_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* bkg_notstereo_notused = tfs->make<TH2F>(bkg_notstereo_notused_name,title,100,-1500,1500,100,-12.5,12.5);

      ce_stereo_used->SetMarkerStyle(kFullTriangleUp);
      ce_stereo_used->SetMarkerColor(kRed);
      ce_stereo_notused->SetMarkerStyle(kOpenTriangleUp);
      ce_stereo_notused->SetMarkerColor(kRed);
      ce_notstereo_used->SetMarkerStyle(kFullCircle);
      ce_notstereo_used->SetMarkerColor(kRed);
      ce_notstereo_notused->SetMarkerStyle(kOpenCircle);
      ce_notstereo_notused->SetMarkerColor(kRed);
      bkg_stereo_used->SetMarkerStyle(kFullTriangleUp);
      bkg_stereo_used->SetMarkerColor(kGreen);
      bkg_stereo_notused->SetMarkerStyle(kOpenTriangleUp);
      bkg_stereo_notused->SetMarkerColor(kGreen);
      bkg_notstereo_used->SetMarkerStyle(kFullCircle);
      bkg_notstereo_used->SetMarkerColor(kGreen);
      bkg_notstereo_notused->SetMarkerStyle(kOpenCircle);
      bkg_notstereo_notused->SetMarkerColor(kGreen);

      for(unsigned ih=0;ih<xyzp.size();++ih){
	if(xyzp[ih].conversion()){
	  if (xyzp[ih].use()) {
	    if (xyzp[ih].stereo()) {
	      ce_stereo_used->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }
	    else {
	      ce_notstereo_used->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }	      
	  }
	  else {
	    if (xyzp[ih].stereo()) {
	      ce_stereo_notused->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }
	    else {
	      ce_notstereo_notused->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }
	  }
	}
	else {
	  if (xyzp[ih].use()) {
	    if (xyzp[ih].stereo()) {
	      bkg_stereo_used->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }
	    else {
	      bkg_notstereo_used->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }	      
	  }
	  else {
	    if (xyzp[ih].stereo()) {
	      bkg_stereo_notused->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }
	    else {
	      bkg_notstereo_notused->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }
	  }
	}
      }
      TF1* line = new TF1("line","[0]+[1]*x",-1500,1500);
      line->SetParameter(0,myhel._fz0);
      line->SetParameter(1,myhel._dfdz);
      line->SetLineColor(kRed);
      TList* flist = ce_stereo_used->GetListOfFunctions();
      flist->Add(line);

      if (mytrk.strawDigiMCCollection() != 0) {
	// Plot the MC true CE hits
	char mctruth_name[100];
	snprintf(mctruth_name,100,"mctshphiz%i",igraph);
	TH2F* mct = tfs->make<TH2F>(mctruth_name,title,100,-1500,1500,100,-12.5,12.5);
	mct->SetMarkerStyle(5);
	mct->SetMarkerColor(kMagenta);
	
	for(std::vector<hitIndex>::const_iterator istr=mytrk.strawHitIndices().begin();
	    istr != mytrk.strawHitIndices().end(); ++istr){

	  StrawDigiMC const& mcdigi = mytrk.strawDigiMCCollection()->at(istr->_index);
	  // use TDC channel 0 to define the MC match
	  StrawDigi::TDCChannel itdc = StrawDigi::zero;
	  if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
	  art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	  art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	  int gid(-1);
	  if(spp->genParticle().isNonnull())
	    gid = spp->genParticle()->generatorId().id();
	
	  bool conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);
	  if (conversion) {
	    mct->Fill(spmcp->position().z(),spmcp->position().phi());
	  }
	}
      }
    }
  }

  bool
  RobustHelixFit::initCircle(XYZPVector const& xyzp,HelixFitResult& myhel) {
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
      if(xyzp[ixyzp].use() && (xyzp[ixyzp].stereo() || (!_stereoinit)) ){
	pos.push_back(xyzp[ixyzp]._pos);
      }
    }
    size_t np = pos.size();
    for(size_t ip=0;ip<np;++ip){
      // pre-compute some values
      double ri2 = pow(pos[ip].x(),2) + pow(pos[ip].y(),2);
      for(size_t jp=ip+1;jp<np; ++jp){
	if(pos[ip].perpPart().diff2(pos[jp].perpPart()) > mind2){
	  double rj2 = pow(pos[jp].x(),2) + pow(pos[jp].y(),2);
	  size_t mink = jp+1;
  // we can force the initialization to use the target position in every triple
	  for(size_t kp=mink;kp<np; ++kp){
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
		double rc = sqrt(cx*cx + cy*cy);
		double rmin = fabs(rc-rho);
		double rmax = rc+rho;
		// test circle parameters for this triple: should be inside the tracker,
		// optionally consistent with the target
		if(rho > _rcmin && rho<_rcmax && rmax < _trackerradius
		  && ( (!_targetinit) || rmin < _targetradius) ) {
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

  void RobustHelixFit::fillXYZP(HelixDef const& mytrk, XYZPVector& xyzp) {
    const Tracker& tracker = getTrackerOrThrow();
    if(mytrk.strawHitPositionCollection() != 0){
    // loop over straw hits, and store their positions
      for(std::vector<hitIndex>::const_iterator istr=mytrk.strawHitIndices().begin();
	  istr != mytrk.strawHitIndices().end(); ++istr){
	StrawHit const& sh = mytrk.strawHitCollection()->at(istr->_index);
	Straw const& straw= tracker.getStraw(sh.strawIndex());
	StrawHitPosition const& shp = mytrk.strawHitPositionCollection()->at(istr->_index);
	XYZP pos(istr->_index,sh,shp,straw);

	// Is this from a conversion hit?
	if(mytrk.strawDigiMCCollection() != 0) {
	  StrawDigiMC const& mcdigi = mytrk.strawDigiMCCollection()->at(istr->_index);
	  // use TDC channel 0 to define the MC match
	  StrawDigi::TDCChannel itdc = StrawDigi::zero;
	  if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
	  art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	  art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	  int gid(-1);
	  if(spp->genParticle().isNonnull())
	    gid = spp->genParticle()->generatorId().id();
	
	  bool conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);
	  if (conversion) {
	    pos.setConversion(true);
	  }
	}

	xyzp.push_back(pos);
      } 
    } else {
      static const double twoinvsqrt12(2.0/sqrt(12.0));
      ConditionsHandle<TrackerCalibrations> tcal("ignored");
      for(std::vector<hitIndex>::const_iterator istr=mytrk.strawHitIndices().begin();
	  istr != mytrk.strawHitIndices().end(); ++istr){
	StrawHit const& sh = mytrk.strawHitCollection()->at(istr->_index);
	Straw const& straw= tracker.getStraw(sh.strawIndex());
	SHInfo shinfo;
	tcal->StrawHitInfo(straw,sh,shinfo);
	xyzp.push_back(XYZP(istr->_index,shinfo._pos,straw.getDirection(),shinfo._tdres,twoinvsqrt12*straw.getRadius()));
      }
    }
  }

  void
  RobustHelixFit::findAGE(XYZPVector const& xyzp, Hep3Vector const& center,double& rmed, double& age,bool useweights) {
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
    // now compute the AGE (Absolute Geometric Error)
    age = 0.0;
    for(unsigned irad=0;irad<radii.size();++irad){
      double wt = useweights ? 1.0/radii[irad]._err : 1.0;
      age += wt*fabs(radii[irad]._val-rmed);
    }
    // normalize
    age *= radii.size()/wtot;
  }

  void
  RobustHelixFit::fillSums(XYZPVector const& xyzp, Hep3Vector const& center,double rmed,SUMS& sums,bool useweights) {
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
  RobustHelixFit::filterDist(XYZPVector& xyzp) {
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
  }

  void
  RobustHelixFit::filterXY(XYZPVector& xyzp, Hep3Vector const& center,double rmed,bool& changed) {
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

  double
  RobustHelixFit::deltaPhi(double phi1, double phi2){
    static const double pi(M_PI);
    static const double twopi(2*pi);
    double dphi = fmod(phi2-phi1,twopi);
    if(dphi>pi)dphi -= twopi;
    if(dphi<-pi)dphi += twopi;
    return dphi;
  }
}


//
// Object to perform helix fit to straw hits
//
// $Id: TrkHelixFit.cc,v 1.9 2012/08/06 16:56:38 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/06 16:56:38 $
//
//
// the following has to come before other BaBar includes
#include "BaBar/BaBar.hh"
#include "TrkPatRec/inc/TrkHelixFit.hh"
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
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
// Root
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TList.h"

namespace mu2e 
{
// comparison functor for ordering points
  struct radcomp : public std::binary_function<RAD, RAD, bool> {
    bool operator()(RAD const& r1, RAD const& r2) { return r1._radius < r2._radius; }
  };

  // comparison functor for sorting by z
  struct zcomp : public std::binary_function<XYZP,XYZP,bool> {
    bool operator()(XYZP const & p1, XYZP const& p2) { return p1._pos.z() < p2._pos.z(); }
  };
  
  void
  XYZP::rinfo(CLHEP::Hep3Vector const& center,RAD& rad) const {
// average the 1-sigma radii to account for non-linear errors
    double rvec = CLHEP::Hep3Vector(_pos - center).perp();
    double rvec1 = CLHEP::Hep3Vector(_pos +_werr*_wdir - center).perp();
    double rvec2 = CLHEP::Hep3Vector(_pos -_werr*_wdir - center).perp();
    rad._radius = 0.5*(rvec1+rvec2);
    rad._rerr = std::max(std::max(fabs(rvec1-rvec),fabs(rvec2-rvec)),_serr);
  }
  
  void
  TrkHelixFit::helixParams(TrkDef const& mytrk, TrkHelix const& helix,CLHEP::HepVector& pvec,CLHEP::HepVector& perr) const {
    static const double pi(M_PI);
    static const double twopi(2*pi);
    static const double halfpi(pi/2.0);
// the helix fit introduces a radial bias due to an asymmetry in the detector (more phase space for
// noise hits outside the circle than inside.  correct for it.
    double radius = helix._radius + _rbias;
    pvec = HepVector(5,0);
// omega is the inverse transverse radius of the particle's circular motion.  Its
// signed by the particle angular momentum about the cirle center.
// This CANNOT be deduced geometrically, so must be supplied as an ad-hoc assumption
    double amsign = copysign(1.0,-mytrk.particle().charge()*_bz);
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
    double dphi = fmod(pvec[HelixTraj::phi0Index]-helix._fz0 - amsign*halfpi,twopi);
// choose z0 (which loop) so that f=0 is as close to z=0 as possible
    if(dphi>pi)dphi -= twopi;
    if(dphi<-pi)dphi += twopi;
    pvec[HelixTraj::z0Index] = dphi*pvec[HelixTraj::tanDipIndex]/pvec[HelixTraj::omegaIndex];
// estimated covariance based on average performance.  These should be parameters, FIXME!!!
    perr = HepVector(5,0);
    perr[HelixTraj::d0Index] = 34.0;
    perr[HelixTraj::phi0Index] = 0.02;
    perr[HelixTraj::omegaIndex]  = 0.0002;
    perr[HelixTraj::tanDipIndex] = 0.05;
    perr[HelixTraj::z0Index] = 15.0;
  }

  TrkHelixFit::TrkHelixFit(fhicl::ParameterSet const& pset) :
  _diag(pset.get<int>("diagLevel",0)),
  _debug(pset.get<int>("debugLevel",0)),
  _mindelta(pset.get<double>("minDelta",5000.0)),
  _minnhit(pset.get<unsigned>("minNHit",10)),
  _maxnhit(pset.get<unsigned>("maxNHit",10)),
  _rfactor(pset.get<double>("rfactor",1.0)),
  _lambda0(pset.get<double>("lambda0",1.0)),
  _lstep(pset.get<double>("lstep",0.2)),
  _minlambda(pset.get<double>("minlambda",0.1)),
  _maxniter(pset.get<unsigned>("maxniter",10)),
  _nsigma(pset.get<double>("nsigma",5)),
  _minzsep(pset.get<double>("minzsep",200)),
  _maxzsep(pset.get<double>("maxzsep",700)),
  _rbias(pset.get<double>("radialBias",-5.0)),
  _sfac(pset.get<double>("strawSizeFactor",2.0)),
  _pmin(pset.get<double>("minP",90)),
  _pmax(pset.get<double>("maxP",115)),
  _tdmin(pset.get<double>("minAbsTanDip",0.5)),
  _tdmax(pset.get<double>("maxAbsTanDip",1.2)),
  _forcep(pset.get<bool>("forceP",true))
    {
    }

  TrkHelixFit::~TrkHelixFit()
    {}

  bool
  TrkHelixFit::findHelix(TrkDef const& mytrk,TrkHelix& myhel) {
    bool retval(false);
// find the magnetic field Z component at the origin
    GeomHandle<BFieldConfig> bfconf;
    _bz = bfconf->getDSUniformValue().z();
//  compute the allowed range in radius for this fit
    double pb = fabs((CLHEP::c_light*1e-3)/(_bz*mytrk.particle().charge()));
    _rmin = _pmin/(pb*sqrt(1.0+_tdmax*_tdmax));
    _rmax = _pmax/(pb*sqrt(1.0+_tdmin*_tdmin));
//  particle charge, field, and direction affect pitch range
    _dfdzsign = copysign(1.0,-mytrk.particle().charge()*mytrk.fitdir().dzdt()*_bz);
// loop over hits, and store the points
    std::vector<XYZP> xyzp;
    fillXYZP(mytrk,xyzp);
// initialize the circle parameters
    if(xyzp.size() > _minnhit){
      if(initCircle(xyzp,myhel)){
// solve for the circle parameters
	retval = findXY(xyzp,myhel);
// fill graphs for display if requested
	if(_diag > 1)plotXY(mytrk,xyzp,myhel);
// extend those into the z direction
	if(retval){
	  retval = findZ(xyzp,myhel);
	  if(_diag > 1)plotZ(mytrk,xyzp,myhel);
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
  TrkHelixFit::findXY(std::vector<XYZP>& xyzp,TrkHelix& myhel) {
    double rmed, age;
    Hep3Vector center = myhel._center;
// first, find the center without weights
    findCenterAGE(xyzp,center,rmed,age,false);
// then, refine that using weights
    bool changed(true);
    unsigned niter(0);
    while(niter < _maxniter && changed){
      findCenterAGE(xyzp,center,rmed,age,true);
      filter(xyzp,center,rmed,changed);
      niter++;
    }
    myhel._center = center;
    myhel._radius = rmed;
    return true;
  }

  void
  TrkHelixFit::plotXY(TrkDef const& mytrk, std::vector<XYZP> const& xyzp,TrkHelix const& myhel) const {
    unsigned igraph = 10*mytrk.eventId()+mytrk.trackId();
    art::ServiceHandle<art::TFileService> tfs;
    //      TGraph* graph = tfs->make<TGraph>(xyzp.size());
    char gname[100];
    snprintf(gname,100,"gshxy%i",igraph);
    char bname[100];
    snprintf(bname,100,"bshxy%i",igraph);
    char title[100];
    snprintf(title,100,"StrawHit XY evt %i trk %i",mytrk.eventId(),mytrk.trackId());
    TH2F* g = tfs->make<TH2F>(gname,title,100,-500,500,100,-500,500);
    TH2F* b = tfs->make<TH2F>(bname,title,100,-500,500,100,-500,500);
    g->SetMarkerStyle(8);
    g->SetMarkerColor(kGreen);
    b->SetMarkerStyle(4);
    b->SetMarkerColor(kBlue);
    //      graph->SetNameTitle(gname,"Straw Hit XY positions");
    for(unsigned ixyzp=0;ixyzp<xyzp.size();++ixyzp){
      //        graph->SetPoint(ixyzp,xyzp[ixyzp]._pos.x()-myhel._center.x(),xyzp[ixyzp]._pos.y()-myhel._center.y());
      if(xyzp[ixyzp]._use)
	g->Fill(xyzp[ixyzp]._pos.x()-myhel._center.x(),xyzp[ixyzp]._pos.y()-myhel._center.y());
      else
	b->Fill(xyzp[ixyzp]._pos.x()-myhel._center.x(),xyzp[ixyzp]._pos.y()-myhel._center.y());
    }
    // need 2 TF1 to model circle as root doesn't support parametric functions
    TF1* circ1 = new TF1("circ1","sqrt([0]*[0]-x*x)",-myhel._radius,myhel._radius);
    circ1->SetParameter(0,myhel._radius);
    circ1->SetLineColor(kRed);
    TF1* circ2 = new TF1("circ2","-sqrt([0]*[0]-x*x)",-myhel._radius,myhel._radius);
    circ2->SetParameter(0,myhel._radius);
    circ2->SetLineColor(kRed);
    TList* flist = g->GetListOfFunctions();
    flist->Add(circ1);
    flist->Add(circ2);
  }

  bool
  TrkHelixFit::findCenterAGE(std::vector<XYZP> const& xyzp,Hep3Vector& center, double& rmed, double& age,bool useweights) {
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
      findAGE(xyzp,cnew,rmed,agenew);
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
// check for errors
        if(agenew > age){
//          std::cout << "error: did not improve AGE!!! Aborting" << std::endl;
//          return false;
        }
      }
// prepare for next iteration
      if(agenew < age){
        center = cnew;
        age = agenew;
      } else
        break;
      ++niter;
    }
// check for convergence
    if(niter > _maxniter){
//      std::cout << "AGE didn't converge!!! " << std::endl;
//      return false;
    }
    return true;
  }
  
  bool
  TrkHelixFit::findZ(std::vector<XYZP>& xyzp,TrkHelix& myhel) {
    using namespace boost::accumulators;
 //
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      xyzp[ixyzp]._phi = atan2((xyzp[ixyzp]._pos.y()-myhel._center.y()),(xyzp[ixyzp]._pos.x()-myhel._center.x()));
    }
// sort these by z
    std::sort(xyzp.begin(),xyzp.end(),zcomp());
// make initial estimate of dfdz using 'nearby' pairs
    std::vector<double> slopes;
    slopes.reserve(4*xyzp.size());
    static const double pi(M_PI);
    static const double twopi(2*pi);
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp]._use){
	for(unsigned jxyzp=ixyzp+1; jxyzp < xyzp.size(); ++jxyzp){
	  if(xyzp[jxyzp]._use){
	    double dz = xyzp[jxyzp]._pos.z() -xyzp[ixyzp]._pos.z();
	    if(fabs(dz) > _minzsep && fabs(dz) < _maxzsep){
	      double dphi = fmod(xyzp[jxyzp]._phi-xyzp[ixyzp]._phi,twopi);
	      if(dphi>pi)dphi -= twopi;
	      if(dphi<-pi)dphi += twopi;
	      slopes.push_back( dphi/dz );
	    }
	  }
        }
      }
    }
    if(slopes.size()>0){
      accumulator_set<double, stats<tag::median(with_p_square_quantile) > > acc;
      acc = std::for_each( slopes.begin(), slopes.end(), acc );
      double dfdz = extract_result<tag::median>(acc);
// if the sign of dfdz disagrees, abort
      if( dfdz * _dfdzsign < 0.0)
	return false;
// if requested, restrict the range
      if(_forcep){
	if(fabs(dfdz) > 1.0/(_rmin*_tdmin))
	  dfdz = copysign(1.0/(_rmin*_tdmin),dfdz);
   	if(fabs(dfdz) < 1.0/(_rmax*_tdmax))
	  dfdz = copysign(1.0/(_rmax*_tdmax),dfdz);
      } else {
	if(fabs(dfdz) > 1.0/(_rmin*_tdmin)) return false;
   	if(fabs(dfdz) < 1.0/(_rmax*_tdmax)) return false;
      }
      
// choose the middle hit Z value to set the convention
      unsigned icomp(0);
      for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
	if(xyzp[ixyzp]._use){
	  icomp = ixyzp;
	  break;
	}
      }
// Use this slope estimate to correct all phi values so that they are are all on the same helix
// iterate over slope and ambiguity resolution
      bool changed(true);
      unsigned niter(0);
      while(changed && niter < _maxniter){
	changed = false;
	++niter;
	for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
	  double dz = xyzp[ixyzp]._pos.z()-xyzp[icomp]._pos.z();
	  double phiex = xyzp[icomp]._phi + dz*dfdz;
	  int nloop = (int)rint((phiex-xyzp[ixyzp]._phi)/twopi);
	  xyzp[ixyzp]._phi += nloop*twopi;
	  changed |= nloop != 0;
	}
	// make a long-range estimate of slope
	slopes.clear();
	for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
	  if(xyzp[ixyzp]._use){
	    for(unsigned jxyzp=ixyzp+1; jxyzp < xyzp.size(); ++jxyzp){
	      if(xyzp[jxyzp]._use){
		double dz = xyzp[jxyzp]._pos.z() -xyzp[ixyzp]._pos.z();
		if(fabs(dz) > _minzsep){
		  slopes.push_back( (xyzp[jxyzp]._phi-xyzp[ixyzp]._phi)/dz );
		}
	      }
	    }
	  }
	}
	accumulator_set<double, stats<tag::median(with_p_square_quantile) > > acc2;
	acc2 = std::for_each( slopes.begin(), slopes.end(), acc2 );
	dfdz = extract_result<tag::median>(acc2);
      }
      myhel._dfdz = dfdz;
// find phi at z intercept
      std::vector<double> inters;
      for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
	if(xyzp[ixyzp]._use)
	  inters.push_back(xyzp[ixyzp]._phi - xyzp[ixyzp]._pos.z()*dfdz);
      }
      accumulator_set<double, stats<tag::median(with_p_square_quantile) > > acc3;
      acc3 = std::for_each( inters.begin(), inters.end(), acc3 );
      double fz0 = fmod(extract_result<tag::median>(acc3),twopi);
      if(fz0>pi)fz0 -= twopi;
      if(fz0<-pi)fz0 += twopi;
      myhel._fz0 = fz0;
// fix the phi for the hit points
      for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
	double phiex = myhel._fz0 + xyzp[ixyzp]._pos.z()*myhel._dfdz;
	int nloop = (int)rint((phiex-xyzp[ixyzp]._phi)/twopi);
	xyzp[ixyzp]._phi += nloop*twopi;
      }
    } else {
      myhel._dfdz = 0.0;
      myhel._fz0 = 0.0;
      return false;
    }
    return true;
  }
  
  void
  TrkHelixFit::plotZ(TrkDef const& mytrk, std::vector<XYZP> const& xyzp,TrkHelix const& myhel) const {
    unsigned igraph = 10*mytrk.eventId()+mytrk.trackId();
    art::ServiceHandle<art::TFileService> tfs;
    char gname[100];
    snprintf(gname,100,"gshphiz%i",igraph);
    char bname[100];
    snprintf(bname,100,"bshphiz%i",igraph);
    char title[100];
    snprintf(title,100,"StrawHit #phi Z evt %i trk %i",mytrk.eventId(),mytrk.trackId());
    TH2F* g = tfs->make<TH2F>(gname,title,50,-1500,1500,50,-5,20);
    TH2F* b = tfs->make<TH2F>(bname,title,50,-1500,1500,50,-5,20);
    g->SetMarkerStyle(8);
    g->SetMarkerColor(kGreen);
    b->SetMarkerStyle(4);
    b->SetMarkerColor(kBlue);
    for(unsigned ixyzp=0;ixyzp<xyzp.size();++ixyzp){
      //        graph->SetPoint(ixyzp,xyzp[ixyzp]._pos.x()-myhel._center.x(),xyzp[ixyzp]._pos.y()-myhel._center.y());
      if(xyzp[ixyzp]._use)
	g->Fill(xyzp[ixyzp]._pos.z(),xyzp[ixyzp]._phi);
      else
	b->Fill(xyzp[ixyzp]._pos.z(),xyzp[ixyzp]._phi);
    }

    TF1* line = new TF1("line","[0]+[1]*x",-1500,1500);
    line->SetParameter(0,myhel._fz0);
    line->SetParameter(1,myhel._dfdz);
    line->SetLineColor(kRed);
    TList* flist = g->GetListOfFunctions();
    flist->Add(line);
  }

  bool
  TrkHelixFit::initCircle(std::vector<XYZP> const& xyzp,TrkHelix& myhel) {
    bool retval(false);
    using namespace boost::accumulators;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accx, accy;
// use a subset of hits
    unsigned istep = std::max((int)ceil(xyzp.size()/_maxnhit),1);
// form all triples, and compute the circle center for unaligned hits.  I can aford to be choosy
    unsigned ntriple(0);
    unsigned nxyzp = xyzp.size();
    for(unsigned ixyzp=0; ixyzp < nxyzp; ixyzp+=istep){
// pre-compute some values
      double ri2 = pow(xyzp[ixyzp]._pos.x(),2) + pow(xyzp[ixyzp]._pos.y(),2);
      for(unsigned jxyzp=ixyzp+1;jxyzp<nxyzp; jxyzp+=istep){
        double rj2 = pow(xyzp[jxyzp]._pos.x(),2) + pow(xyzp[jxyzp]._pos.y(),2);
        for(unsigned kxyzp=jxyzp+1;kxyzp<nxyzp; kxyzp+=istep){
// this effectively measures the slope difference
          double delta = (xyzp[kxyzp]._pos.x() - xyzp[jxyzp]._pos.x())*(xyzp[jxyzp]._pos.y() - xyzp[ixyzp]._pos.y()) - 
            (xyzp[jxyzp]._pos.x() - xyzp[ixyzp]._pos.x())*(xyzp[kxyzp]._pos.y() - xyzp[jxyzp]._pos.y());
          if(fabs(delta) > _mindelta){
            double rk2 = pow(xyzp[kxyzp]._pos.x(),2) + pow(xyzp[kxyzp]._pos.y(),2);
// find circle center for this triple
            double cx = 0.5* (
              (xyzp[kxyzp]._pos.y() - xyzp[jxyzp]._pos.y())*ri2 + 
              (xyzp[ixyzp]._pos.y() - xyzp[kxyzp]._pos.y())*rj2 + 
              (xyzp[jxyzp]._pos.y() - xyzp[ixyzp]._pos.y())*rk2 ) / delta;

            double cy = -0.5* (
              (xyzp[kxyzp]._pos.x() - xyzp[jxyzp]._pos.x())*ri2 + 
              (xyzp[ixyzp]._pos.x() - xyzp[kxyzp]._pos.x())*rj2 + 
              (xyzp[jxyzp]._pos.x() - xyzp[ixyzp]._pos.x())*rk2 ) / delta;
// accumulate 
            ++ntriple;
	    accx(cx);
	    accy(cy);
          }
        }
      }
    }
    if(ntriple > _minnhit){
      double centx = extract_result<tag::median>(accx);
      double centy = extract_result<tag::median>(accy);
      myhel._center = CLHEP::Hep3Vector(centx,centy,0.0);
// use the center to estimate the radius
      accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accr;
      for(unsigned ixyzp=0; ixyzp<nxyzp; ++ixyzp){
        accr((xyzp[ixyzp]._pos - myhel._center).perpPart().mag());
      }
      myhel._radius = extract_result<tag::median>(accr);
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
  TrkHelixFit::fillXYZP(TrkDef const& mytrk, std::vector<XYZP>& xyzp) {
// calibration and tracker information
    const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
// loop over straw hits, and store their positions
    for(std::vector<hitIndex>::const_iterator istr=mytrk.strawHitIndices().begin();
    istr != mytrk.strawHitIndices().end(); ++istr){
      StrawHit const& sh = mytrk.strawHitCollection()->at(istr->_index);
      CLHEP::Hep3Vector wpos;
      double wtime,wtimeres,tdres;
      tcal->StrawHitInfo(sh,wpos,wtime,tdres,wtimeres);
      const Straw& straw = tracker.getStraw(sh.strawIndex());
      xyzp.push_back(XYZP(wpos,straw.getDirection(),tdres,_sfac*straw.getRadius()));
    } 
  }
  
  void
  TrkHelixFit::findAGE(std::vector<XYZP> const& xyzp, Hep3Vector const& center,double& rmed, double& age,bool useweights) {
// protection against empty data
    if(xyzp.size() == 0)return;
// fill radial information for all points, given this center
    std::vector<RAD> radii;
    unsigned nxyzp = xyzp.size();
    double wtot(0.0);
    for(unsigned ixyzp=0; ixyzp < nxyzp; ++ixyzp){
      if(xyzp[ixyzp]._use){
// find radial information for this point
        RAD rad;
        xyzp[ixyzp].rinfo(center,rad);
        radii.push_back(rad);
// compute the normalization too
        wtot += useweights ? 1.0/rad._rerr : 1.0;
      }
    }
// sort these by radius
    std::sort(radii.begin(),radii.end(),radcomp());
// find the median radius.  Use the weights to interpolate
    double wmid = 0.5*wtot;
    double wsum(0.0);
    rmed = -1.0;
    for(unsigned irad=0;irad<radii.size();++irad){
      double wt = useweights ? 1.0/radii[irad]._rerr : 1.0;
      if(wsum + wt > wmid){
        if(irad >0 ){
          rmed = (radii[irad-1]._radius*(wsum+wt-wmid) + radii[irad]._radius*(wmid-wsum))/wt;
        } else {
// degenerate case; all the weight is in the first entry!  Just use that as the median
          rmed =radii[irad]._radius;
        }
        break;
      }
      wsum += wt;
    }
// if requested force radius into range
    if(_forcep)
      rmed = std::max(std::min(rmed,_rmax),_rmin);
// now compute the AGE
    age = 0.0;
    for(unsigned irad=0;irad<radii.size();++irad){
      double wt = useweights ? 1.0/radii[irad]._rerr : 1.0;
      age += wt*fabs(radii[irad]._radius-rmed);
    }
// normalize
    age *= radii.size()/wtot;
  }
  
  void
  TrkHelixFit::fillSums(std::vector<XYZP> const& xyzp, Hep3Vector const& center,double rmed,SUMS& sums,bool useweights) {
// initialize sums
    sums.clear();
// compute the transverse sums
    double wtot(0.0);
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp]._use){
// find radial information for this point
        RAD rad;
        xyzp[ixyzp].rinfo(center,rad);
        double rerr = useweights ? _rfactor*rad._rerr : 1.0;
        double wt = useweights ? 1.0/rad._rerr : 1.0;
	wtot += wt;
// now x,y projections
        double pcos = (xyzp[ixyzp]._pos.x()-center.x())/rad._radius;
        double psin = (xyzp[ixyzp]._pos.y()-center.y())/rad._radius;
// 3 conditions: either the radius is inside the median, outside the median, or 'on' the median.  We define 'on'
// in terms of the error
        if(fabs(rad._radius -rmed) < rerr ){
          sums._scc += wt*fabs(pcos);
          sums._ssc += wt*fabs(psin);
          ++sums._nc;
        } else if (rad._radius > rmed) {
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
  TrkHelixFit::filter(std::vector<XYZP>& xyzp, Hep3Vector const& center,double rmed,bool& changed) {
    changed = false;
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      RAD rad;
      xyzp[ixyzp].rinfo(center,rad);
      bool use = fabs(rad._radius -rmed) < _nsigma*rad._rerr;
      changed |= use != xyzp[ixyzp]._use;
      xyzp[ixyzp]._use = use;
    }
  }
}


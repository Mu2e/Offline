//
// Object to perform helix fit to straw hits
//
// $Id: TrkHelixFit.cc,v 1.1 2011/09/06 23:38:05 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/09/06 23:38:05 $
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
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TList.h"

namespace mu2e 
{
// comparison functor for ordering points
  struct radcomp : public binary_function<RAD, RAD, bool> {
    bool operator()(RAD const& r1, RAD const& r2) { return r1._radius < r2._radius; }
  };

  // comparison functor for sorting by z
  typedef std::pair<double,double> FZP;
  struct zcomp : public binary_function<FZP,FZP,bool> {
    bool operator()(FZP const & p1, FZP const& p2) { return p1.second < p2.second; }
  };
  
  void
  XYZP::rinfo(CLHEP::Hep3Vector const& center,RAD& rad) const {
    Hep3Vector rvec = (_pos - center).perpPart();
    rad._radius = rvec.perp();
// compute the angle between the radius and the wire direction to get the radial error
    Hep3Vector rhat = rvec.unit();
    double wcos = rhat.dot(_wdir);
    double scos = rhat.dot(_sdir);
    rad._rerr = sqrt(pow(wcos*_werr,2) + pow(scos*_serr,2));
// for now, set this to a constant; descent method can't follow the changes in the real errors.
// I need a more powerful non-linear optimizer, FIXME!!!
//    rad._rerr = 2.0;
  }
  
  void
  TrkHelix::helixParams(CLHEP::HepVector& pvec,CLHEP::HepVector& perr) const {
// fit biases should be parameters, FIXME!!!
    const double rbias(-13);
    const double d0bias(20);
    static const double pi(M_PI);
    static const double twopi(2*pi);
    static const double halfpi(0.5*pi);
// should sign omega by the sign of Bz(0), FIXME!!!!!
    double radius = _radius + rbias;
    pvec = HepVector(5,0);
    pvec[HelixTraj::omegaIndex] = 1.0/radius;
    pvec[HelixTraj::d0Index] = _center.perp() - radius + d0bias;
// account for the convention difference
    pvec[HelixTraj::phi0Index] = atan2(-_center.x(),_center.y());
    pvec[HelixTraj::tanDipIndex] = 1.0/(radius*_dfdz);
    double dphi = pvec[HelixTraj::phi0Index] - _fz0 + 3*halfpi;
    int nloop = (int)rint(dphi/twopi);
    dphi -= nloop*twopi;
    pvec[HelixTraj::z0Index] = radius*pvec[HelixTraj::tanDipIndex]*dphi;
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
  _maxniter(pset.get<unsigned>("maxniter",100)),
  _nsigma(pset.get<double>("nsigma",10)),
  _minzsep(pset.get<double>("minzsep",200)),
  _maxzsep(pset.get<double>("maxzsep",700))
    {}

  TrkHelixFit::~TrkHelixFit()
    {}

  bool
  TrkHelixFit::findHelix(TrkDef const& mytrk,TrkHelix& myhel) {
    bool retval(false);
// loop over hits, and store the points
    std::vector<XYZP> xyzp;
    fillXYZP(mytrk,xyzp);
// initialize the circle parameters
    if(xyzp.size() > _minnhit && initCircle(xyzp,myhel)){
// solve for the circle parameters
      retval = findXY(xyzp,myhel);
// extend those into the z direction
      if(retval) retval = findZ(xyzp,myhel);
// set the success
      if(retval)
	myhel._fit = TrkErrCode(TrkErrCode::succeed);
      else
	myhel._fit = TrkErrCode(TrkErrCode::fail);
//      retval = true;
    }
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
// fill graphs for display if requested
    if(_diag > 1){
      static unsigned igraph(0);
      art::ServiceHandle<art::TFileService> tfs;
//      TGraph* graph = tfs->make<TGraph>(xyzp.size());
      char gname[100];
      snprintf(gname,100,"shxy%i",++igraph);
      TH2F* g = tfs->make<TH2F>(gname,"Straw Hit XY positions",100,-500,500,100,-500,500);
      g->SetMarkerStyle(8);
//      graph->SetNameTitle(gname,"Straw Hit XY positions");
      for(unsigned ixyzp=0;ixyzp<xyzp.size();++ixyzp){
//        graph->SetPoint(ixyzp,xyzp[ixyzp]._pos.x()-myhel._center.x(),xyzp[ixyzp]._pos.y()-myhel._center.y());
        g->Fill(xyzp[ixyzp]._pos.x()-myhel._center.x(),xyzp[ixyzp]._pos.y()-myhel._center.y());
      }
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
    return true;
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
  TrkHelixFit::findZ(std::vector<XYZP> const& xyzp,TrkHelix& myhel) {
  //
    std::vector<FZP > fz;
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp]._use){
        double phi = atan2((xyzp[ixyzp]._pos.y()-myhel._center.y()),(xyzp[ixyzp]._pos.x()-myhel._center.x()));
        fz.push_back(make_pair(phi,xyzp[ixyzp]._pos.z()));
      }
    }
// sort these by z
    std::sort(fz.begin(),fz.end(),zcomp());
// make initial estimate of tandip using 'nearby pairs   
    std::vector<double> slopes;
    slopes.reserve(2*fz.size());
    static const double pi(3.1415965);
    static const double twopi(2*pi);
    for(unsigned ifz=0;ifz<fz.size();++ifz){
      for(unsigned jfz=ifz+1;jfz<fz.size();++jfz){
        double dz = fz[jfz].second -fz[ifz].second;
        if(fabs(dz) > _minzsep && fabs(dz) < _maxzsep){
// take care of phi mapping around 0; this assumes points are never more than 1 loop apart
          if(fabs(fz[ifz].first-fz[jfz].first) > pi){
            if(fz[ifz].first>fz[jfz].first)
              fz[jfz].first += twopi;
            else
              fz[jfz].first -= twopi;
          }
          slopes.push_back( (fz[jfz].first-fz[ifz].first)/dz );
        }
      }
    }
    if(slopes.size()>0){
      std::sort(slopes.begin(),slopes.end());
// find the median could interpolate, but it's false accuracy
      unsigned ihalf = (unsigned)rint(slopes.size()/2.0);
      double dfdz = slopes[ihalf];
// check: if the slope is less than that of a conversion, force it to be the average of a conversion
// this only works for conversions and so should be an option, FIXME!!!
      if(dfdz < 0.004)dfdz = 0.0048; // average of all helices
// Use this estimate to correct all points phi so that they are are all on the same helix
      for(unsigned ifz=1;ifz<fz.size();++ifz){
        double phiexp = fz[0].first + (fz[ifz].second-fz[0].second)*dfdz;
        int nloop = (int)rint((phiexp-fz[ifz].first)/twopi);
        fz[ifz].first += nloop*twopi;
      }
// make a long-range estimate of slope
      slopes.clear();
      for(unsigned ifz=0;ifz<fz.size();++ifz){
        for(unsigned jfz=ifz+1;jfz<fz.size();++jfz){
          double dz = fz[jfz].second -fz[ifz].second;
          if(fabs(dz) > _minzsep){
            slopes.push_back( (fz[jfz].first-fz[ifz].first)/dz );
          }
        }
      }
      std::sort(slopes.begin(),slopes.end());
      ihalf = (unsigned)rint(slopes.size()/2.0);
      dfdz = slopes[ihalf];
      myhel._dfdz = dfdz;
// find phi at z intercept
      std::vector<double> inters;
      for(unsigned ifz=0;ifz<fz.size();++ifz){
        inters.push_back(fz[ifz].first - fz[ifz].second*dfdz);
      }
      std::sort(inters.begin(),inters.end());
      unsigned jhalf = (unsigned)rint(inters.size()/2.0);
      myhel._fz0 = inters[jhalf];
    } else {
      myhel._dfdz = 0.0;
      myhel._fz0 = 0.0;
      return false;
    }
// fill graphs for display if requested
    if(_diag > 1){
      static unsigned igraph(0);
      art::ServiceHandle<art::TFileService> tfs;
//      TGraph* graph = tfs->make<TGraph>(fz.size());
      char gname[100];
			snprintf(gname,100,"shphiz%i",++igraph);
      TH2F* graph = tfs->make<TH2F>(gname,"Straw Hit phi vs Z positions",50,-1500,1500,50,-5,20);
      graph->SetMarkerStyle(8);
//        graph->SetNameTitle(gname,"Straw Hit phi vs Z positions");
      for(unsigned ifz=0;ifz<fz.size();++ifz){
//        graph->SetPoint(ifz,fz[ifz].second,fz[ifz].first);
        graph->Fill(fz[ifz].second,fz[ifz].first);
      }
      TF1* line = new TF1("line","[0]+[1]*x",-1500,1500);
      line->SetParameter(0,myhel._fz0);
      line->SetParameter(1,myhel._dfdz);
      line->SetLineColor(kRed);
      TList* flist = graph->GetListOfFunctions();
      flist->Add(line);
    }
    return true;
  }

  bool
  TrkHelixFit::initCircle(std::vector<XYZP> const& xyzp,TrkHelix& myhel) {
// use a subset of hits
    unsigned istep = max((int)ceil(xyzp.size()/_maxnhit),1);
// form all triples, and compute the circle center for unaligned hits.  I can aford to be choosy
    unsigned ntriple(0);
    double cxsum(0.0), cysum(0.0);
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
// increment sums              
// might need to do an iterative outlier search here someday?
            ++ntriple;
            cxsum += cx;
            cysum += cy;
          }
        }
      }
    }
    if(ntriple > _minnhit){
      myhel._center = CLHEP::Hep3Vector(cxsum/ntriple,cysum/ntriple,0.0);
// use the center to estimate the radius
      double rsum(0.0);
      for(unsigned ixyzp=0; ixyzp<nxyzp; ++ixyzp){
        rsum += (xyzp[ixyzp]._pos - myhel._center).perpPart().mag();
      }
// set a maximum size; this only works for conversions and should be an option, FIXME!!!
      double radius = rsum/xyzp.size();
      if(radius > 300.0) radius = 300.0;
      myhel._radius = radius;
      return true;
    } else
      return false;
  }

  void
  TrkHelixFit::fillXYZP(TrkDef const& mytrk, std::vector<XYZP>& xyzp) {
// convenience factor
    static const double sfac(2.0/sqrt(12.0));
// calibration and tracker information
    const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
// loop over straw hits, and store their XY projection
    for(std::vector<size_t>::const_iterator istr=mytrk.strawHitIndices().begin();
    istr != mytrk.strawHitIndices().end(); ++istr){
      StrawHit const& sh = mytrk.strawHitCollection()->at(*istr);
      CLHEP::Hep3Vector wpos;
      double wtime,wtimeres,tdres;
      tcal->StrawHitInfo(sh,wpos,wtime,tdres,wtimeres);
      const Straw& straw = tracker.getStraw(sh.strawIndex());
      xyzp.push_back(XYZP(wpos,straw.getDirection(),tdres,sfac*straw.getRadius()));
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
// check for other degenerate case: all the weight in the last entry
    if(rmed < 0.0)
      rmed = 270.0;
// force radius into physical range; this is only for conversion electrons, FIXME!!!
    if(rmed > 300.0) rmed = 300.0;
// now compute the AGE
    age = 0.0;
    for(unsigned irad=0;irad<radii.size();++irad){
      double wt = useweights ? 1.0/radii[irad]._rerr : 1.0;
      age += wt*fabs(radii[irad]._radius-rmed);
    }    
  }
  
  void
  TrkHelixFit::fillSums(std::vector<XYZP> const& xyzp, Hep3Vector const& center,double rmed,SUMS& sums,bool useweights) {
// initialize sums
    sums.clear();
// compute the transverse sums
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp]._use){
// find radial information for this point
        RAD rad;
        xyzp[ixyzp].rinfo(center,rad);
        double rerr = useweights ? _rfactor*rad._rerr : 1.0;
        double wt = useweights ? 1.0/rad._rerr : 1.0;
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


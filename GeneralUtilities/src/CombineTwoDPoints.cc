#include "Offline/GeneralUtilities/inc/CombineTwoDPoints.hh"
#include "TMath.h"
#include <stdexcept>
namespace mu2e {
  using SVEC = ROOT::Math::SVector<double,2>;
  using SMAT = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>;

  CombineTwoDPoints::CombineTwoDPoints(std::vector<TwoDPoint> const& points, float ivar) : ivar_(ivar) {
    for(size_t key=0; key < points.size(); ++key) {
      auto const& point = points[key];
      addPoint(point,key);
    }
  }

  TwoDPoint const& CombineTwoDPoints::point() const {
    if(!ptcalc_){
      ptcalc_ = true;
      point_ = wt_.point();
    }
    return point_;
  }

  void CombineTwoDPoints::addPoint(TwoDPoint const& point,size_t key) {
    auto wt = TwoDWeight(point,ivar_);
    double dchi0 = ROOT::Math::Similarity(point.pos(),wt.weight());
    wts_.emplace(std::make_pair(key,CWT(wt,dchi0)));
    wt_ += wt;
    chi0_ += dchi0;
    chicalc_ = false;
    ptcalc_ = false;
  }

  void CombineTwoDPoints::removePoint(size_t key) {
    auto ifnd = wts_.find(key);
    if(ifnd != wts_.end()){
      CWT const& cwt = ifnd->second;
      wt_ -= cwt.wt_;
      chi0_ -= cwt.dchi0_;
      chicalc_ = false;
      ptcalc_ = false;
      wts_.erase(ifnd);
    }
  }

  double CombineTwoDPoints::dChi2(size_t key) const {
    auto const& cwt = wts_.at(key);
    double chi0 = chi0_ - cwt.dchi0_;
    auto const& ptwt = cwt.wt_;
    auto wt = wt_; wt -= ptwt;
    int ifail(0);
    auto cov = wt.wt().InverseFast(ifail);
    if(ifail != 0)throw std::invalid_argument( "Inversion Failure" );
    double chisq = chi0 - ROOT::Math::Similarity(wt.wtPos(),cov);
    return chisquared() - chisq;
  }

  double CombineTwoDPoints::dChi2(TwoDPoint const& pt) const {
    TwoDWeight pwt(pt,ivar_);
    double chi0 = chi0_ + ROOT::Math::Similarity(pt.pos(),pwt.weight());
    auto wt = wt_; wt += pwt;
    int ifail(0);
    auto cov = wt.wt().InverseFast(ifail);
    if(ifail != 0)throw std::invalid_argument( "Inversion Failure" );
    return chi0 - ROOT::Math::Similarity(wt.wtPos(),cov);
  }

  double CombineTwoDPoints::chisquared() const {
    if(!chicalc_){
      chicalc_ = true;
      chisq_ = chi0_ - ROOT::Math::Similarity(wt_.wtPos(),point().cov());
    }
    return chisq_;
  }

  double CombineTwoDPoints::consistency() const{
    auto ndof = nDOF();
    if( ndof > 0 && chisquared() > 0.0)
      return TMath::Prob(chisquared(),ndof);
    else
      return -1.0;
  }

  void CombineTwoDPoints::print(std::ostream& os) const {
    os << "CombineTwoDPoints with " << nPoints() << " points, " << point() << " chisquared " << chisquared() << std::endl;
  }
}

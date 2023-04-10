#include "Offline/GeneralUtilities/inc/CombineTwoDPoints.hh"
#include "TMath.h"
#include <stdexcept>
namespace mu2e {
  CombinedTwoDPoints::CombinedTwoDPoints(TwoDPoint const& point, float ivar) : ivar_(ivar), ptcalc_(true), chicalc_(true), point_(point), wt_(point),  chisq_(0.0) {
    chi0_ = ROOT::Math::Similarity(point.pos(),wt_.weight());
    wts_.emplace(std::make_pair(0,CWT(wt_,chi0_)));
  }

  CombinedTwoDPoints::CombinedTwoDPoints(std::vector<TwoDPoint> points, float ivar) : ivar_(ivar), ptcalc_(false), chicalc_(false), chi0_(0.0), chisq_(0.0){
//    wts_.reserve(points.size());
    for(size_t key=0; key < points.size(); ++key) {
      auto const& point = points[key];
      auto wt = TwoDWeight(point);
      double dchi0 = ROOT::Math::Similarity(point.pos(),wt.weight());
      wts_.emplace(std::make_pair(key,CWT(wt,dchi0)));
      wt_ += wt;
      chi0_ += dchi0;
    }
  }

  TwoDPoint const& CombinedTwoDPoints::point() const {
    if(!ptcalc_){
      ptcalc_ = true;
      point_ = wt_.point();
    }
    return point_;
  }

  size_t CombinedTwoDPoints::addPoint(TwoDPoint const& point) {
    auto wt = TwoDWeight(point);
    double dchi0 = ROOT::Math::Similarity(point.pos(),wt.weight());
    auto key = wts_.size();
    wts_.emplace(std::make_pair(key,CWT(wt,dchi0)));
    wt_ += wt;
    chi0_ += dchi0;
    chicalc_ = false;
    ptcalc_ = false;
    return key;
  }

  void CombinedTwoDPoints::removePoint(size_t key) {
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

  double CombinedTwoDPoints::dChi2(size_t key) const {
    auto const& cwt = wts_.at(key);
    double chi0 = chi0_ - cwt.dchi0_;
    auto const& ptwt = cwt.wt_;
    auto wt = wt_; wt -= ptwt;
    int ifail(0);
    auto cov = wt.wt().Inverse(ifail);
    if(ifail != 0)throw std::invalid_argument( "Inversion Failure" );
    double chisq = chi0 - ROOT::Math::Similarity(wt.wtPos(),cov);
    return chisquared() - chisq;
  }

  double CombinedTwoDPoints::dChi2(TwoDPoint const& pt) const {
    auto dpt = point().pos() - pt.pos();
    auto cov = point().cov() + pt.cov();
    TwoDPoint dpoint(dpt,cov);
    TwoDWeight wt(dpoint,ivar_);
    return ROOT::Math::Similarity(dpoint.pos(),wt.wt());
  }

  double CombinedTwoDPoints::chisquared() const {
    if(!chicalc_){
      chicalc_ = true;
      chisq_ = chi0_ - ROOT::Math::Similarity(wt_.wtPos(),point().cov());
    }
    return chisq_;
  }

  double CombinedTwoDPoints::consistency() const{
    auto ndof = nDOF();
    if( ndof > 0 && chisquared() > 0.0)
      return TMath::Prob(chisquared(),ndof);
    else
      return -1.0;
  }

  void CombinedTwoDPoints::print(std::ostream& os) const {
    os << "CombinedTwoDPoints with " << nPoints() << " points, " << point() << " chisquared " << chisquared() << std::endl;
  }
}

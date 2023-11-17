#include "Offline/TrkHitReco/inc/CombineStereoPoints.hh"
#include "TMath.h"
#include <stdexcept>
namespace mu2e {
  using SVEC = ROOT::Math::SVector<double,2>;
  using SMAT = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>;

  StereoPoint CombineStereoPoints::point(size_t key) const {
    auto const& wt = wts_.at(key);
    return StereoPoint(wt.wt_.point(ivar_),wt.z_);
  }

  StereoPoint const& CombineStereoPoints::point() const {
    if(!ptcalc_){
      ptcalc_ = true;
      // compute the Z position as the trace-weighted average
      double z(0.0), trsum(0.0);
      for (auto const& wt : wts_){
        auto const& wm = wt.second.wt_.weight();
        double tr = wm.Trace();
        z += wt.second.z_*tr;
        trsum += tr;
      }
      z /= trsum;
      auto twodpoint = wt_.point(); // transverse coordinates; add linear fit in 2 directions TODO
      point_ = StereoPoint(twodpoint,z);
    }
    return point_;
  }

  void CombineStereoPoints::addPoint(StereoPoint const& point,size_t key) {
    auto wt = TwoDWeight(point.point(),ivar_);
    double dchi0 = ROOT::Math::Similarity(point.point().pos(),wt.weight());
    wts_.emplace(std::make_pair(key,CWT(wt,point.z(),dchi0)));
    wt_ += wt;
    chi0_ += dchi0;
    chicalc_ = false;
    ptcalc_ = false;
  }

  void CombineStereoPoints::removePoint(size_t key) {
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

  double CombineStereoPoints::dChi2(size_t key) const {
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

  double CombineStereoPoints::dChi2(StereoPoint const& pt) const {
    TwoDWeight pwt(pt.point(),ivar_);
    double chi0 = chi0_ + ROOT::Math::Similarity(pt.point().pos(),pwt.weight());
    auto wt = wt_; wt += pwt;
    int ifail(0);
    auto cov = wt.wt().InverseFast(ifail);
    if(ifail != 0)throw std::invalid_argument( "Inversion Failure" );
    return chi0 - ROOT::Math::Similarity(wt.wtPos(),cov);
  }

  double CombineStereoPoints::chisquared() const {
    if(!chicalc_){
      chicalc_ = true;
      chisq_ = chi0_ - ROOT::Math::Similarity(wt_.wtPos(),point().point().cov());
    }
    return chisq_;
  }

  double CombineStereoPoints::consistency() const{
    auto ndof = nDOF();
    if( ndof > 0 && chisquared() > 0.0)
      return TMath::Prob(chisquared(),ndof);
    else
      return -1.0;
  }

  void CombineStereoPoints::print(std::ostream& os) const {
    os << "CombineStereoPoints with " << nPoints() << " points, " << point() << " chisquared " << chisquared() << std::endl;
  }
}

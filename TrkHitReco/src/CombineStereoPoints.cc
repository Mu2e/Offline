#include "Offline/TrkHitReco/inc/CombineStereoPoints.hh"
#include "TMath.h"
#include <stdexcept>
namespace mu2e {
  StereoPoint CombineStereoPoints::point(size_t key) const {
    auto const& wt = wts_.at(key);
    return StereoPoint(wt.pt_);
  }

  StereoPoint const& CombineStereoPoints::point() const {
    if(!ptcalc_){
      ptcalc_ = true;
      // compute the Z position as the trace-weighted average
      double z(0.0), trsum(0.0);
      for (auto const& wt : wts_){
        auto const& wm = wt.second.wt_.weight();
        double tr = wm.Trace();
        z += wt.second.pt_.z()*tr;
        trsum += tr;
      }
      z /= trsum;
      auto twodpoint = wt_.point(); // transverse coordinates; add linear fit in 2 directions TODO
      point_ = StereoPoint(twodpoint,z);
    }
    return point_;
  }

  void CombineStereoPoints::addPoint(StereoPoint const& point,size_t key) {
    CWT cwt(point,ivar_);
    wts_.emplace(std::make_pair(key,cwt));
    wt_ += cwt.wt_;
    chi0_ += cwt.dchi0_;
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

  bool CombineStereoPoints::stereoLine(StereoLine& sline) const {
    // 2-D weighted linear regression
    bool retval(false);
    // first solve for z0
    double z0(0.0), trsum(0.0);
    for (auto const& wt : wts_){
      auto const& wm = wt.second.wt_.weight();
      double tr = wm.Trace();
      trsum += tr;
      z0 += wt.second.pt_.z()*tr;
    }
    sline.z0_ = z0/trsum;
    // now accumulate the z-dependent sums
    double alpha(0.0);
    StereoLine::SVEC beta;
    StereoLine::SMAT gamma;
    const static StereoLine::SMAT ident =  ROOT::Math::SMatrixIdentity();
    for (auto const& wt : wts_){
      auto const& spt = wt.second.pt_;
      auto const& twodpt = spt.point();
      auto pos = twodpt.pos2();
      auto udir = twodpt.udir();
      auto vdir = twodpt.vdir();
      auto uwt = 1.0/twodpt.uvar();
      auto vwt = 1.0/twodpt.vvar();
      double pu = pos.Dot(udir);
      double pv = pos.Dot(vdir);
      alpha += uwt*pu*pu + vwt*pv*pv;
      double dz = spt.z()-sline.z0_;
      StereoLine::SVEC etau(udir.X(),udir.Y(),dz*udir.X(),dz*udir.Y());
      StereoLine::SVEC etav(vdir.X(),vdir.Y(),dz*vdir.X(),dz*vdir.Y());
      beta += uwt*pu*etau + vwt*pv*etav;
      // cannot use TensorProd as that results in a non-symmetric matrix
      // The following is clumsy but works
      StereoLine::SMAT etau_mat;
      etau_mat.SetDiagonal(etau);
      StereoLine::SMAT etav_mat;
      etav_mat.SetDiagonal(etav);
      gamma += uwt*ROOT::Math::Similarity(etau_mat,ident) + vwt*ROOT::Math::Similarity(etav_mat,ident);
    }
    int inv;
    sline.cov_ = gamma.Inverse(inv);
    if(inv == 0){
      retval = true;
      sline.pars_ = sline.cov_*beta;
      sline.chisq_ = std::max(alpha - ROOT::Math::Similarity(beta,sline.cov_), 0.0);
      sline.ndof_ = wts_.size()*2 - 4;
    }
    return retval;
  }

  void CombineStereoPoints::print(std::ostream& os) const {
    os << "CombineStereoPoints with " << nPoints() << " points, " << point() << " chisquared " << chisquared() << std::endl;
  }
}

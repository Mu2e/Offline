#include "Offline/GeneralUtilities/inc/CombineTwoDPoints.hh"
#include "TMath.h"
#include <stdexcept>
namespace mu2e {
  using SVEC = ROOT::Math::SVector<float,2>;
  using SMAT = ROOT::Math::SMatrix<float,2,2,ROOT::Math::MatRepSym<float,2>>;
  TwoDWeight::TwoDWeight(TwoDPoint const& pos,float ivar) {
    wt_ = pos.weight(ivar);
    wtpos_ = wt_*pos.pos();
  }

  TwoDPoint TwoDWeight::point() const {
    int ifail(0);
    SMAT cov = wt_.Inverse(ifail);
    if(ifail != 0)throw std::invalid_argument( "Inversion Failure" );
    auto pos = cov*wtpos_;
    return TwoDPoint(pos,cov);
  }
  TwoDWeight& TwoDWeight::operator +=(TwoDWeight const& other) {
    wt_ += other.wt();
    wtpos_ += other.wtPos();
    return *this;
  }

  CombinedTwoDPoints::CombinedTwoDPoints(TwoDPoint const& point, float ivar) : ivar_(ivar), point_(point), wt_(point), chisq_(0.0) {
    wts_.push_back(wt_);
    chi0_ = ROOT::Math::Similarity(point.pos(),wts_.back().wt());
  }

  CombinedTwoDPoints::CombinedTwoDPoints(std::vector<TwoDPoint> points, float ivar) : ivar_(ivar) {
    wts_.reserve(points.size());
    for (auto const& point : points ){
      wts_.push_back(TwoDWeight(point));
      chi0_ += ROOT::Math::Similarity(point.pos(),wts_.back().wt());
    }
    for(auto const& wt : wts_){
      wt_ += wt;
    }
    point_ = wt_.point();
    chisq_ = chi0_ - ROOT::Math::Similarity(wt_.wtPos(),point_.cov());
  }

  void CombinedTwoDPoints::addPoint(TwoDPoint const& point) {
    wts_.push_back(TwoDWeight(point,ivar_));
    chi0_ += ROOT::Math::Similarity(point.pos(),wts_.back().wt());
    wt_ += wts_.back();
    point_ = wt_.point();
    chisq_ += - ROOT::Math::Similarity(wt_.wtPos(),point_.cov());
  }

  float CombinedTwoDPoints::dChi2(TwoDPoint const& point) const {
    auto dpt = point_.pos() - point.pos();
    SMAT cov = point_.cov() + point.cov();
    TwoDPoint dpoint(dpt,cov);
    TwoDWeight wt(dpoint,ivar_);
    return ROOT::Math::Similarity(dpoint.pos(),wt.wt());
  }

  double CombinedTwoDPoints::probability() const{
    auto ndof = nDOF();
    if( ndof > 0 && chisq_ > 0.0)
      return TMath::Prob(chisq_,ndof);
    else
      return -1.0;
  }
}

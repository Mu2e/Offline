#include "Offline/GeneralUtilities/inc/CombineTwoDPoints.hh"
#include "TMath.h"
#include <stdexcept>
namespace mu2e {
  using SVEC = ROOT::Math::SVector<double,2>;
  using SMAT = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>;
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

  void TwoDWeight::print(std::ostream& os) const {
    os << *this << std::endl;
  }

  CombinedTwoDPoints::CombinedTwoDPoints(TwoDPoint const& point, float ivar) : ivar_(ivar), ptcalc_(true), chicalc_(true), point_(point), wt_(point),  chisq_(0.0) {
    wts_.push_back(wt_);
    chi0_ = ROOT::Math::Similarity(point.pos(),wts_.back().weight());
  }

  CombinedTwoDPoints::CombinedTwoDPoints(std::vector<TwoDPoint> points, float ivar) : ivar_(ivar), ptcalc_(false), chicalc_(false), chi0_(0.0), chisq_(0.0) {
    wts_.reserve(points.size());
    for (auto const& point : points ){
      wts_.push_back(TwoDWeight(point));
      chi0_ += ROOT::Math::Similarity(point.pos(),wts_.back().weight());
      wt_ += wts_.back();
    }
  }

  TwoDPoint const& CombinedTwoDPoints::point() const {
    if(!ptcalc_){
      ptcalc_ = true;
      point_ = wt_.point();
    }
    return point_;
  }

  void CombinedTwoDPoints::addPoint(TwoDPoint const& point) {
    wts_.push_back(TwoDWeight(point,ivar_));
    wt_ += wts_.back();
    chi0_ += ROOT::Math::Similarity(point.pos(),wts_.back().weight());
    chicalc_ = false;
    ptcalc_ = false;
  }

  double CombinedTwoDPoints::dChi2(TwoDPoint const& pt) const {
    auto dpt = point().pos() - pt.pos();
    SMAT cov = point().cov() + pt.cov();
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

std::ostream& operator << (std::ostream& ost, mu2e::TwoDWeight const& wt) {
  ost << "TwoDWeight position (" << wt.weightPosition()[0] << "," << wt.weightPosition()[1] << "), upper weight matrix ("
    << wt.weight()(0,0) << "," <<  wt.weight()(0,1) << "," << wt.weight()(1,1) << ") ";
  return ost;
}

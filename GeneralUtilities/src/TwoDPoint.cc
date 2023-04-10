#include "Offline/GeneralUtilities/inc/TwoDPoint.hh"
#include <stdexcept>
#include <limits>
namespace mu2e {
  using SMAT = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>;
  TwoDPoint::TwoDPoint(VEC2 const& pos,float ucos, float uvar, float vvar) : pos_(pos.X(),pos.Y()) {
    if(fabs(ucos)>1.0) throw std::invalid_argument( "Unphysical covariance orientation" );
    if( uvar < std::numeric_limits<float>::min() || vvar < std::numeric_limits<float>::min())
      throw std::invalid_argument( "Unphysical variance" );
    double ucos2 = ucos*ucos;
    double usin2 = 1.0-ucos2;
    double usin = sqrt(usin2);
    double cov[3] = {ucos2*uvar + usin2*vvar, ucos*usin*(uvar -vvar),usin2*uvar + ucos2*vvar};
    cov_ = SMAT(cov,cov+3);
  }

  TwoDPoint::TwoDPoint(VEC3 const& pos,VEC3 const& udir, float uvar, float vvar) : TwoDPoint(VEC2(pos.X(),pos.Y()),udir.Dot(VEC3(1.0,0.0,0.0)),uvar,vvar){}

  void TwoDPoint::print(std::ostream& os) const {
    os << *this << std::endl;
  }
}

std::ostream& operator << (std::ostream& ost, mu2e::TwoDPoint const& pt) {
  ost << "TwoDPoint position (" << pt.pos()[0] << "," << pt.pos()[1] << "), upper covariance ("
    << pt.cov()(0,0) << "," << pt.cov()(0,1) << "," << pt.cov()(1,1) << ") ";
  return ost;
}

#include "Offline/GeneralUtilities/inc/TwoDPoint.hh"
#include <stdexcept>
#include <limits>
namespace mu2e {
  using SMAT = ROOT::Math::SMatrix<float,2,2,ROOT::Math::MatRepSym<float,2>>;
  TwoDPoint::TwoDPoint(VEC2 const& pos,float ucos, float uvar, float vvar) : pos_(pos.X(),pos.Y()) {
    if(fabs(ucos)>1.0) throw std::invalid_argument( "Unphysical covariance orientation" );
    if( uvar < std::numeric_limits<float>::min() || vvar < std::numeric_limits<float>::min())
      throw std::invalid_argument( "Unphysical variance" );
    float ucos2 = ucos*ucos;
    float usin2 = 1.0-ucos2;
    float usin = sqrt(usin2);
    float cov[3] = {ucos2*uvar + usin2*vvar, ucos*usin*(uvar -vvar),usin2*uvar + ucos2*vvar};
    cov_ = SMAT(cov,cov+3);
  }

  TwoDPoint::TwoDPoint(VEC3 const& pos,VEC3 const& udir, float uvar, float vvar) : TwoDPoint(VEC2(pos.X(),pos.Y()),udir.Dot(VEC3(1.0,0.0,0.0)),uvar,vvar){}

  SMAT TwoDPoint::weight(float ivar) const {
    SMAT cov(cov_);
    if(ivar > 0.0){
      cov(0,0) += ivar;
      cov(1,1) += ivar;
    }
    int ifail(0);
    SMAT retval = cov.Inverse(ifail);
    if(ifail != 0)throw std::invalid_argument( "Inversion Failure" );
    return retval;
  }
}

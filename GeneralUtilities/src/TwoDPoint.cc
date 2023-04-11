#include "Offline/GeneralUtilities/inc/TwoDPoint.hh"
#include <stdexcept>
#include <limits>
namespace mu2e {
  using SMAT = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>;
  using VEC2 = ROOT::Math::XYVectorF;

  TwoDPoint::TwoDPoint(VEC2 const& pos,VEC2 const& udir, float uvar, float vvar) : pos_(pos.X(),pos.Y()) {
    if( uvar < std::numeric_limits<float>::min() || vvar < std::numeric_limits<float>::min())
      throw std::invalid_argument( "Unphysical variance" );
    auto ud = udir.Unit(); // unit vector
    static const VEC2 xdir(1.0,0.0);
    static const VEC2 ydir(0.0,1.0);
    double ucos = ud.Dot(xdir);
    double ucos2 = ucos*ucos;
    double usin = ud.Dot(ydir);
    double usin2 = usin*usin;
    double cov[3] = {ucos2*uvar + usin2*vvar, ucos*usin*(uvar -vvar),usin2*uvar + ucos2*vvar};
    cov_ = SMAT(cov,cov+3);
  }

  TwoDPoint::TwoDPoint(VEC2 const& pos,UVRes const& uvres) : TwoDPoint(pos, uvres.udir_, uvres.uvar(),uvres.vvar()) {}
  TwoDPoint::TwoDPoint(VEC3 const& pos,VEC3 const& udir, float uvar, float vvar) :
    TwoDPoint(VEC2(pos.X(),pos.Y()), VEC2(udir.X(),udir.Y()),uvar,vvar){}

  UVRes TwoDPoint::uvRes() const {
    // compute eignvalues
    double ldiff = cov_(0,0)-cov_(1,1);
    double lr = sqrt(ldiff*ldiff + 4*cov_(0,1)*cov_(0,1));
    double lsum = cov_(0,0)+cov_(1,1);
    if(lr > lsum) throw std::runtime_error( "Unphysical covariance" );
    UVRes retval;
    retval.ures_ = sqrt(0.5*(lsum+lr));// larger error along the u direction by convention
    retval.vres_ = sqrt(0.5*(lsum-lr));
    // compute angle: this can be degenerate
    if(lr/lsum > 1e-6 ){
      double ux2 = std::min(1.0,0.5*(1.0 + ldiff/lr));
      float ux = sqrt(std::max(0.0,ux2));
      float uy = copysign(sqrt(1.0 -ux2),cov_(0,1));
      retval.udir_ = VEC2(ux,uy);
    } else
      retval.udir_ = VEC2(1.0,0.0); // circular errors, set an arbitrary angle
    return retval;
  }

  void TwoDPoint::print(std::ostream& os) const {
    os << *this << std::endl;
  }
}

std::ostream& operator << (std::ostream& ost, mu2e::TwoDPoint const& pt) {
  ost << "TwoDPoint position (" << pt.pos()[0] << "," << pt.pos()[1] << "), upper covariance ("
    << pt.cov()(0,0) << "," << pt.cov()(0,1) << "," << pt.cov()(1,1) << ") ";
  return ost;
}

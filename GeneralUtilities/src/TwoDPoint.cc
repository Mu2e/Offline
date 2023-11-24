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

  TwoDPoint::TwoDPoint(VEC2 const& pos,UVVar const& uvvar) : TwoDPoint(pos, uvvar.udir_, uvvar.uvar(),uvvar.vvar()) {}
  TwoDPoint::TwoDPoint(VEC3 const& pos,VEC3 const& udir, float uvar, float vvar) :
    TwoDPoint(VEC2(pos.X(),pos.Y()), VEC2(udir.X(),udir.Y()),uvar,vvar){}

  UVVar TwoDPoint::uvRes() const {
    // diagonalize the covariance matrix.  First, compute the eignvalues
    UVVar retval;
    double det;
    double halftr = 0.5*cov_.Trace();
    // protect against degenerate matrices
    bool ok = cov_.Det2(det);
    if(ok){
      double delta =  halftr*halftr - det;
      ok = delta > 0.0;
      if(ok){
        double sqrtdelta = sqrt(delta);
        retval.uvar_ = halftr + sqrtdelta; // by convention the larger eigenvalue is 'u'
        retval.vvar_ = halftr - sqrtdelta;
        // compute u eigenvector
        retval.udir_ = VEC2(retval.uvar_ - cov_(1,1),cov_(2,2)).Unit();
      }
    }
    if(!ok){
      //fall back to circular errors
      if(cov_(0,0)> cov_(1,1)){
        retval.uvar_ = cov_(0,0);
        retval.vvar_ = cov_(1,1);
        retval.udir_ = VEC2(1.0,0.0);
      } else {
        retval.uvar_ = cov_(1,1);
        retval.vvar_ = cov_(0,0);
        retval.udir_ = VEC2(0.0,1.0);
      }
    }
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

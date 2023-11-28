#include "Offline/GeneralUtilities/inc/TwoDPoint.hh"
#include <stdexcept>
#include <limits>
namespace mu2e {
  using SMAT = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>;
  using VEC2 = ROOT::Math::XYVectorF;

  TwoDPoint::TwoDPoint(SVEC const& pos, SMAT const& cov) : pos_(pos), cov_(cov) {
    // To compute the U direction from the covariance we have to find the eigenvalues and vectors
    double det;
    double halftr = 0.5*cov_.Trace();
    // protect against degenerate matrices
    bool ok = cov_.Det2(det);
    if(ok){
      // check for degenerate eigenvalues
      double delta =  halftr*halftr - det;
      if(delta > std::numeric_limits<float>::min()){
        double sqrtdelta = sqrt(delta);
        double uvar = halftr + sqrtdelta; // by convention the larger eigenvalue is 'u'
        // compute u eigenvector
        udir_ = VEC2(uvar - cov_(1,1),cov_(0,1)).Unit();
      } else {
        udir_ = VEC2(1.0,0.0);
      }
    }
    if(!ok)throw std::invalid_argument( "Unphysical covariance matrix" );
  }

  TwoDPoint::TwoDPoint(VEC2 const& pos,VEC2 const& udir, float uvar, float vvar) : pos_(pos.X(),pos.Y()), udir_(udir) {
    if( uvar < std::numeric_limits<float>::min() || vvar < std::numeric_limits<float>::min())
      throw std::invalid_argument( "Unphysical variance" );
    auto ud = udir.Unit(); // unit vector
    double ucos = ud.X();
    double ucos2 = ucos*ucos;
    double usin = ud.Y();
    double usin2 = usin*usin;
    double cov[3] = {ucos2*uvar + usin2*vvar, ucos*usin*(uvar -vvar),usin2*uvar + ucos2*vvar};
    cov_ = SMAT(cov,cov+3);
  }

  TwoDPoint::TwoDPoint(VEC3 const& pos,VEC3 const& udir, float uvar, float vvar) :
    TwoDPoint(VEC2(pos.X(),pos.Y()), VEC2(udir.X(),udir.Y()),uvar,vvar){}

  float TwoDPoint::uvar() const {
    SVEC uv(udir_.X(),udir_.Y());
    return ROOT::Math::Similarity(uv,cov_);
  }

  float TwoDPoint::vvar() const {
    auto vd = vdir();
    SVEC vv(vd.X(),vd.Y());
    return ROOT::Math::Similarity(vv,cov_);
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

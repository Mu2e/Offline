#include "Offline/TrkHitReco/inc/StereoLine.hh"

namespace mu2e {
  using VEC3 = ROOT::Math::XYZVectorF;
  VEC3 StereoLine::pos(float zval) const {
    double dz = zval - z0_;
    return VEC3(pars_[posx] + pars_[dxdz]*dz, pars_[posy] + pars_[dydz]*dz, zval);
  }

  VEC3 StereoLine::dir() const {
    return VEC3(pars_[dxdz], pars_[dydz],1.0).Unit();
  }

}
std::ostream& operator << (std::ostream& ost, mu2e::StereoLine const& sline) {
  ost << "StereoLine z0 " << sline.z0() << " chisquared " << sline.chisq() << " NDOF " << sline.ndof() << std::endl;
  ost << "Position at z0 " << sline.pos(sline.z0()) << " direction " << sline.dir() << std::endl;
  ost << "Covariance " << sline.cov();
  return ost;
}

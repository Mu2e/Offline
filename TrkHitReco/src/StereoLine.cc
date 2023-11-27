#include "Offline/TrkHitReco/inc/StereoLine.hh"

namespace mu2e {
  using VEC3 = ROOT::Math::XYZVectorF;
  VEC3 StereoLine::pos(float zval) const {
    double dz = zval - z0_;
    return VEC3(pars_[px] + pars_[rx]*dz, pars_[py] + pars_[ry]*dz, zval);
  }

  VEC3 StereoLine::dir() const {
    return VEC3(pars_[rx], pars_[ry],1.0).Unit();
  }

}
std::ostream& operator << (std::ostream& ost, mu2e::StereoLine const& sline) {
  ost << "StereoLine z0 " << sline.z0() << " chisquared " << sline.chisq() << " NDOF " << sline.ndof()
    << " position at z0 " << sline.pos() << " direction " << sline.dir() << " covariance " << sline.cov();
  return ost;
}

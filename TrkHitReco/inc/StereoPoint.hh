//
//  class used in stereo hit reconstruction, essentially a TwoDPoint with a (non-parametric) Z position
//  Original author: Dave Brown (LBNL) 11/15/2023
//
#ifndef TrkHitReco_StereoPoint_hh
#define TrkHitReco_StereoPoint_hh
#include "Offline/GeneralUtilities/inc/TwoDPoint.hh"
namespace mu2e {
  class StereoPoint {
    public:
      using VEC3 = ROOT::Math::XYZVectorF; // spatial only 3D vector
      StereoPoint(TwoDPoint const& pt, double z) : point_(pt), z_(z) {}
      StereoPoint() : z_(0.0) {}
      StereoPoint(VEC3 const& pos,VEC3 const& udir, float uvar, float vvar) : point_(pos, udir, uvar, vvar),z_(pos.Z()) {}
      TwoDPoint const& point() const { return point_; }
      VEC3 pos3() const { return point_.pos3(z_); }
      double z() const { return z_; }
      auto udir() const { return point_.udir(); }
      auto vdir() const { return point_.vdir(); }
      auto uvar() const { return point_.uvar(); }
      auto vvar() const { return point_.vvar(); }
    private:
      TwoDPoint point_;
      double z_; // z position of this point
  };
}
std::ostream& operator << (std::ostream& ost, mu2e::StereoPoint const& pt);
#endif

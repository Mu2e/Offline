//
//  Track position error on a plane reached by a straight step, shared by TrkCrvMatcher and its truth module.
//
#ifndef CRVReco_TrkCrvPlaneError_hh
#define CRVReco_TrkCrvPlaneError_hh
#include "KinKal/General/Vectors.hh"
#include <algorithm>
#include <cmath>

namespace mu2e {
  // Track position error along the in-plane direction a, on the plane with normal w reached by a straight step of
  // signed path length s from a state (position, momentum p) with 6x6 covariance C (position first). With
  // d = p/|p| and P = I - d w^T/(d.w) (move along d onto the plane), the gradient of a.(position on the plane) is
  // P^T a for the position block and (s/|p|)(I - d d^T) P^T a for the momentum block. s = 0 gives the error of the
  // state itself, measured in the plane rather than in 3D. Returns 0 for a track parallel to the plane.
  template<class COV> double trkCrvPlaneError(COV const& C, KinKal::VEC3 const& p, KinKal::VEC3 const& a,
                                              KinKal::VEC3 const& w, double s) {
    double pm = std::sqrt(p.Mag2());
    if(pm <= 0) return 0;
    KinKal::VEC3 d = p / pm;
    double dw = d.Dot(w);
    if(std::fabs(dw) < 1e-3) return 0;
    KinKal::VEC3 gp = a - w * (d.Dot(a) / dw);
    KinKal::VEC3 gm = (s / pm) * (gp - d * d.Dot(gp));
    double g[6] = {gp.X(), gp.Y(), gp.Z(), gm.X(), gm.Y(), gm.Z()};
    double var = 0;
    for(int i = 0; i < 6; ++i) for(int j = 0; j < 6; ++j) var += g[i] * C(i, j) * g[j];
    return std::sqrt(std::max(var, 0.0));
  }
}
#endif

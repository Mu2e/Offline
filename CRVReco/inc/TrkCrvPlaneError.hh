//
//  Temporary stand-in for a KinKal projection of a track state's spatial covariance onto a plane; replace with the
//  KinKal implementation when it exists (both modules call only trkCrvPlaneCovariance).
//
#ifndef CRVReco_TrkCrvPlaneError_hh
#define CRVReco_TrkCrvPlaneError_hh
#include "KinKal/General/Vectors.hh"
#include "Math/SMatrix.h"
#include <cmath>

namespace mu2e {
  using TrkCrvCov2 = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>;

  // Covariance of the track position on the plane with normal w, in the in-plane axes (a, b), after a straight step of
  // signed path length s from a state (position, momentum p) with 6x6 covariance C (position first). With d = p/|p|
  // and P = I - d w^T/(d.w) (move along d onto the plane), the gradient of e.(position on the plane) is P^T e for the
  // position block and (s/|p|)(I - d d^T) P^T e for the momentum block. s = 0 gives the covariance of the state itself,
  // measured in the plane. Returns a zero matrix for a track parallel to the plane.
  template<class COV> TrkCrvCov2 trkCrvPlaneCovariance(COV const& C, KinKal::VEC3 const& p, KinKal::VEC3 const& a,
                                                       KinKal::VEC3 const& b, KinKal::VEC3 const& w, double s) {
    TrkCrvCov2 out;
    double pm = std::sqrt(p.Mag2());
    if(pm <= 0) return out;
    KinKal::VEC3 d = p / pm;
    double dw = d.Dot(w);
    if(std::fabs(dw) < 1e-3) return out;
    double g[2][6];
    KinKal::VEC3 const* axes[2] = {&a, &b};
    for(int k = 0; k < 2; ++k) {
      KinKal::VEC3 gp = *axes[k] - w * (d.Dot(*axes[k]) / dw);
      KinKal::VEC3 gm = (s / pm) * (gp - d * d.Dot(gp));
      double gk[6] = {gp.X(), gp.Y(), gp.Z(), gm.X(), gm.Y(), gm.Z()};
      for(int i = 0; i < 6; ++i) g[k][i] = gk[i];
    }
    for(int k = 0; k < 2; ++k)
      for(int l = k; l < 2; ++l) {
        double v = 0;
        for(int i = 0; i < 6; ++i) for(int j = 0; j < 6; ++j) v += g[k][i] * C(i, j) * g[l][j];
        out(k, l) = v;
      }
    return out;
  }
}
#endif

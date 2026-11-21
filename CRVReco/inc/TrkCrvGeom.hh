//
//  CRV geometry used by the track-CRV association and its truth: the sector planes (KinKalGeom CRV rectangles, per
//  CRSScintillatorShield index) and the dominant shield of a coincidence cluster.
//  Sector frame: u along the bars (uDirection), v across the bars (vDirection), w the sector normal, all in the
//  detector frame; TrkCrvMatch residuals and covariances are expressed in (u, v).
//
#ifndef CRVReco_TrkCrvGeom_hh
#define CRVReco_TrkCrvGeom_hh
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/RecoDataProducts/inc/KalIntersection.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Geometry/Plane.hh"
#include "KinKal/Geometry/Rectangle.hh"
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace mu2e {
  class CosmicRayShield;
  class CrvCoincidenceCluster;

  struct TrkCrvSector {
    bool valid = false;                         // the shield has a KinKalGeom CRV rectangle
    std::string name;                           // CRSScintillatorShield name
    SurfaceId sid;                              // its SurfaceId (unknown if the name has none)
    std::shared_ptr<KinKal::Rectangle> rect;    // mid-plane rectangle, detector frame
    KinKal::VEC3 center, u, v, w;               // rectangle centre and axes
    double uHalf = 0, vHalf = 0;                // half lengths along u and v [mm]
  };

  // one entry per CRSScintillatorShield, in shield order (the sector_ index of TrkCrvMatch); call once per run
  std::vector<TrkCrvSector> trkCrvSectors();
  // the shield with the largest PE sum among the cluster's pulses, -1 if none; multiSector set if pulses span shields
  int trkCrvDominantShield(CrvCoincidenceCluster const& cluster, CosmicRayShield const& crs, bool* multiSector = nullptr);

  constexpr double trkCrvMinCosToNormal = 1e-3;   // |direction . plane normal| below which a track is parallel to the plane

  // Covariance of the track position on the plane through point with the sector axes (u, v), for the straight line from a KalSeed
  // intersection reaching that plane at time (the intersection state's covariance carried over the straight step; time = ki.time()
  // gives the covariance of the state itself, measured in the plane). 0 for a line parallel to the plane.
  inline KinKal::PMAT trkCrvPlaneCovariance(KalIntersection const& ki, TrkCrvSector const& g, KinKal::VEC3 const& point, double time) {
    if(std::fabs(ki.momentum3().Unit().Dot(g.w)) < trkCrvMinCosToNormal) return KinKal::PMAT();
    return ki.kinematicLine().planeCovariance(time, KinKal::Plane(g.w, g.u, point), trkCrvMinCosToNormal);
  }

  // Covariance of the track position at a KalSeed intersection with the sector, in the sector axes (u, v): the KinKal projected
  // covariance of the intersection (filled by the extrapolation, rotated from the intersection's local plane axes) if present,
  // otherwise that of the intersection state measured in the plane.
  inline KinKal::PMAT trkCrvXingCovariance(KalIntersection const& ki, TrkCrvSector const& g) {
    auto const& inter = ki.intersection();
    if(!inter.hascov_) return trkCrvPlaneCovariance(ki, g, ki.position3(), ki.time());
    ROOT::Math::SMatrix<double,2,2> rot;
    KinKal::VEC3 const* axes[2] = {&g.u, &g.v};
    KinKal::VEC3 const* idirs[2] = {&inter.udir_, &inter.vdir_};
    for(int i = 0; i < 2; ++i) for(int j = 0; j < 2; ++j) rot(i, j) = axes[i]->Dot(*idirs[j]);
    return ROOT::Math::Similarity(rot, inter.pcov_);
  }
}
#endif

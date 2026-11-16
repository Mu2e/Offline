//
//  CRV geometry used by the track-CRV association and its truth: the sector planes (KinKalGeom CRV rectangles, per
//  CRSScintillatorShield index) and the dominant shield of a coincidence cluster.
//  Sector frame: u along the bars (uDirection), v across the bars (vDirection), w the sector normal, all in the
//  detector frame; TrkCrvMatch residuals and covariances are expressed in (u, v).
//
#ifndef CRVReco_TrkCrvGeom_hh
#define CRVReco_TrkCrvGeom_hh
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Geometry/Rectangle.hh"
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
}
#endif

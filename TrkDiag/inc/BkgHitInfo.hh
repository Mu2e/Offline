#ifndef TrkDiag_BkgHitInfo_HH
#define TrkDiag_BkgHitInfo_HH
#include "Offline/TrkDiag/inc/StrawHitInfo.hh"
namespace mu2e {
  // extend StrawHitInfo struct with additional info on background hits
  struct BkgHitInfo : public StrawHitInfo {
    float _gdist = 0; // generalized distance
    float _rerr = 0; // estimated error on transverse distance
    int _index = 0; // index into StrawHit collection
    XYZVectorF _rpos; // relative position of hit WRT cluster center
    bool _active = false; // hit is active in cluster
    bool _cbkg = false; // hit is classified backg
  };
}
#endif

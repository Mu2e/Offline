//
//  MC truth for a TrkCrvMatch, by SimParticle identity: the true CRV sector crossings of the track's particle and its
//  ancestors, the truth of each reco crossing, and the relation of each candidate cluster to the track's particle.
//  One TrkCrvMatchMC per TrkCrvMatch, in the same order; candidates and reco crossings are parallel to the reco record.
//
#ifndef MCDataProducts_TrkCrvMatchMC_hh
#define MCDataProducts_TrkCrvMatchMC_hh
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/RecoDataProducts/inc/TrkCrvMatch.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <cstdint>
#include <vector>

namespace mu2e {

  struct TrkCrvTruth {
    // relation of a cluster's energy deposits to the track's particle, best first
    enum Relation : uint8_t { SameParticle = 0, SameFamily, SamePrimary, Unrelated, Background, NoMC };
  };

  // one true crossing of a CRV sector by the track's particle (generation 0) or one of its ancestors
  struct TrkCrvXingMC {
    SurfaceId sid_;              // sector surface, as the reco crossings
    int16_t   sector_ = -1;      // CRSScintillatorShield index
    uint8_t   generation_ = 0;   // 0 = the track's particle, 1 = its parent, ...
    art::Ptr<SimParticle> sim_;
    bool      geom_ = false;     // the particle's trajectory crosses the sector mid-plane rectangle (gap-inclusive)
    XYZVectorF pos_;             // trajectory crossing of the mid-plane, detector frame [mm] (geom_ only)
    XYZVectorF dir_;             // unit direction of the trajectory chord at the crossing (geom_ only)
    float     time_ = 0;         // time of the mid-plane crossing [ns] (geom_ only)
    float     mom_ = 0;          // momentum magnitude at the mid-plane crossing [MeV/c] (geom_ only)
    float     dEdgeU_ = 0;       // distance of that crossing inside the sector edge along the bars [mm] (geom_ only)
    float     dEdgeV_ = 0;       // distance of that crossing inside the sector edge across the bars [mm] (geom_ only)
    float     edep_ = 0;         // visible energy the particle deposited in the sector's bars [MeV]
    uint8_t   nLayers_ = 0;      // layers of the sector with a deposit from the particle
    XYZVectorF depPos_;          // energy-weighted position of those deposits, detector frame [mm]
    int16_t   trueCluster_ = -1; // cluster (index in the matched collection) with the most energy from the particle here
    int16_t   trueCand_ = -1;    // index into TrkCrvMatch::cands_ of trueCluster_, -1 if it is not a candidate
    int16_t   recoXing_ = -1;    // index of the reco crossing of the same sector, -1 = not extrapolated there
  };

  struct TrkCrvCandidateMC {
    uint8_t rel_ = TrkCrvTruth::NoMC;  // a TrkCrvTruth::Relation
    float edepSame_ = 0;         // cluster energy from the track's particle and its descendants [MeV]
    float edepFamily_ = 0;       // cluster energy from the track particle's ancestors and their other descendants [MeV]
    float edepTotal_ = 0;        // all cluster energy with MC information [MeV]
  };

  struct TrkCrvMatchMC {
    art::Ptr<TrkCrvMatch> match_;           // the reco record this truth belongs to
    art::Ptr<SimParticle> trackSim_;        // the KalSeedMC primary contributor
    std::vector<TrkCrvXingMC> xings_;        // true crossings
    std::vector<int16_t> recoXingTruth_;     // per reco crossing: index into xings_ of the same sector, -1 = no true crossing
    std::vector<int16_t> recoXingNeighbor_;  // per reco crossing without one: true crossing of a coplanar sector by the track's particle, -1 = none
    std::vector<TrkCrvCandidateMC> cands_;   // parallel to TrkCrvMatch::cands_
    auto const& xings() const { return xings_; }
    auto const& candidates() const { return cands_; }
  };
  using TrkCrvMatchMCCollection = std::vector<TrkCrvMatchMC>;
}
#endif

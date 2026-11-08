//
//  MC truth for a TrkCrvMatch: the true CRV sector crossings of the track's particle and its ancestors, the truth of
//  each reco crossing, and the relation of each candidate cluster to the track's particle.
//  One TrkCrvMatchMC per TrkCrvMatch, in the same order; candidates and reco crossings are parallel to the reco record.
//
#ifndef MCDataProducts_TrkCrvMatchMC_hh
#define MCDataProducts_TrkCrvMatchMC_hh
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/MCDataProducts/inc/KalSeedMC.hh"
#include "Offline/MCDataProducts/inc/MCRelationship.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/RecoDataProducts/inc/TrkCrvMatch.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <cstdint>
#include <vector>

namespace mu2e {

  // one true crossing of a CRV sector by the track's particle or one of its ancestors (sim_); the
  // position, direction, time and momentum are those of sim_'s MC trajectory at the sector mid-plane
  struct TrkCrvXingMC {
    SurfaceId sid_;              // sector surface, as the reco crossings
    int16_t   sector_ = -1;      // CRSScintillatorShield index
    MCRelationship rel_;         // sim_ relative to the track's particle: same, mother (its parent) or umother (further
                                 // up); removal() counts the generations, not counting simulation-stage copies
    art::Ptr<SimParticle> sim_;  // the crossing particle
    bool      geom_ = false;     // sim_'s trajectory crosses the sector mid-plane rectangle (gap-inclusive)
    XYZVectorF pos_;             // sim_'s trajectory at the mid-plane, detector frame [mm] (geom_ only)
    XYZVectorF dir_;             // unit direction of sim_'s trajectory there (geom_ only)
    float     time_ = 0;         // sim_'s time at the mid-plane [ns] (geom_ only)
    float     mom_ = 0;          // sim_'s momentum magnitude at the mid-plane [MeV/c] (geom_ only)
    float     dEdgeU_ = 0;       // distance of that point inside the sector edge along the bars [mm] (geom_ only)
    float     dEdgeV_ = 0;       // distance of that point inside the sector edge across the bars [mm] (geom_ only)
    float     edep_ = 0;         // visible energy sim_ deposited in the sector's bars [MeV] (CrvSteps; 0 without them)
    uint8_t   nLayers_ = 0;      // layers of the sector with a deposit from sim_
    XYZVectorF depPos_;          // energy-weighted position of those deposits, detector frame [mm]
    int16_t   trueCluster_ = -1; // cluster (index in the matched collection) with the most energy from sim_ here
    int16_t   trueCand_ = -1;    // index into TrkCrvMatch::cands_ of trueCluster_, -1 if it is not a candidate
    int16_t   recoXing_ = -1;    // index of the reco crossing of the same sector, -1 = not extrapolated there
  };

  struct TrkCrvCandidateMC {
    MCRelationship rel_;         // closest relation of the cluster's deposit particles to the track's particle (none: unrelated or no MC)
    bool  hasMC_ = false;        // the cluster has MC information (false: background or no MC)
    float edepSame_ = 0;         // cluster energy from the track's particle (relation same) [MeV]
    float edepRelated_ = 0;      // cluster energy from any related particle, the track's particle included [MeV]
    float edepTotal_ = 0;        // all cluster energy with MC information [MeV]
  };

  struct TrkCrvMatchMC {
    art::Ptr<TrkCrvMatch> match_;           // the reco record this truth belongs to
    art::Ptr<KalSeedMC> kalSeedMC_;         // the track's KalSeedMC; its primary contributor is the track's particle
    bool hasDeposits_ = false;              // CrvSteps were available: edep_, nLayers_ and depPos_ are filled
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

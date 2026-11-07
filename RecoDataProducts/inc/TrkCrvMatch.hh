//
//  Association of one reconstructed track (KalSeed) with CRV coincidence clusters: the track's CRV crossings and
//  every cluster inside the loose time or geometry storage windows, with signed residuals, ranks and flags.
//  One TrkCrvMatch per input KalSeed, in the KalSeed collection's order.
//  Sector frame of the residuals and covariances: u along the bars, v across the bars (the uDirection and vDirection
//  of the sector's KinKalGeom CRV rectangle, detector frame); trkCrvSectors() in CRVReco/inc/TrkCrvGeom.hh returns
//  them per sector_ index.
//
#ifndef RecoDataProducts_TrkCrvMatch_hh
#define RecoDataProducts_TrkCrvMatch_hh
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "Math/SMatrix.h"
#include <cmath>
#include <cstdint>
#include <vector>

namespace mu2e {

  // clusters a track fit used as constraints, for a future fit with CRV constraints (nothing produces it yet); the
  // matcher flags those candidates usedInFit
  using TrkCrvFitClusterAssns = art::Assns<KalSeed, CrvCoincidenceCluster>;
  // symmetric covariance in the sector plane, (u, v) order [mm^2]
  using TrkCrvPlaneCov = ROOT::Math::SMatrix<float,2,2,ROOT::Math::MatRepSym<float,2>>;

  // one CRV layer on the track's straight line through a crossed sector, with the best pulse near the predicted bar
  struct TrkCrvLayerProbe {
    int16_t bar_ = -1;           // predicted CRSScintillatorBarIndex in this layer, -1 if the line misses the layer
    int8_t  dBar_ = INT8_MAX;    // bar offset of the selected pulse from the predicted bar, INT8_MAX if none in the windows
    float   uPred_ = 0;          // predicted position along the bar from its centre, along the Mu2e length axis [mm]
    float   tPred_ = 0;          // track time at the layer [ns]
    float   pe_ = 0;             // PEs of the in-window pulses of the selected bar (both ends); kept since unclustered pulses are not persisted
    float   dt_ = 0;             // pulse time - predicted pulse time (track time plus fibre propagation) [ns]
    bool    inCluster_ = false;  // an in-window pulse of the selected bar belongs to a coincidence cluster
    float   uTdiv_ = NAN;        // position along the bar from the two-end time difference [mm], NaN if one-ended
  };

  // one CRV sector crossing of the track's extrapolation (a CRV KalIntersection)
  struct TrkCrvXing {
    SurfaceId sid_;              // CRV sector surface (CRV_EX ... CRV_M8)
    int16_t   inter_ = -1;       // index into KalSeed::intersections()
    int16_t   sector_ = -1;      // CRSScintillatorShield index
    XYZVectorF pos_;             // crossing position on the sector mid-plane, detector frame [mm]
    XYZVectorF dir_;             // unit direction of the track at the crossing
    float     time_ = 0;         // crossing time [ns]
    TrkCrvPlaneCov cov_;         // track position covariance in the sector plane (u, v), error model applied [mm^2]
    float     sigT_ = 0;         // track time error [ns]
    float     mom_ = 0;          // track momentum magnitude at the crossing [MeV/c]
    float     momErr_ = 0;       // its error [MeV/c]
    bool      inBounds_ = false; // crossing inside the sector rectangle
    float     dEdgeU_ = 0;       // distance inside the sector edge along the bars [mm], negative outside
    float     dEdgeV_ = 0;       // distance inside the sector edge across the bars [mm], negative outside
    float     edgeSig_ = 0;      // min(dEdgeU_/sigU(), dEdgeV_/sigV()): edge distance in units of the track error
    int16_t   best_ = -1;        // index into TrkCrvMatch::cands_ of the best geometric candidate in this crossing's sector, -1 if none
    std::vector<TrkCrvLayerProbe> layers_; // one per layer of the sector, in layer order; empty unless the layer probe ran
    float sigU() const { return std::sqrt(cov_(0,0)); }
    float sigV() const { return std::sqrt(cov_(1,1)); }
  };

  // one CRV coincidence cluster associated with the track
  struct TrkCrvCandidate {
    enum Flag : uint16_t { timeCand = 1, geomCand = 2, singleEnded = 4, multiSector = 8, bestForTrack = 16,
                           bestForCluster = 32, exclusive = 64, outOfBounds = 128, usedInFit = 256,
                           sameSector = 512, bestInSector = 1024 };
    art::Ptr<CrvCoincidenceCluster> cluster_;
    int16_t  xing_ = -1;         // reference crossing (nearest in path, projected to the cluster plane), -1 = time only
    int16_t  sector_ = -1;       // dominant shield of the cluster (largest PE sum)
    float    dtRef_ = 0;         // t_track(reference surface) - cluster avgHitTime [ns]
    float    dtRefToF_ = 0;      // dtRef - |x_ref - x_cluster|/c (beta = 1 straight path) [ns]
    float    dtXing_ = 0;        // cluster avgHitTime - t_track at the cluster plane [ns]
    float    du_ = 0;            // cluster - track along the bars, at the cluster plane [mm]
    float    dv_ = 0;            // cluster - track across the bars, at the cluster plane [mm]
    TrkCrvPlaneCov cov_;         // track covariance of (du, dv), carried from the reference crossing over the step, error model applied [mm^2]
    float    chi2_ = -1;         // space + time chi2 (geometric candidates only), with the full (u, v) covariance
    uint8_t  ndof_ = 0;          // chi2 degrees of freedom: 2 for single-ended clusters (no along-bar term), else 3
    uint8_t  nShare_ = 0;        // tracks of this collection holding this cluster as a candidate
    uint16_t flags_ = 0;
    bool hasFlag(Flag f) const { return flags_ & f; }
    void setFlag(Flag f) { flags_ |= f; }
    float sigU() const { return std::sqrt(cov_(0,0)); }
    float sigV() const { return std::sqrt(cov_(1,1)); }
  };

  struct TrkCrvMatch {
    enum Probe : uint8_t { probeOff = 0, probeNoPulses, probeClusteredPulses, probeAllPulses };
    art::Ptr<KalSeed> track_;
    SurfaceId refSurface_;       // surface that defined tRef_ (unknown = KalSeed t0); other surface times are in the KalSeed
    float tRef_ = 0;             // reference track time used by dtRef [ns]
    float timeOffset_ = 0;       // CRV - tracker time offset applied to the cluster times [ns]
    uint8_t probe_ = probeOff;   // pulses the layer probe saw: none (predicted bars only), clustered only, or all
    std::vector<TrkCrvXing> xings_;
    std::vector<TrkCrvCandidate> cands_; // geometric candidates by chi2, then time-only by |dtRef - centre|
    int16_t bestGeom_ = -1;      // lowest-chi2 geometric candidate, -1 if none
    int16_t bestTime_ = -1;      // candidate closest in dtRef to the configured centre, -1 if none
    auto const& xings() const { return xings_; }
    auto const& candidates() const { return cands_; }
  };
  using TrkCrvMatchCollection = std::vector<TrkCrvMatch>;
}
#endif

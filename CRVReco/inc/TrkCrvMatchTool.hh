//
//  Track-CRV matching shared by the TrkCrvMatcher producer and any fit that uses CRV clusters as constraints, so that
//  both apply one definition of a match. See TrkCrvMatch.hh for the output.
//
#ifndef CRVReco_TrkCrvMatchTool_hh
#define CRVReco_TrkCrvMatchTool_hh
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkCrvMatch.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "KinKal/General/Vectors.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace mu2e {
  class CosmicRayShield;
  class DetectorSystem;

  class TrkCrvMatchTool {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Sequence<std::string> refSurfaces{Name("ReferenceSurfaces"), Comment("track reference-time surfaces in order of preference; KalSeed t0 if none present"), std::vector<std::string>{"TT_Front","TT_Mid"}};
        fhicl::Sequence<double> dtRefWindow{Name("DtRefWindow"), Comment("storage window on t_track(ref) - t_cluster [ns]"), std::vector<double>{-350.0, 350.0}};
        fhicl::Atom<double> dtRefCenter{Name("DtRefCenter"), Comment("dtRef value the time ranking is centred on [ns]"), 75.0};
        fhicl::Sequence<double> dtXingWindow{Name("DtXingWindow"), Comment("storage window on t_cluster - t_track at the cluster plane [ns]"), std::vector<double>{-100.0, 100.0}};
        fhicl::Atom<double> maxDv{Name("MaxDv"), Comment("storage window on |across-bar residual| [mm]"), 500.0};
        fhicl::Atom<double> maxDu{Name("MaxDu"), Comment("storage window on |along-bar residual| for two-ended clusters [mm]"), 3000.0};
        fhicl::Atom<bool> doGeometric{Name("DoGeometric"), Comment("form geometric candidates from CRV KalIntersections"), true};
        fhicl::Atom<bool> requireSameSector{Name("RequireSameSector"), Comment("project only from a crossing of the cluster's own sector"), false};
        fhicl::Sequence<double> sigmaV{Name("SigmaV"), Comment("CRV across-bar resolution [mm], per cluster sector type 1-4 (EX, T1, T2, muon taggers; last entry for higher types); MC fit on 50 MDC2025bb files"), std::vector<double>{7.0, 13.0, 7.0, 17.0}};
        fhicl::Sequence<double> sigmaU{Name("SigmaU"), Comment("CRV along-bar resolution of two-ended clusters [mm], per sector type"), std::vector<double>{300.0}};
        fhicl::Sequence<double> sigmaT{Name("SigmaT"), Comment("CRV time resolution of two-ended clusters [ns], per sector type"), std::vector<double>{5.0}};
        fhicl::Atom<double> sigmaTSingleEnded{Name("SigmaTSingleEnded"), Comment("CRV time resolution of single-ended clusters [ns]"), 15.0};
        fhicl::Atom<double> trackErrorScale{Name("TrackErrorScale"), Comment("scale on the propagated track position error; sigma = sqrt((scale*err)^2 + TrackPosErrQuad^2)"), 1.5};
        fhicl::Atom<double> trackPosErrQuad{Name("TrackPosErrQuad"), Comment("added in quadrature to the scaled track position error [mm]"), 25.0};
        fhicl::Atom<double> minTrackTimeErr{Name("MinTrackTimeErr"), Comment("floor on the (unscaled) track time error [ns]"), 1.0};
        fhicl::Atom<double> timeOffset{Name("TimeOffset"), Comment("added to cluster and pulse times [ns]; test knob, data offsets belong in conditions"), 0.0};
        fhicl::Atom<bool> keepAll{Name("KeepAllClusters"), Comment("ignore every storage window (commissioning)"), false};
        fhicl::Atom<bool> layerProbe{Name("LayerProbe"), Comment("probe each layer of every crossed sector for pulses near the predicted bar; needs the reco pulses"), true};
        fhicl::Atom<int> probeBarWindow{Name("ProbeBarWindow"), Comment("layer probe: bars searched on each side of the predicted bar"), 1};
        fhicl::Sequence<double> probeDtWindow{Name("ProbeDtWindow"), Comment("layer probe: window on pulse time - predicted pulse time [ns]"), std::vector<double>{-50.0, 50.0}};
        fhicl::Atom<double> fiberSignalSpeed{Name("FiberSignalSpeed"), Comment("layer probe: signal speed in the fibres [mm/ns], as in CrvCoincidenceFinder"), 140.0};
        fhicl::Atom<double> pulseDelay{Name("ProbePulseDelay"), Comment("layer probe: electronics and scintillation delay subtracted from pulse times [ns], as the CrvCoincidenceFinder timeOffset"), 34.5};
      };
      explicit TrkCrvMatchTool(Config const& conf);

      void beginRun();   // caches the CRV geometry
      // per-event inputs; pulses may be null (no layer probe)
      void setEvent(std::vector<art::Ptr<CrvCoincidenceCluster>> const& clusters, CrvRecoPulseCollection const* pulses);
      // one track; usedInFit lists the clusters its fit used as constraints
      TrkCrvMatch match(art::Ptr<KalSeed> const& seed, std::vector<art::Ptr<CrvCoincidenceCluster>> const& usedInFit = {}) const;
      // cross-track pass over one collection: nShare, bestForCluster, exclusive
      void resolve(TrkCrvMatchCollection& matches) const;
      bool probeActive() const { return layerProbe_; }

    private:
      struct SectorGeom {         // CRV sector mid-plane, detector frame
        bool valid = false;
        KinKal::VEC3 center, u, v, w;
        double uHalf = 0, vHalf = 0;
      };
      struct ClusterInfo {        // per-event cluster quantities
        KinKal::VEC3 pos;         // detector frame
        double time = 0;
        int shield = -1;
        int sectorType = 0;
        bool multiSector = false;
        bool singleEnded = false;
      };
      double resolution(std::vector<double> const& tab, int sectorType) const;
      double trkPosErr(double err) const { return std::sqrt(trkErrScale_ * trkErrScale_ * err * err + posErrQuad_ * posErrQuad_); }
      void probe(TrkCrvXing& x, double speed) const;
      static bool sameCluster(art::Ptr<CrvCoincidenceCluster> const& a, art::Ptr<CrvCoincidenceCluster> const& b);

      std::vector<SurfaceId> refSurfaces_;
      double dtRefLo_, dtRefHi_, dtRefCenter_, dtXingLo_, dtXingHi_, maxDv_, maxDu_;
      bool doGeometric_, requireSameSector_;
      std::vector<double> sigmaV_, sigmaU_, sigmaT_;
      double sigmaTSingle_, trkErrScale_, posErrQuad_, minTimeErr_, timeOffset_;
      bool keepAll_, layerProbe_;
      int probeBarWindow_;
      double probeDtLo_, probeDtHi_, fiberSpeed_, pulseDelay_;

      CosmicRayShield const* crs_ = nullptr;
      DetectorSystem const* det_ = nullptr;
      std::map<SurfaceId,int> sidToShield_;
      std::vector<SectorGeom> sectors_;   // indexed by shield

      std::vector<art::Ptr<CrvCoincidenceCluster>> clusters_;
      std::vector<ClusterInfo> cinfo_;
      CrvRecoPulseCollection const* pulses_ = nullptr;
      std::unordered_map<int, std::vector<size_t>> pulsesByBar_;   // bar index -> pulses
      std::set<std::tuple<int,int,double>> clusteredPulses_;        // (bar, SiPM, time) of the pulses in clusters

      static constexpr double cLight_ = 299.792458; // mm/ns
  };
}
#endif

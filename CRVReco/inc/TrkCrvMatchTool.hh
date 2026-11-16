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
#include "Offline/CRVReco/inc/TrkCrvGeom.hh"
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
      // physics settings have no C++ default: CRVReco/fcl/prolog_v12.fcl (TrkCrvMatcher) sets them with their sources
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Sequence<std::string> refSurfaces{Name("ReferenceSurfaces"), Comment("surfaces whose KalSeed intersection gives the track reference time; KalSeed t0 if none present")};
        fhicl::Atom<bool> earliestReference{Name("EarliestReference"), Comment("true: the earliest intersection on any of ReferenceSurfaces (the tracker entry); false: the first surface in list order present, nearest in time to t0")};
        fhicl::Sequence<double> dtRefWindow{Name("DtRefWindow"), Comment("storage window on t_track(ref) - t_cluster [ns]")};
        fhicl::Atom<double> dtRefCenter{Name("DtRefCenter"), Comment("dtRef value the time ranking is centred on [ns]")};
        fhicl::Sequence<double> dtXingWindow{Name("DtXingWindow"), Comment("storage window on t_cluster - t_track at the cluster plane [ns]")};
        fhicl::Atom<double> maxDv{Name("MaxDv"), Comment("storage window on |across-bar residual| [mm]")};
        fhicl::Atom<double> maxDu{Name("MaxDu"), Comment("storage window on |along-bar residual| for two-ended clusters [mm]")};
        fhicl::Atom<bool> doGeometric{Name("DoGeometric"), Comment("form geometric candidates from CRV KalIntersections")};
        fhicl::Atom<bool> requireSameSector{Name("RequireSameSector"), Comment("project only from a crossing of the cluster's own sector (A/B tests)"), false};
        fhicl::Sequence<double> sigmaV{Name("SigmaV"), Comment("CRV across-bar resolution [mm], indexed by the cluster's sector type (1 first; the last entry serves higher types); sector type numbering depends on the geometry")};
        fhicl::Sequence<double> sigmaU{Name("SigmaU"), Comment("CRV along-bar resolution of two-ended clusters [mm], by sector type")};
        fhicl::Sequence<double> sigmaT{Name("SigmaT"), Comment("CRV time resolution of two-ended clusters [ns], by sector type")};
        fhicl::Atom<double> sigmaTSingleEnded{Name("SigmaTSingleEnded"), Comment("CRV time resolution of single-ended clusters [ns]")};
        fhicl::Atom<double> trackErrorScale{Name("TrackErrorScale"), Comment("scale k on the propagated track position covariance: k^2 C + TrackPosErrQuad^2 I")};
        fhicl::Atom<double> trackPosErrQuad{Name("TrackPosErrQuad"), Comment("added in quadrature to the scaled track position errors [mm]")};
        fhicl::Atom<double> minTrackTimeErr{Name("MinTrackTimeErr"), Comment("floor on the (unscaled) track time error [ns]")};
        fhicl::Atom<double> timeOffset{Name("TimeOffset"), Comment("added to cluster and pulse times [ns]; test knob, data offsets belong in conditions"), 0.0};
        fhicl::Atom<bool> keepAll{Name("KeepAllClusters"), Comment("ignore every storage window (commissioning)"), false};
        fhicl::Atom<bool> layerProbe{Name("LayerProbe"), Comment("probe each layer of every crossed sector for pulses near the predicted bar; needs the reco pulses")};
        fhicl::Atom<bool> probeAllPulses{Name("ProbeAllPulses"), Comment("the pulse collection holds every reco pulse (true) or only the clustered ones, e.g. SelectReco (false)")};
        fhicl::Atom<int> probeBarWindow{Name("ProbeBarWindow"), Comment("layer probe: bars searched on each side of the predicted bar")};
        fhicl::Sequence<double> probeDtWindow{Name("ProbeDtWindow"), Comment("layer probe: window on pulse time - predicted pulse time [ns]")};
        fhicl::Atom<double> fiberSignalSpeed{Name("FiberSignalSpeed"), Comment("layer probe: signal speed in the fibres [mm/ns]")};
        fhicl::Atom<double> pulseDelay{Name("ProbePulseDelay"), Comment("layer probe: electronics and scintillation delay subtracted from pulse times [ns]")};
      };
      explicit TrkCrvMatchTool(Config const& conf);

      void beginRun();   // caches the CRV geometry
      // per-event inputs; pulses may be null (the layer probe then records predicted bars only)
      void setEvent(std::vector<art::Ptr<CrvCoincidenceCluster>> const& clusters, CrvRecoPulseCollection const* pulses);
      // one track; usedInFit lists the clusters its fit used as constraints
      TrkCrvMatch match(art::Ptr<KalSeed> const& seed, std::vector<art::Ptr<CrvCoincidenceCluster>> const& usedInFit = {}) const;
      // cross-track pass over one collection: nShare, bestForCluster, exclusive
      void resolve(TrkCrvMatchCollection& matches) const;
      bool probeActive() const { return layerProbe_; }
      std::vector<TrkCrvSector> const& sectors() const { return sectors_; }   // per shield index, after beginRun

    private:
      struct ClusterInfo {        // per-event cluster quantities
        KinKal::VEC3 pos;         // detector frame
        double time = 0;
        int shield = -1;
        int sectorType = 0;
        bool multiSector = false;
        bool singleEnded = false;
      };
      double resolution(std::vector<double> const& tab, int sectorType) const;
      TrkCrvPlaneCov errorModel(ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>> const& cov) const;   // k^2 C + a^2 I
      void probe(TrkCrvXing& x, double speed) const;
      static bool sameCluster(art::Ptr<CrvCoincidenceCluster> const& a, art::Ptr<CrvCoincidenceCluster> const& b);

      std::vector<SurfaceId> refSurfaces_;
      bool earliestRef_;
      double dtRefLo_, dtRefHi_, dtRefCenter_, dtXingLo_, dtXingHi_, maxDv_, maxDu_;
      bool doGeometric_, requireSameSector_;
      std::vector<double> sigmaV_, sigmaU_, sigmaT_;
      double sigmaTSingle_, trkErrScale_, posErrQuad_, minTimeErr_, timeOffset_;
      bool keepAll_, layerProbe_, probeAllPulses_;
      int probeBarWindow_;
      double probeDtLo_, probeDtHi_, fiberSpeed_, pulseDelay_;

      CosmicRayShield const* crs_ = nullptr;
      DetectorSystem const* det_ = nullptr;
      std::map<SurfaceId,int> sidToShield_;
      std::vector<TrkCrvSector> sectors_;   // indexed by shield

      std::vector<art::Ptr<CrvCoincidenceCluster>> clusters_;
      std::vector<ClusterInfo> cinfo_;
      CrvRecoPulseCollection const* pulses_ = nullptr;
      std::unordered_map<int, std::vector<size_t>> pulsesByBar_;   // bar index -> pulses
      std::set<std::tuple<int,int,double>> clusteredPulses_;        // (bar, SiPM, time) of the pulses in clusters
  };
}
#endif

//
//  Associate reconstructed tracks (KalSeeds) with CRV coincidence clusters, one TrkCrvMatch per KalSeed, using
//  TrkCrvMatchTool:
//  - time: every cluster with t_track(reference) - t_cluster inside DtRefWindow is a time candidate; the reference
//    is by default the track's entry into the tracker (earliest of TT_Front, TT_Back, TT_Outer);
//  - space: from the track's CRV crossing nearest to the cluster's sector plane, a straight-line projection onto that
//    plane gives du, dv and dtXing, and a space+time chi2 with the full (u, v) covariance; clusters inside MaxDu,
//    MaxDv and DtXingWindow are geometric candidates; a straight-line fit with no CRV crossing (its line passes
//    outside every sector) is projected instead from its intersection nearest the cluster's sector onto that
//    sector's extended plane, recorded as a crossing with inBounds false;
//  - every candidate is kept: a track can have several, and a cluster can be a candidate of several tracks. They are
//    ranked by chi2, then by time, with bestGeom and bestTime per track and bestInSector per sector; across the
//    collection nShare, bestForCluster and a greedy one-to-one 'exclusive' assignment resolve shared clusters.
//  Storage windows are loose on purpose: selections belong to the analysis.
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/TableFragment.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/CRVReco/inc/TrkCrvMatchTool.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkCrvMatch.hh"

#include <map>
#include <memory>
#include <vector>

namespace mu2e {

  class TrkCrvMatcher : public art::EDProducer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> kalSeedTag{Name("KalSeedCollection"), Comment("KalSeed collection; one TrkCrvMatch per KalSeed")};
        fhicl::Atom<art::InputTag> clusterTag{Name("CrvCoincidenceClusterCollection"), Comment("CRV clusters; must be the persisted collection"), art::InputTag("SelectReco:CrvCoincidenceClusterFinder")};
        fhicl::Atom<art::InputTag> pulseTag{Name("CrvRecoPulseCollection"), Comment("all CRV reco pulses, for the layer probe; empty turns the probe off"), art::InputTag("CrvRecoPulses")};
        fhicl::Atom<art::InputTag> usedInFitTag{Name("UsedInFitAssns"), Comment("KalSeed-cluster associations of the clusters each fit used as constraints; empty: none"), art::InputTag()};
        fhicl::Atom<int> diag{Name("diagLevel"), Comment("diagnostic printout level"), 0};
        fhicl::TableFragment<TrkCrvMatchTool::Config> tool;
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit TrkCrvMatcher(const Parameters& conf);
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;

    private:
      art::ProductToken<KalSeedCollection> kalSeedToken_;
      art::ProductToken<CrvCoincidenceClusterCollection> clusterToken_;
      art::InputTag pulseTag_, usedInFitTag_;
      int diag_;
      TrkCrvMatchTool tool_;
      bool warnedNoPulses_ = false;
  };

  TrkCrvMatcher::TrkCrvMatcher(const Parameters& conf) :
    art::EDProducer(conf),
    kalSeedToken_(consumes<KalSeedCollection>(conf().kalSeedTag())),
    clusterToken_(consumes<CrvCoincidenceClusterCollection>(conf().clusterTag())),
    pulseTag_(conf().pulseTag()), usedInFitTag_(conf().usedInFitTag()), diag_(conf().diag()),
    tool_(conf().tool())
  {
    if(tool_.probeActive() && !pulseTag_.empty()) mayConsume<CrvRecoPulseCollection>(pulseTag_);
    if(!usedInFitTag_.empty()) consumes<TrkCrvFitClusterAssns>(usedInFitTag_);
    produces<TrkCrvMatchCollection>();
  }

  void TrkCrvMatcher::beginRun(art::Run&) { tool_.beginRun(); }

  void TrkCrvMatcher::produce(art::Event& event) {
    auto const& ksH = event.getValidHandle(kalSeedToken_);
    auto const& ccH = event.getValidHandle(clusterToken_);
    std::vector<art::Ptr<CrvCoincidenceCluster>> clusters;
    clusters.reserve(ccH->size());
    for(size_t ic = 0; ic < ccH->size(); ++ic) clusters.emplace_back(ccH, ic);

    CrvRecoPulseCollection const* pulses = nullptr;
    if(tool_.probeActive() && !pulseTag_.empty()) {
      auto pH = event.getHandle<CrvRecoPulseCollection>(pulseTag_);
      if(pH.isValid()) pulses = pH.product();
      else if(!warnedNoPulses_) {
        mf::LogWarning("TrkCrvMatcher") << "no CrvRecoPulseCollection " << pulseTag_ << ": the layer probe records predicted bars only";
        warnedNoPulses_ = true;
      }
    }
    tool_.setEvent(clusters, pulses);

    std::map<size_t, std::vector<art::Ptr<CrvCoincidenceCluster>>> used;   // KalSeed key -> clusters its fit used
    if(!usedInFitTag_.empty()) {
      auto const& aH = event.getValidHandle<TrkCrvFitClusterAssns>(usedInFitTag_);
      for(auto const& a : *aH) if(a.first.id() == ksH.id()) used[a.first.key()].push_back(a.second);
    }

    auto out = std::make_unique<TrkCrvMatchCollection>();
    out->reserve(ksH->size());
    static const std::vector<art::Ptr<CrvCoincidenceCluster>> none;
    for(size_t is = 0; is < ksH->size(); ++is) {
      auto iu = used.find(is);
      out->push_back(tool_.match(art::Ptr<KalSeed>(ksH, is), iu != used.end() ? iu->second : none));
    }
    tool_.resolve(*out);

    if(diag_ > 0) {
      for(auto const& m : *out) {
        mf::LogInfo log("TrkCrvMatcher");
        log << "event " << event.id() << " track " << m.track_.key() << " tRef " << m.tRef_ << " (" << m.refSurface_ << ")"
            << " xings " << m.xings_.size() << " cands " << m.cands_.size() << " bestGeom " << m.bestGeom_ << " bestTime " << m.bestTime_;
        if(diag_ > 1) {
          for(auto const& c : m.cands_)
            log << "\n  cluster " << c.cluster_.key() << " shield " << c.sector_ << " xing " << c.xing_ << " dtRef " << c.dtRef_
                << " dtXing " << c.dtXing_ << " du " << c.du_ << " dv " << c.dv_ << " chi2 " << c.chi2_ << " flags " << c.flags_;
          for(auto const& x : m.xings_)
            for(size_t L = 0; L < x.layers_.size(); ++L)
              log << "\n  probe shield " << x.sector_ << " layer " << L << " bar " << x.layers_[L].bar_ << " dBar " << int(x.layers_[L].dBar_)
                  << " uPred " << x.layers_[L].uPred_ << " dt " << x.layers_[L].dt_ << " pe " << x.layers_[L].pe_
                  << " inCluster " << x.layers_[L].inCluster_ << " uTdiv " << x.layers_[L].uTdiv_;
        }
      }
    }
    event.put(std::move(out));
  }
}

DEFINE_ART_MODULE(mu2e::TrkCrvMatcher)

// Fit the track dt_{hit} / dt_{track fit poca time} distribution
//
// Michael MacKenzie, 2026

#include <iostream>
#include <string>
#include <vector>

#include "cetlib_except/exception.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeedDtDt.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"

#include "TH1.h"

namespace mu2e {

  //================================================================
  class TrackDtDt : public art::EDProducer {

    struct Config
    {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Sequence<std::string> kalSeeds  { Name("kalSeeds") , Comment("KalSeed (Ptr) collection names") };
      fhicl::Atom<int>             diagLevel { Name("diagLevel"), Comment("Diag Level")                    ,0 };
    };

    std::vector<std::string> kalSeeds_;
    int diagLevel_;

  public:
    explicit TrackDtDt(const art::EDProducer::Table<Config>& config);
    virtual void produce(art::Event& event);
    KalSeedDtDt fitTrackDtDt(const KalSeed& seed);
  };

  //================================================================
  TrackDtDt::TrackDtDt(const art::EDProducer::Table<Config>& config)
    : art::EDProducer{config}
    , kalSeeds_(config().kalSeeds())
    , diagLevel_(config().diagLevel())
  {
    for(const auto& name : kalSeeds_) {
      produces<mu2e::KalSeedDtDtCollection>(name);
    }
  }

  //================================================================
  KalSeedDtDt TrackDtDt::fitTrackDtDt(const KalSeed& seed) {
    LsqSums2 fitter;
    for(const auto& hit : seed.hits()) {
      const double dt_hit = hit.time() - hit.TOTDriftTime();
      const double dt_trk = hit.particleToca();
      const double t_unc  = hit._udresid / hit._dvel;
      const double weight = 1./t_unc;
      fitter.addPoint(dt_trk, dt_hit, weight);
    }
    KalSeedDtDt result;
    result.slope_     = fitter.dydx();
    result.offset_    = fitter.y0();
    result.slopeUnc_  = fitter.dydxErr();
    result.dof_       = fitter.qn() - 2; // two parameters fitted
    result.chisq_     = fitter.chi2Dof()*result.dof_;
    return result;
  }

  //================================================================
  void TrackDtDt::produce(art::Event& event) {

    // Loop over all KalSeed collections
    for(const auto& name : kalSeeds_) {

      // Retrieve the KalSeed collection from the event, checking if it is a collection of KalSeed or KalSeedPtr
      art::Handle<KalSeedCollection> seedHandle;
      event.getByLabel(name, seedHandle);
      art::Handle<KalSeedPtrCollection> seedPtrHandle;
      const bool isSeedCollection = seedHandle.isValid();
      if(!isSeedCollection) {
        event.getByLabel(name, seedPtrHandle);
        if(!seedPtrHandle.isValid()) {
          throw cet::exception("RECO") << "TrackDtDt: No KalSeed or KalSeedPtr collection with label " << name << std::endl;
        }
      }

      // Create the output collection
      std::unique_ptr<KalSeedDtDtCollection> dtDtCol(new KalSeedDtDtCollection());

      // Loop over all seeds in the collection and fit the dt_{hit} / dt_{track fit poca time} distribution
      const auto nseeds = (isSeedCollection) ? seedHandle->size() : seedPtrHandle->size();
      if(diagLevel_ > 0) std::cout << "[TrackDtDt::" << __func__ << "] Processing " << nseeds << " seeds from collection " << name << std::endl;
      for(size_t iseed = 0; iseed < nseeds; ++iseed) {
        const auto& seed = (isSeedCollection) ? seedHandle->at(iseed) : *seedPtrHandle->at(iseed);
        if(diagLevel_ > 1) std::cout << "[TrackDtDt::" << __func__ << "] Processing seed " << iseed << " with " << seed.hits().size() << " hits" << std::endl;
        dtDtCol->emplace_back(fitTrackDtDt(seed));
      }

      // Put the results into the event
      event.put(std::move(dtDtCol), name);
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::TrackDtDt)

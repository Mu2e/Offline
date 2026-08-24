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
#include "TGraph.h"

namespace mu2e {

  //================================================================
  class TrackDtDt : public art::EDProducer {

    struct Config
    {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Sequence<std::string> kalSeeds    { Name("kalSeeds")    , Comment("KalSeed (Ptr) collection names") };
      fhicl::Atom<bool>            doHistograms{ Name("doHistograms"), Comment("Make histograms")               , false };
      fhicl::Atom<int>             diagLevel   { Name("diagLevel")   , Comment("Diag Level")                    ,0 };
    };

    std::vector<std::string> kalSeeds_;
    bool doHistograms_;
    int diagLevel_;

    // Histograms
    TH1* h_slope_ = nullptr; // fit results
    TH1* h_offset_ = nullptr;
    TH1* h_slopeUnc_ = nullptr;
    TH1* h_chisq_ = nullptr;
    TH1* h_dof_ = nullptr;
    TH1* h_chisqDof_ = nullptr;
    TH1* h_trk_chisq_ = nullptr; // track parameters
    TH1* h_trk_fitcon_ = nullptr;
    TGraph* g_t_hit_vs_z = nullptr;
    TGraph* g_t_poca_vs_z = nullptr;
    TGraph* g_t_hit_vs_t_poca = nullptr;

  public:
    explicit TrackDtDt(const art::EDProducer::Table<Config>& config);
    virtual void produce(art::Event& event);
    KalSeedDtDt fitTrackDtDt(const KalSeed& seed);

    void bookHistograms();
    void fillHistograms(const KalSeedDtDt& dtDt, const KalSeed& seed);
  };

  //================================================================
  TrackDtDt::TrackDtDt(const art::EDProducer::Table<Config>& config)
    : art::EDProducer{config}
    , kalSeeds_(config().kalSeeds())
    , doHistograms_(config().doHistograms())
    , diagLevel_(config().diagLevel())
  {
    for(const auto& name : kalSeeds_) {
      produces<mu2e::KalSeedDtDtCollection>(name);
    }
    if(doHistograms_) bookHistograms();
  }

  //================================================================
  void TrackDtDt::bookHistograms() {
    art::ServiceHandle<art::TFileService> tfs;
    h_slope_     = tfs->make<TH1D>("h_slope",     "Slope of dt_{hit} vs dt_{track fit poca time};Slope;Entries", 100, -2., 4.);
    h_offset_    = tfs->make<TH1D>("h_offset",    "Offset of dt_{hit} vs dt_{track fit poca time};Offset [ns];Entries", 100, -100., 100.);
    h_slopeUnc_  = tfs->make<TH1D>("h_slopeUnc",  "Uncertainty of slope;Slope Uncertainty;Entries", 100, 0., 0.1);
    h_chisq_     = tfs->make<TH1D>("h_chisq",     "Chi2 of fit;Chi2;Entries", 100, 0., 100.);
    h_dof_      = tfs->make<TH1D>("h_dof",      "Degrees of freedom of fit;DOF;Entries", 100, 0., 100.);
    h_chisqDof_   = tfs->make<TH1D>("h_chisqDof",   "Chi2/DOF of fit;Chi2/DOF;Entries", 100, 0., 5.);
    h_trk_chisq_ = tfs->make<TH1D>("h_trk_chisq", "Track chi2;Chi2;Entries", 100, 0., 100.);
    h_trk_fitcon_ = tfs->make<TH1D>("h_trk_fitcon", "Track fit convergence;Convergence;Entries", 100, 0., 1.);
    g_t_hit_vs_z = tfs->makeAndRegister<TGraph>("t_hit_vs_z", "t_{hit} vs. z;z [mm];t_{hit} [ns]");
    g_t_poca_vs_z = tfs->makeAndRegister<TGraph>("t_poca_vs_z", "t_{POCA} vs. z;z [mm];t_{POCA} [ns]");
    g_t_hit_vs_t_poca = tfs->makeAndRegister<TGraph>("t_hit_vs_t_poca", "t_{hit} vs. t_{POCA};t_{POCA} [ns];t_{hit} [ns]");
  }

  //================================================================
  void TrackDtDt::fillHistograms(const KalSeedDtDt& dtDt, const KalSeed& seed) {
    if(!doHistograms_) return;
    h_slope_->Fill(dtDt.slope_);
    h_offset_->Fill(dtDt.offset_);
    h_slopeUnc_->Fill(dtDt.slopeUnc_);
    h_chisq_->Fill(dtDt.chisq_);
    h_dof_->Fill(dtDt.dof_);
    h_chisqDof_->Fill((dtDt.dof_ > 0) ? dtDt.chisq_/dtDt.dof_ : 0.);
    h_trk_chisq_->Fill(seed.chisquared());
    h_trk_fitcon_->Fill(seed.fitConsistency());

    // Fill one example event
    if(g_t_hit_vs_z->GetN() == 0) {
      for(const auto& hit : seed.hits()) {
        const double t_trk = hit.particleToca(); // estimated time of the particle
        const double t_hit = t_trk + hit.fitDt(); // add dt to get the hit time
        const double z_trk = hit._upoca.z();
        g_t_hit_vs_z     ->AddPoint(z_trk, t_hit);
        g_t_poca_vs_z    ->AddPoint(z_trk, t_trk);
        g_t_hit_vs_t_poca->AddPoint(t_trk, t_hit);
      }
    }
  }

  //================================================================
  KalSeedDtDt TrackDtDt::fitTrackDtDt(const KalSeed& seed) {
    LsqSums2 fitter;
    for(const auto& hit : seed.hits()) {
      const double t_trk = hit.particleToca(); // estimated time of the particle
      const double t_hit = t_trk + hit.fitDt(); // add dt to get the hit time
      const double t_var  = hit.fitTocaVar(); // variance of the track poca time and the hit time
      const double weight = 1./t_var; // inverse variance weighting
      fitter.addPoint(t_trk, t_hit, weight);
    }
    KalSeedDtDt result;
    result.slope_     = fitter.dydx();
    result.offset_    = fitter.y0();
    result.slopeUnc_  = fitter.dydxErr();
    result.dof_       = fitter.qn() - 2; // two parameters fitted
    result.chisq_     = (result.dof_ > 0) ? fitter.chi2Dof()*result.dof_ : 0.f; // chi2Dof() not defined for dof <= 0

    if(diagLevel_ > 1) {
      std::cout << "[TrackDtDt::" << __func__ << "] Fit results for seed with " << seed.hits().size() << " hits:\n"
                << "  slope = " << result.slope_ << " +/- " << result.slopeUnc_ << "\n"
                << "  offset = " << result.offset_ << "\n"
                << "  chi2 = " << result.chisq_ << "\n"
                << "  dof = " << result.dof_ << "\n";
    }
    if(doHistograms_) fillHistograms(result, seed);
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

// Momentum distribution of our pions
//
// Andrei Gaponenko, 2018

#include <string>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"

#include "TH1.h"
#include "TH2.h"

namespace mu2e {

  class PionMomentumAnalyzer : public art::EDAnalyzer {
  public:

    struct ConfigStruct {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> inputs{Name("inputs"), Comment("Input SimParticleCollection")};

      fhicl::Atom<int>    nPBins{Name("nPBins"), Comment("number of bins for momentum histograms")};
      fhicl::Atom<double> pmin{Name("pmin"), Comment("Momentum histogram lower limit")};
      fhicl::Atom<double> pmax{Name("pmax"), Comment("Momentum histogram upper limit")};
    };

    typedef art::EDAnalyzer::Table<ConfigStruct> Conf;

    //----------------------------------------------------------------

    explicit PionMomentumAnalyzer(const Conf& config);
    explicit PionMomentumAnalyzer(const Conf& config, art::TFileDirectory tfdir);

    virtual void beginJob() override;
    virtual void analyze(const art::Event&) override;

  private:
    Conf conf_;

    TH1* h_p_all_;
    TH2* h_p_by_parent_;
    TH2* h_p_by_process_;

    const ParticleDataTable *particleTable_;

    static bool is_muon_daughter(const SimParticle& p);
  };

  //================================================================
  PionMomentumAnalyzer::PionMomentumAnalyzer(const Conf& config)
    : PionMomentumAnalyzer(config, *art::ServiceHandle<art::TFileService>())
  {}

  PionMomentumAnalyzer::PionMomentumAnalyzer(const Conf& c, art::TFileDirectory tf)
    : art::EDAnalyzer(c)
    , conf_{c}

    , h_p_all_{tf.make<TH1D>("p_pion", "Pion production momentum", c().nPBins(), c().pmin(), c().pmax())}
    , h_p_by_parent_{tf.make<TH2D>("p_pion_by_parent", "Pion production momentum vs parent PID",  1, 0., 0., c().nPBins(), c().pmin(), c().pmax())}
    , h_p_by_process_{tf.make<TH2D>("p_pion_by_process", "Pion production momentum vs production process",  1, 0., 0., c().nPBins(), c().pmin(), c().pmax())}
    , particleTable_{nullptr}
  {
    h_p_by_parent_->SetOption("colz");
    h_p_by_process_->SetOption("colz");
  }

  //================================================================
  void PionMomentumAnalyzer::beginJob() {
    GlobalConstantsHandle<ParticleDataTable> ph;
    particleTable_ = &(*ph);
  }

  //================================================================
  void PionMomentumAnalyzer::analyze(const art::Event& evt) {
    const auto sc = evt.getValidHandle<SimParticleCollection>(conf_().inputs());
    for(const auto& spe: *sc) {
      const SimParticle& p = spe.second;

      if(p.pdgId() == PDGCode::pi_minus) {

        const double momentum = p.startMomentum().vect().mag();
        h_p_all_->Fill(momentum);

        const SimParticle& parent{*p.parent()};

        const auto pref = particleTable_->particle(parent.pdgId()).ref();
        std::string parentName = pref.name();
        h_p_by_parent_->Fill(parentName.c_str(), momentum, 1.0);

        std::string codename = p.creationCode().name(); // need to bind for the c_str() call below
        h_p_by_process_->Fill(codename.c_str(), momentum, 1.0);
      }
    }
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::PionMomentumAnalyzer);

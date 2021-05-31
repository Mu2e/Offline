// Andrei Gaponenko, 2015

#include <string>
#include <vector>

#include "TDirectory.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace mu2e {

  //================================================================
  class CollectionSizeAnalyzer : public art::EDAnalyzer {
    TProfile *hStepPointSize_;
    TH1D *hStepPointSum_;

    TH2D *hStepPointDistLin_;
    TH2D *hStepPointDistLog_;

    TH2D *hSimParticleDistLin_;
    TH2D *hSimParticleDistLog_;

    bool useModuleLabel_;
    bool useInstanceName_;
    bool useProcessName_;

    template<class C>
    std::string getCollectionName(const art::Handle<C>& c) {
      std::ostringstream collName;
      if(useModuleLabel_) {
        collName<<c.provenance()->moduleLabel();
      }
      if(useInstanceName_) {
        if(!collName.str().empty()) {
          collName<<":";
        }
        collName<<c.provenance()->productInstanceName();
      }
      if(useProcessName_) {
        if(!collName.str().empty()) {
          collName<<":";
        }
        collName<<c.provenance()->processName();
      }

      return collName.str();
    }

    void doStepPoints(const art::Event& event);
    void doSimParticles(const art::Event& event);

  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<bool> useModuleLabel{Name("useModuleLabel"), Comment("Include input module label into collection name in the plots") };
      fhicl::Atom<bool> useInstanceName{Name("useInstanceName"), Comment("Include input instance name into collection name in the plots")};
      fhicl::Atom<bool> useProcessName{Name("useProcessName"), Comment("Include input process name into collection name in the plots")};
    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit CollectionSizeAnalyzer(const Parameters& conf);

    virtual void analyze(const art::Event& event);
  };

  //================================================================
  CollectionSizeAnalyzer::CollectionSizeAnalyzer(const Parameters& conf)
    : art::EDAnalyzer(conf)
    , hStepPointSize_(0)
    , useModuleLabel_(conf().useModuleLabel())
    , useInstanceName_(conf().useInstanceName())
    , useProcessName_(conf().useProcessName())
  {
    art::ServiceHandle<art::TFileService> tfs;
    hStepPointSize_ = tfs->make<TProfile>("avgStepPointsSize", "Average step point collection size", 1, 0., 0.);
    hStepPointSum_ = tfs->make<TH1D>("stepPointsSum", "Sum of step point collection entries", 1, 0., 0.);

    // We use automatic binning on the X axis to accommodate an
    // apriory unknown set of collections.  ROOT's kCanRebin is a
    // property of a histogram, not an axis, therefore the other axis
    // will also be auto-rebinned, whether we want it or not.

    hStepPointDistLin_ = tfs->make<TH2D>("stepPointsSizeDistributionLiny", "Step point multiplicity vs collection name",
                                         1, 0., 0., 1000, -0.5, 999.5);
    hStepPointDistLin_->SetOption("colz");

    hStepPointDistLog_ = tfs->make<TH2D>("stepPointsSizeDistributionLogy", "log10(step point multiplicity) vs collection name",
                                      1, 0., 0., 100, 0., 10.);
    hStepPointDistLog_->SetOption("colz");

    hSimParticleDistLin_ = tfs->make<TH2D>("simParticlesSizeDistributionLin", "SimParticle multiplicity vs collection name",
                                        1, 0., 0., 1000, -0.5, 999.5);
    hSimParticleDistLin_->SetOption("colz");

    hSimParticleDistLog_ = tfs->make<TH2D>("simParticlesSizeDistributionLog", "log10(SimParticle multiplicity) vs collection name",
                                        1, 0., 0., 100, 0., 10.);
    hSimParticleDistLog_->SetOption("colz");
  }

  //================================================================
  void CollectionSizeAnalyzer::doStepPoints(const art::Event& event) {
    std::vector<art::Handle<StepPointMCCollection> > stepHandles = event.getMany<StepPointMCCollection>();
    for(const auto& c: stepHandles) {
      const std::string cn = getCollectionName(c);
      hStepPointSize_->Fill(cn.c_str(), double(c->size()));
      hStepPointSum_->Fill(cn.c_str(), double(c->size()));
      hStepPointDistLin_->Fill(cn.c_str(), double(c->size()), 1.);
      hStepPointDistLog_->Fill(cn.c_str(), log10(std::max(1., double(c->size()))), 1.);
    }
  }

  //================================================================
  void CollectionSizeAnalyzer::doSimParticles(const art::Event& event) {
    std::vector<art::Handle<SimParticleCollection> > stepHandles = event.getMany<SimParticleCollection>();
    for(const auto& c: stepHandles) {
      const std::string cn = getCollectionName(c);
      hSimParticleDistLin_->Fill(cn.c_str(), double(c->size()), 1.);
      hSimParticleDistLog_->Fill(cn.c_str(), log10(std::max(1., double(c->size()))), 1.);
    }
  }

  //================================================================
  void CollectionSizeAnalyzer::analyze(const art::Event& event) {
    doStepPoints(event);
    doSimParticles(event);
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CollectionSizeAnalyzer);

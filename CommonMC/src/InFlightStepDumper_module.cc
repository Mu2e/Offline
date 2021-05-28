// Write out info from a StepPointMC collection for re-sampling in the following stage.
//
// Andrei Gaponenko, 2015

#include <string>
#include <algorithm>

#include "cetlib_except/exception.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "art_root_io/TFileService.h"

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "GeneralUtilities/inc/RSNTIO.hh"

#include "TTree.h"

#include <algorithm>
#include <iterator>

namespace mu2e {

  //================================================================
  class InFlightStepDumper : public art::EDAnalyzer {
  public:
    explicit InFlightStepDumper(fhicl::ParameterSet const& pset);
    void beginJob() override;
    void analyze(const art::Event& evt) override;
  private:
    art::InputTag input_;
    TTree *nt_;
    int pie_; // particle number in the current event
    IO::InFlightParticleD data_;
  };

  //================================================================
  InFlightStepDumper::InFlightStepDumper(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , input_(pset.get<std::string>("inputCollection"))
    , nt_()
    , pie_()
  {}

  //================================================================
  void InFlightStepDumper::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    nt_ = tfs->make<TTree>( "particles", "In-flight particles ntuple");
    nt_->Branch("particles", &data_, IO::InFlightParticleD::branchDescription());
    nt_->Branch("pie", &pie_);
  }

  //================================================================
  void InFlightStepDumper::analyze(const art::Event& event) {
    auto ih = event.getValidHandle<StepPointMCCollection>(input_);
    pie_ = 0;
    for(const auto& hit : *ih) {
      data_.x = hit.position().x();
      data_.y = hit.position().y();
      data_.z = hit.position().z();

      data_.time = hit.time();

      data_.px = hit.momentum().x();
      data_.py = hit.momentum().y();
      data_.pz = hit.momentum().z();

      data_.pdgId = hit.simParticle()->pdgId();

      nt_->Fill();
      ++pie_;
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::InFlightStepDumper);

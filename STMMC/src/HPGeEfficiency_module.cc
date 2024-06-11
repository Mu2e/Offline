#include <iostream>
#include <string>
#include <cmath>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"

#include "canvas/Utilities/InputTag.h"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"

#include "TNtuple.h"

#include "art_root_io/TFileService.h"
using namespace std;

namespace mu2e{
  class HPGeEfficiency : public art::EDAnalyzer
  {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config
    {
      fhicl::Atom<bool> v{Name("v"), Comment("verbose")};
    };
    using Parameters=art::EDAnalyzer::Table<Config>;

    explicit HPGeEfficiency(const Parameters& pset);
    void analyze(art::Event const& event) override;
    void beginJob() override;

  private:
    bool verbose = false;
    TNtuple* _nt = nullptr;
  };
  // ===================================================
  HPGeEfficiency::HPGeEfficiency(const Parameters& config):
    art::EDAnalyzer{config},
    verbose(config().v())
    {};
  // ===================================================
  void HPGeEfficiency::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    _nt = tfs->make<TNtuple>("nt", "Energy Deposit", "E");
    return;
  };
  // ===================================================
  void HPGeEfficiency::analyze(art::Event const& event)
  {
    art::Handle<StepPointMCCollection> _inputStepPointMCs;
    event.getByLabel("g4run:STMDet", _inputStepPointMCs);

    double E = 0.0;
    for (const StepPointMC &step : *_inputStepPointMCs)
    {
      E += step.ionizingEdep();
    };
    _nt->Fill(E);
    return;
  };
}

DEFINE_ART_MODULE(mu2e::HPGeEfficiency)

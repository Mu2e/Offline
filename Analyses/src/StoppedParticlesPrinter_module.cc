// Print out particles from a SimParticlePtrCollections into a text file.
//
// Andrei Gaponenko, 2013

#include <string>
#include <fstream>
#include <iomanip>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"

namespace mu2e {

  //================================================================
  class StoppedParticlesPrinter : public art::EDAnalyzer {
  public:
    explicit StoppedParticlesPrinter(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt) override;
  private:
    art::InputTag input_;
    std::ofstream outFile_;
  };

  //================================================================
  StoppedParticlesPrinter::StoppedParticlesPrinter(const fhicl::ParameterSet& pset):
    art::EDAnalyzer(pset),
    input_(pset.get<std::string>("inputCollection"))
  {
    const std::string outFileName(pset.get<std::string>("outFileName"));
    outFile_.open(outFileName.c_str());
    if(pset.get<bool>("writeBeginDataMarker", true)) {
      outFile_<< "begin data" << std::endl;
    }
  }

  //================================================================
  void StoppedParticlesPrinter::analyze(const art::Event& event) {
    auto ih = event.getValidHandle<SimParticlePtrCollection>(input_);
    for(const auto& p : *ih) {
      outFile_<< std::setprecision(8)
              << p->endPosition().x() << " "
              << p->endPosition().y() << " "
              << p->endPosition().z() << " "
              << p->endGlobalTime()
              << std::endl;
    }
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedParticlesPrinter);

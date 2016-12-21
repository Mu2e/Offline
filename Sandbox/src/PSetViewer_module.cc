//
// Look at the content of the parameter set registry.
// and the ProcessHistoryRegistry.
//
// Contact person Rob Kutschke
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSetRegistry.h"
#include "art/Persistency/Provenance/ProcessHistoryRegistry.h"

#include <iostream>

namespace mu2e {

  class PSetViewer : public art::EDAnalyzer {
  public:

    explicit PSetViewer(fhicl::ParameterSet const& pset);

    virtual void analyze(const art::Event& e);

  };

  PSetViewer::PSetViewer(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset){
  }

  void PSetViewer::analyze( const art::Event& event ) {
    std::cout << "Event: " << event.id() << std::endl;
    std::cout << "Parameter set registry size: " << fhicl::ParameterSetRegistry::get().size() << std::endl;
    int n(-1);
    for ( auto const& i : fhicl::ParameterSetRegistry::get() ){
      std::cout << "\nParameterSet: " << ++n << ": ID: " << i.first << std::endl;
      std::cout << i.second.to_indented_string() << std::endl;
    }

    n=-1;
    std::cout << "Process history registry size: " << art::ProcessHistoryRegistry::get().size() << std::endl;
    for ( auto const& i : art::ProcessHistoryRegistry::get() ){
      std::cout << "\nProcessHistory: " << ++n << ": ID: " << i.first << std::endl;
      std::cout << i.second << std::endl;
    }

  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::PSetViewer);

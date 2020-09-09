//
//  Look at all provenances stored in the registry
//
//
//  Original author Rob Kutschke
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include <iostream>

using namespace std;

namespace mu2e {

  class ViewAllProvenances : public art::EDAnalyzer {

  public:
    explicit ViewAllProvenances(fhicl::ParameterSet const& pset) : art::EDAnalyzer(pset) {}

    void analyze(const art::Event& event) override;

  private:

  };

  void ViewAllProvenances::analyze(const art::Event& event){

    cout << "Number of enteries in the parameter set registry: " <<  fhicl::ParameterSetRegistry::get().size() << endl;

    // The type of tmp is some sort of key value pair.
    for ( auto const& tmp : fhicl::ParameterSetRegistry::get() ){

      fhicl::ParameterSet const& pset = tmp.second;

      cout << "\nParameter set content: " << endl;
      cout << "--------------------------------------" << endl;
      cout << pset.to_indented_string();
      cout << "--------------------------------------" << endl;
    }


  } // end analyze

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ViewAllProvenances);

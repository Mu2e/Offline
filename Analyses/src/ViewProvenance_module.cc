//
//  A module to look at the provenance of a product.
//
//
//  Original author Rob Kutschke
//

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

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

  class ViewProvenance : public art::EDAnalyzer {

  public:
    explicit ViewProvenance(fhicl::ParameterSet const& pset);

    void analyze(const art::Event& event) override;

  private:

    void printProvenance( art::Provenance const& );

  };

  ViewProvenance::ViewProvenance(fhicl::ParameterSet const& pset)
    : art::EDAnalyzer(pset){
  }

  void ViewProvenance::analyze(const art::Event& event){

    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel("generate", gensHandle );
    printProvenance( *gensHandle.provenance( ) );

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel("g4run", simsHandle );
    printProvenance( *simsHandle.provenance( ) );

  } // end analyze

  void ViewProvenance::printProvenance( art::Provenance const& prov ) {

    cout << "\nProvenance information for product name: "  << prov.branchName()  << endl;
    cout << "Creator module label: " << prov.moduleLabel() << endl;

    std::set<fhicl::ParameterSetID> const& ids = prov.productDescription().psetIDs();
    if ( ids.size() == 1 ){
      cout << "Parameter set used to configure this module: " << endl;
    }else{
      cout << "Parameter sets used to configure this module: " << endl;
    }

    for ( auto const& id : ids ){
      if ( ids.size() > 1 ){
        cout << "\nParameter Set content: " << endl;
      }
      cout << "--------------------------------------" << endl;
      fhicl::ParameterSet const& pset = fhicl::ParameterSetRegistry::get(id);
      cout << pset.to_indented_string();
      cout << "--------------------------------------" << endl;
    }

  }

} // end namespace mu2e

using mu2e::ViewProvenance;
DEFINE_ART_MODULE(ViewProvenance);

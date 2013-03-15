//
//  A module to look at the provenance of a product.
//
//  $Id: ViewProvenance_module.cc,v 1.2 2013/03/15 15:52:03 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2013/03/15 15:52:03 $
//
//  Original author Rob Kutschke
//

#include "MCDataProducts/inc/GenParticleCollection.hh"

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
    explicit ViewProvenance(fhicl::ParameterSet const& pset){}

    void analyze(const art::Event& event);

  private:

  };

  void ViewProvenance::analyze(const art::Event& event){

    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel("generate", gensHandle );

    art::Provenance const * prov = gensHandle.provenance( );

    if ( prov->isPresent() ) {
      cout << "Product name: "
           << prov->branchName()
           << endl;
      cout << "Creator module label: " << prov->moduleLabel() << endl;
      cout << "Parameter set used to configure this module: " << endl;
      // Fixme: this is broken in art /v1_03_08
      //fhicl::ParameterSet const& pset =
      //fhicl::ParameterSetRegistry::get( prov->productDescription().parameterSetID());
      //cout << pset.to_indented_string()
      //   << endl;
    } else{
      cout << "This product has no stored provenance ... " << endl;
    }

  } // end analyze

} // end namespace mu2e

using mu2e::ViewProvenance;
DEFINE_ART_MODULE(ViewProvenance);

//
// Find a product, modify it and rewrite it.
//
// Original author Rob Kutschke
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

using namespace std;

namespace mu2e {

  class ModifyTestProduct : public art::EDProducer {
  public:

    explicit ModifyTestProduct(fhicl::ParameterSet const& pset);

    void produce( art::Event& e) override;

  private:
    art::InputTag tag_;

  };

  ModifyTestProduct::ModifyTestProduct(fhicl::ParameterSet const& pset ):
    art::EDProducer{pset},
    tag_(pset.get<std::string>("productTag")){
    produces<int>( tag_.instance() );
  }

  void ModifyTestProduct::produce(art::Event& event) {

    auto testp = event.getValidHandle<int>(tag_);

    unique_ptr<int> prod = std::make_unique<int>(*testp+29);
    event.put(std::move(prod));

  } // end ModifyTestProduct::analyze

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ModifyTestProduct)

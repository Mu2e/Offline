//
// Create a product in an output file.
//
// Original author Rob Kutschke
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

using namespace std;

namespace mu2e {

  class MakeTestProduct : public art::EDProducer {
  public:

    explicit MakeTestProduct(fhicl::ParameterSet const& pset);

    void produce( art::Event& e) override;

  private:

  };

  MakeTestProduct::MakeTestProduct(fhicl::ParameterSet const& pset): 
    art::EDProducer{pset}
  {
    produces<int>();
  }

  void MakeTestProduct::produce(art::Event& event) {

    unique_ptr<int> prod = std::make_unique<int>(event.id().event());
    event.put(std::move(prod));

  } // end MakeTestProduct::analyze

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::MakeTestProduct)

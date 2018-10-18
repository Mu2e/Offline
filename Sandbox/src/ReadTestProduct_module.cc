//
// Create a product in an output file.
//
// Original author Rob Kutschke
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

using namespace std;

namespace mu2e {

  class ReadTestProduct : public art::EDAnalyzer {
  public:

    explicit ReadTestProduct(fhicl::ParameterSet const& pset);

    void analyze( art::Event const& e) override;

  private:
    art::InputTag tag_;

  };

  ReadTestProduct::ReadTestProduct(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    tag_(pset.get<std::string>("productTag"))
  {}

  void ReadTestProduct::analyze(art::Event const& event) {

    auto testp = event.getValidHandle<int>(tag_);
    cout << "Event: " << event.id() << "  Value: " << *testp << endl;


  } // end ReadTestProduct::analyze

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ReadTestProduct)

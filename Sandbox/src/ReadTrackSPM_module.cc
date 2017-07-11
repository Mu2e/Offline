//
// Create a product in an output file.
//
// Original author Rob Kutschke
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"

using namespace std;

namespace mu2e {

  class ReadTrackSPM : public art::EDAnalyzer {
  public:

    explicit ReadTrackSPM(fhicl::ParameterSet const& pset);

    void analyze( art::Event const& e) override;

  private:
    art::InputTag tag_;

  };

  ReadTrackSPM::ReadTrackSPM(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    tag_(pset.get<std::string>("productTag"))
  {}

  void ReadTrackSPM::analyze(art::Event const& event) {

    auto testp = event.getValidHandle<std::vector<mu2e::StepPointMC> >(tag_);

    cout << "Event: " << event.id() << " contains the following:  " << endl;

    int count = 0;
    for ( auto theTSPM : *testp ) {
      cout << ++count << theTSPM << endl;
    }


  } // end ReadTrackSPM::analyze

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ReadTrackSPM)

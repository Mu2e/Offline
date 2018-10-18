//
// Read and report a list of products from an input file.
//
// Original author Rob Kutschke
// Modified by David (Lou) Brown
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"

#include <vector>
using namespace std;

namespace mu2e {

  class ReadTrackSPM : public art::EDAnalyzer {
  public:

    explicit ReadTrackSPM(fhicl::ParameterSet const& pset);

    void analyze( art::Event const& e) override;

  private:
    vector<string> tags_;

  };

  ReadTrackSPM::ReadTrackSPM(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    tags_(pset.get<vector<string> >("productTags",vector<string>()))
  {}

  void ReadTrackSPM::analyze(art::Event const& event) {

    cout << "Event: " << event.id() << " contains the following:  " << endl;

    for ( auto aTag : tags_ ) { 
      art::InputTag theTag(aTag);
      auto testp = event.getValidHandle<vector<mu2e::StepPointMC> >(theTag);
      if ( 0 != testp ) {  // This is probably redundant.  Think about
	cout << "  For Tag " << aTag << ": " << endl;
	int count = 0;
	for ( auto theTSPM : *testp ) {
	  cout << ++count << theTSPM << endl;
	} // end for loop over step points
      } // end if 0 != testp
    } // end of loop over tags
  } // end ReadTrackSPM::analyze

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ReadTrackSPM)

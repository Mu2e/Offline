//
// Read the data product read in by Source00_module and print out some information.
//
// Original author Rob Kutschke.
//

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// C++ includes.
#include <iostream>

using namespace std;

namespace mu2e {

  class ReadSourceTest : public art::EDAnalyzer {
  public:
    explicit ReadSourceTest(fhicl::ParameterSet const& );

    void analyze( art::Event const& e);

  private:

  };

  ReadSourceTest::ReadSourceTest(fhicl::ParameterSet const& pset )
    : art::EDAnalyzer(pset){
  }

  void
  ReadSourceTest::analyze(art::Event const& event) {

    //art::Handle<int> handle;
    //event.getByLabel("Source00",handle);

    auto productHandle = event.getValidHandle<int>("Source00");

    cout << "Read event: " << event.id() << " " << *productHandle << endl;

  } // end of ::analyze.

}

using mu2e::ReadSourceTest;
DEFINE_ART_MODULE(ReadSourceTest)

//
//  The HelloWorld plugin; the first example of a module.
//
//
//  Add a TODO
//
//  Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"

using namespace std;

namespace mu2e {

  class HelloWorld : public art::EDAnalyzer {

  public:
    explicit HelloWorld(fhicl::ParameterSet const& pset);

    void analyze(const art::Event& event);

  private:

  };

  HelloWorld::HelloWorld(fhicl::ParameterSet const& pset)
    : art::EDAnalyzer(pset){
  }

  void HelloWorld::analyze(const art::Event& event){
    cerr << "Hello, world.  From analyze: "
         << event.id()
         << endl;
  }

} // end namespace mu2e

using mu2e::HelloWorld;
DEFINE_ART_MODULE(HelloWorld)

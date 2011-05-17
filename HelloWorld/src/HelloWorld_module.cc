//
//  The HelloWorld plugin; the first example of a module.
//
//  $Id: HelloWorld_module.cc,v 1.2 2011/05/17 22:22:46 wb Exp $
//  $Author: wb $
//  $Date: 2011/05/17 22:22:46 $
//   
//  Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

using namespace std;

namespace mu2e {

  class HelloWorld : public art::EDAnalyzer {

  public:
    explicit HelloWorld(fhicl::ParameterSet const& pset){}

    void analyze(const art::Event& event);

  private:

  };

  void HelloWorld::analyze(const art::Event& event){
    cerr << "Hello, world.  From analyze: "
         << event.id()
         << endl;
  }

} // end namespace mu2e

using mu2e::HelloWorld;
DEFINE_ART_MODULE(HelloWorld);

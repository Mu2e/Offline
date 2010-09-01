//
//  The HelloWorld plugin; the first example of a module.
//
//  $Id: HelloWorld_plugin.cc,v 1.3 2010/09/01 18:55:51 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2010/09/01 18:55:51 $
//   
//  Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

namespace mu2e {

  class HelloWorld : public edm::EDAnalyzer {

  public:
    explicit HelloWorld(edm::ParameterSet const& pset){}

    void analyze(const edm::Event& event, edm::EventSetup const&);

  private:

  };

  void HelloWorld::analyze(const edm::Event& event, edm::EventSetup const&){
    cerr << "Hello, world.  From analyze: "
         << event.id()
         << endl;
  }

} // end namespace mu2e

using mu2e::HelloWorld;
DEFINE_FWK_MODULE(HelloWorld);

//
//  The hello world plugin, for a very basic introduction.
//
//  $Id: HelloWorld_plugin.cc,v 1.1 2010/04/16 15:13:00 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2010/04/16 15:13:00 $
//   
//  Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

namespace mu2e {

  class HelloWorld : public edm::EDAnalyzer {

  public:
    explicit HelloWorld(edm::ParameterSet const& pset){
      cerr << "Hello, world.  From constructor. "
           << "My magic number is: "
           << pset.pset.getUntrackedParameter<int>("magicNumber",-1)
           << endl;
    }
    virtual ~HelloWorld() { }

    virtual void beginJob(edm::EventSetup const&);

    virtual void beginRun(edm::Run const &run, 
                          edm::EventSetup const& eSetup );
 
    void analyze(const edm::Event& e, edm::EventSetup const&);

  private:

  };

  void HelloWorld::beginJob(edm::EventSetup const& ){
    cerr << "Hello, world.  From beginJob." << endl;
  }

  void HelloWorld::beginRun(edm::Run const& run,
                            edm::EventSetup const& eSetup ){

    cerr << "Hello, world.  From beginRun: "
         << run.id().run()
         << endl;
  }

  void HelloWorld::analyze(const edm::Event& event, edm::EventSetup const&) {
    cerr << "Hello, world.  From analyze: "
         << event.id()
         << endl;
  }

} // end namespace mu2e

using mu2e::HelloWorld;
DEFINE_FWK_MODULE(HelloWorld);

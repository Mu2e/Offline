//
//  A second hello world plugin, with a little more detail.
//
//  $Id: HelloWorld2_plugin.cc,v 1.1 2010/09/01 18:55:36 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2010/09/01 18:55:36 $
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

  class HelloWorld2 : public edm::EDAnalyzer {

  public:
    explicit HelloWorld2(edm::ParameterSet const& pset):
      _magicNumber(pset.getUntrackedParameter<int>("magicNumber",-1)){

      cerr << "Hello, world.  From constructor. "
           << "My magic number is: "
           << _magicNumber
           << endl;
    }

    ~HelloWorld2() { }

    void beginJob(edm::EventSetup const&);

    void beginRun(edm::Run const &run, 
                          edm::EventSetup const& eSetup );
 
    void analyze(const edm::Event& event, edm::EventSetup const&);

    void endJob();

  private:

    int _magicNumber;

  };

  void HelloWorld2::beginJob(edm::EventSetup const& ){
    cerr << "Hello, world.  From beginJob. "
         << "  Magic number: " 
         << _magicNumber
         << endl;
  }

  void HelloWorld2::beginRun(edm::Run const& run,
                              edm::EventSetup const& eSetup ){

    cerr << "Hello, world.  From beginRun: "
         << run.id().run()
         << "  Magic number: " 
         << _magicNumber
         << endl;
  }

  void HelloWorld2::analyze(const edm::Event& event, edm::EventSetup const&) {
    cerr << "Hello, world.  From analyze: "
         << event.id()
         << "  Magic number: " 
         << _magicNumber
         << endl;
  }

  void HelloWorld2::endJob(){
    cerr << "Hello, world.  From endJob. "
         << "  Magic number: " 
         << _magicNumber
         << endl;
  }


} // end namespace mu2e

using mu2e::HelloWorld2;
DEFINE_FWK_MODULE(HelloWorld2);

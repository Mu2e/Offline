//
//  A second hello world plugin, with a little more detail.
//
//  $Id: HelloWorld2_module.cc,v 1.1 2011/05/17 16:30:14 greenc Exp $
//  $Author: greenc $
//  $Date: 2011/05/17 16:30:14 $
//   
//  Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"

using namespace std;

namespace mu2e {

  class HelloWorld2 : public art::EDAnalyzer {

  public:
    explicit HelloWorld2(fhicl::ParameterSet const& pset):
      _magicNumber(pset.get<int>("magicNumber",-1)){

      cerr << "Hello, world.  From constructor. "
           << "My magic number is: "
           << _magicNumber
           << endl;
    }

    ~HelloWorld2() { }

    void beginJob(art::EventSetup const&);

    void beginRun(art::Run const &run, 
                          art::EventSetup const& eSetup );
 
    void analyze(const art::Event& event, art::EventSetup const&);

    void endJob();

  private:

    int _magicNumber;

  };

  void HelloWorld2::beginJob(art::EventSetup const& ){
    cerr << "Hello, world.  From beginJob. "
         << "  Magic number: " 
         << _magicNumber
         << endl;
  }

  void HelloWorld2::beginRun(art::Run const& run,
                              art::EventSetup const& eSetup ){

    cerr << "Hello, world.  From beginRun: "
         << run.id().run()
         << "  Magic number: " 
         << _magicNumber
         << endl;
  }

  void HelloWorld2::analyze(const art::Event& event, art::EventSetup const&) {
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
DEFINE_ART_MODULE(HelloWorld2);

//
//  A second hello world plugin, with a little more detail.
//
//
//  Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

using namespace std;

namespace mu2e {

  class HelloWorld2 : public art::EDAnalyzer {

  public:
    explicit HelloWorld2(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset),
      _magicNumber(pset.get<int>("magicNumber",-1)){

      cerr << "Hello, world.  From constructor. "
           << "My magic number is: "
           << _magicNumber
           << endl;
    }

    ~HelloWorld2() { }

    void beginJob();

    void beginRun(art::Run const &run);

    void analyze(const art::Event& event);

    void endJob();

  private:

    int _magicNumber;

  };

  void HelloWorld2::beginJob(){
    cerr << "Hello, world.  From beginJob. "
         << "  Magic number: "
         << _magicNumber
         << endl;
  }

  void HelloWorld2::beginRun(art::Run const& run){

    cerr << "Hello, world.  From beginRun: "
         << run.id().run()
         << "  Magic number: "
         << _magicNumber
         << endl;
  }

  void HelloWorld2::analyze(const art::Event& event) {
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

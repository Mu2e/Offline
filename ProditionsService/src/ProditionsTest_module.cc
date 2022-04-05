///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

// C++ includes
#include <iostream>
#include <stdexcept>
#include <string>

// Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"


namespace mu2e {

  class ProditionsTest : public art::EDAnalyzer {

  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<int> verbose{Name("verbose"),
          Comment("verbose flag, 0 to 10"),1};

    };

    // this line is required by art to allow the command line help print
    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit ProditionsTest(const Parameters& conf):
      art::EDAnalyzer(conf),_conf(conf()) {}

    ~ProditionsTest() {}
    void analyze(const art::Event& event) override;

  private:

    Config _conf;

    ProditionsHandle<StrawResponse> _srh;
  };

//-----------------------------------------------------------------------------
  void ProditionsTest::analyze(const art::Event& event) {

    std::cout << "ProditionsTest::analyze  " << event.id() << std::endl;

    _srh.get(event.id());

  }
};

DEFINE_ART_MODULE(mu2e::ProditionsTest);

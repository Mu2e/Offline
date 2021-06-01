///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

// C++ includes
#include <iostream>
#include <stdexcept>
#include <string>

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
//#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "fhiclcpp/ParameterSet.h"
#include "DbService/inc/DbHandle.hh"
#include "DbTables/inc/TstCalib1.hh"

namespace mu2e {

  class DbServiceTest : public art::EDAnalyzer {

  public:

    explicit DbServiceTest(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset) {}

    ~DbServiceTest() {}
    void beginJob() override;
    void beginSubRun(const art::SubRun& subrun) override;
    void analyze(const art::Event& event) override;

  private:
    mu2e::DbHandle<mu2e::TstCalib1> _testCalib1;
  };

//-----------------------------------------------------------------------------
// Get access to the TFile service and book histograms
//-----------------------------------------------------------------------------
  void DbServiceTest::beginJob(){
  }

  void DbServiceTest::beginSubRun(const art::SubRun& subrun){
  }

//-----------------------------------------------------------------------------
  void DbServiceTest::analyze(const art::Event& event) {

    std::cout << "DbServiceTest::analyze" << std::endl;

    auto const& myTable = _testCalib1.get(event.id());
    std::cout << myTable.csv();
  };
};

DEFINE_ART_MODULE(mu2e::DbServiceTest);

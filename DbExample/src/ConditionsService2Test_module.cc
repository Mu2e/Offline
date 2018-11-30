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

#include "fhiclcpp/ParameterSet.h"
#include "DbExample/inc/ConditionsHandle2.hh"
#include "DbExample/inc/DetData1.hh"
#include "DbExample/inc/DetData2.hh"
#include "DbExample/inc/DetData3.hh"

namespace mu2e {

  class ConditionsService2Test : public art::EDAnalyzer {

  public:

    explicit ConditionsService2Test(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset) {}

    ~ConditionsService2Test() {}
    //void beginJob() override;
    //void beginSubRun(const art::SubRun& subrun) override;
    void analyze(const art::Event& event) override;

  private:
    ConditionsHandle2<DetData1> _detData1_h;
    ConditionsHandle2<DetData2> _detData2_h;
    ConditionsHandle2<DetData3> _detData3_h;
  };

//-----------------------------------------------------------------------------
// Get access to the TFile service and book histograms
//-----------------------------------------------------------------------------
//  void ConditionsService2Test::beginJob(){
//  }
//
//  void ConditionsService2Test::beginSubRun(const art::SubRun& subrun){
//  }

//-----------------------------------------------------------------------------
  void ConditionsService2Test::analyze(const art::Event& event) {

    std::cout << "ConditionsService2Test::analyze" << std::endl;

    auto const& d1 = _detData1_h.get(event.id());
    auto const& d2 = _detData2_h.get(event.id());
    auto const& d3 = _detData3_h.get(event.id());
    
    std::cout 
      << "r/s/e"  << std::setw(5) << event.run()
      << "/" << std::setw(5) << event.subRun()
      << "/" << std::setw(5) << event.event()
      << "   data1: " << std::setw(6) << std::setprecision(2) 
      << d1.getData()
      << "   data2: " ;
    for(auto x : d2.getData()) {
      std::cout 
	<< std::setw(6) << std::setprecision(2) 
	<< x;
    }
    std::cout  << "   data3: " ;
    for(auto x : d3.getData()) {
      std::cout 
	<< std::setw(6) << std::setprecision(2) 
	<< x;
    }
    std::cout << std::endl;

  };
};

DEFINE_ART_MODULE(mu2e::ConditionsService2Test);

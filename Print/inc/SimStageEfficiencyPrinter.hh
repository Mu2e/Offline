//
//  Utility class to print SimStageEfficiency
// 
#ifndef Print_inc_SimStageEfficiencyPrinter_hh
#define Print_inc_SimStageEfficiencyPrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/SimStageEfficiency.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class SimStageEfficiencyPrinter : public ProductPrinter {
  public:

    struct Config : public ProductPrinter::Config {
      fhicl::Atom<bool> forDb{ fhicl::Name("forDb"), 
	  fhicl::Comment("true/false whether to print SimStageEfficiencies in format for Db"),false};
    };

    SimStageEfficiencyPrinter(): _forDb(false) { }
    SimStageEfficiencyPrinter(const Config& conf):ProductPrinter(conf) { 
      _forDb = conf.forDb();
    }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override {};
    void PrintRun(art::Run const& run,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<SimStageEfficiency>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<SimStageEfficiency>& handle, 
	       std::ostream& os = std::cout);
    void Print(const SimStageEfficiency& obj, 
	       std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);

  private:
    bool _forDb;
  };

}
#endif

//
//  Utility class to print StrawDigiMC
// 
#ifndef Print_inc_StrawDigiMCPrinter_hh
#define Print_inc_StrawDigiMCPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class StrawDigiMCPrinter : public ProductPrinter {
  public:

    StrawDigiMCPrinter() { }
    StrawDigiMCPrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<StrawDigiMCCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<StrawDigiMCCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const StrawDigiMCCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<StrawDigiMC>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::StrawDigiMC& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif

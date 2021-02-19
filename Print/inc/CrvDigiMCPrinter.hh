//
//  Utility class to print CrvDigiMC
// 
#ifndef Print_inc_CrvDigiMCPrinter_hh
#define Print_inc_CrvDigiMCPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/CrvDigiMC.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class CrvDigiMCPrinter : public ProductPrinter {
  public:

    CrvDigiMCPrinter() { }
    CrvDigiMCPrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<CrvDigiMCCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<CrvDigiMCCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const CrvDigiMCCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<CrvDigiMC>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::CrvDigiMC& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif

//
//  Utility class to print CrvDigi
// 
#ifndef Print_inc_CrvDigiPrinter_hh
#define Print_inc_CrvDigiPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class CrvDigiPrinter : public ProductPrinter {
  public:

    CrvDigiPrinter() { }
    CrvDigiPrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<CrvDigiCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<CrvDigiCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const CrvDigiCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<CrvDigi>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::CrvDigi& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif

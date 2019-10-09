//
//  Utility class to print CaloDigi
// 
#ifndef Print_inc_CaloDigiPrinter_hh
#define Print_inc_CaloDigiPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class CaloDigiPrinter : public ProductPrinter {
  public:

    CaloDigiPrinter() { }
    CaloDigiPrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<CaloDigiCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<CaloDigiCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const CaloDigiCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<CaloDigi>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::CaloDigi& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif

//
//  Utility class to print StatusG4
// 
#ifndef Print_inc_StatusG4Printer_hh
#define Print_inc_StatusG4Printer_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class StatusG4Printer : public ProductPrinter {
  public:

    StatusG4Printer() {}
    StatusG4Printer(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<StatusG4>& handle, 
               std::ostream& os = std::cout);
    void Print(const art::ValidHandle<StatusG4>& handle, 
               std::ostream& os = std::cout);
    void Print(const mu2e::StatusG4& obj, 
	       int ind = -1, std::ostream& os = std::cout);
    void PrintHeader(const std::string& tag, 
                     std::ostream& os = std::cout);

  };

}
#endif

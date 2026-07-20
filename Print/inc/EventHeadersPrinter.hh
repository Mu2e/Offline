//
//  Utility class to print EventHeaders
//
#ifndef Print_inc_EventHeadersPrinter_hh
#define Print_inc_EventHeadersPrinter_hh

#include <iostream>

#include "Offline/Print/inc/ProductPrinter.hh"
#include "artdaq-core-mu2e/Data/EventHeader.hh"
#include "art/Framework/Principal/Handle.h"

namespace mu2e {

class EventHeadersPrinter : public ProductPrinter {
 public:
  EventHeadersPrinter() {}
  EventHeadersPrinter(const Config& conf) : ProductPrinter(conf) {}

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void Print(const art::Handle<EventHeaders>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<EventHeaders>& handle,
             std::ostream& os = std::cout);
  void Print(const mu2e::EventHeaders& obj,
             std::ostream& os = std::cout);
  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
};

}  // namespace mu2e
#endif

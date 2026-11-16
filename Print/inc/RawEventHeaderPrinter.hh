//
//  Utility class to print RawEventHeader
//
#ifndef Print_inc_RawEventHeaderPrinter_hh
#define Print_inc_RawEventHeaderPrinter_hh

#include <iostream>

#include "Offline/Print/inc/ProductPrinter.hh"
#include "artdaq-core/Data/RawEvent.hh"
#include "art/Framework/Principal/Handle.h"

namespace mu2e {

class RawEventHeaderPrinter : public ProductPrinter {
 public:
  RawEventHeaderPrinter() {}
  RawEventHeaderPrinter(const Config& conf) : ProductPrinter(conf) {}

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void Print(const art::Handle<artdaq::detail::RawEventHeader>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<artdaq::detail::RawEventHeader>& handle,
             std::ostream& os = std::cout);
  void Print(const artdaq::detail::RawEventHeader& obj,
             std::ostream& os = std::cout);
  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
};

}  // namespace mu2e
#endif

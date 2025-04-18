//
//  Utility class to print an EventHeader with conversions to physical units.
//
#ifndef Print_inc_EventHeaderPrinter_hh
#define Print_inc_EventHeaderPrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/Print/inc/ProductPrinter.hh"
#include "artdaq-core-mu2e/Data/EventHeader.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

class EventHeaderPrinter : public ProductPrinter {
 public:
  EventHeaderPrinter() {}
  EventHeaderPrinter(const Config& conf) : ProductPrinter(conf) {}

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void Print(const art::Handle<EventHeader>& handle, std::ostream& os = std::cout);
  void Print(const art::ValidHandle<EventHeader>& handle,
             std::ostream& os = std::cout);
  void Print(const mu2e::EventHeader& obj, int ind = -1,
             std::ostream& os = std::cout);
  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
};

}  // namespace mu2e
#endif

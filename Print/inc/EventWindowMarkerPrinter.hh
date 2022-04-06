//
//  Utility class to print EventWindowMarker
//
#ifndef Print_inc_EventWindowMarkerPrinter_hh
#define Print_inc_EventWindowMarkerPrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/Print/inc/ProductPrinter.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

class EventWindowMarkerPrinter : public ProductPrinter {
 public:
  EventWindowMarkerPrinter() {}
  EventWindowMarkerPrinter(const Config& conf) : ProductPrinter(conf) {}

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void Print(const art::Handle<EventWindowMarker>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<EventWindowMarker>& handle,
             std::ostream& os = std::cout);
  void Print(const mu2e::EventWindowMarker& obj, int ind = -1,
             std::ostream& os = std::cout);
  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
};

}  // namespace mu2e
#endif

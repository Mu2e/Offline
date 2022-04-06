//
//  Utility class to print TriggerInfo
//
#ifndef Print_inc_TriggerInfoPrinter_hh
#define Print_inc_TriggerInfoPrinter_hh

#include <cstring>
#include <iostream>

#include "Offline/Print/inc/ProductPrinter.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

class TriggerInfoPrinter : public ProductPrinter {
 public:
  TriggerInfoPrinter() {}
  TriggerInfoPrinter(const Config& conf) : ProductPrinter(conf) {}

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void Print(const art::Handle<mu2e::TriggerInfo>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<mu2e::TriggerInfo>& handle,
             std::ostream& os = std::cout);
  void Print(const mu2e::TriggerInfo& obj, int ind = -1,
             std::ostream& os = std::cout);
  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
};

}  // namespace mu2e
#endif

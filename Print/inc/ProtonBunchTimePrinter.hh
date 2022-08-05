//
//  Utility class to print ProtonBunchTime
//
#ifndef Print_inc_ProtonBunchTimePrinter_hh
#define Print_inc_ProtonBunchTimePrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/Print/inc/ProductPrinter.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

class ProtonBunchTimePrinter : public ProductPrinter {
 public:
  ProtonBunchTimePrinter() {}
  ProtonBunchTimePrinter(const Config& conf) : ProductPrinter(conf) {}

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void Print(const art::Handle<ProtonBunchTime>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<ProtonBunchTime>& handle,
             std::ostream& os = std::cout);
  void Print(const mu2e::ProtonBunchTime& obj, int ind = -1,
             std::ostream& os = std::cout);
  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
};

}  // namespace mu2e
#endif

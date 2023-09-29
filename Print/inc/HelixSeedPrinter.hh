//
//  Utility class to print HelixSeed
//
#ifndef Print_inc_HelixSeedPrinter_hh
#define Print_inc_HelixSeedPrinter_hh

#include <cstring>
#include <iostream>

#include "Offline/Print/inc/ProductPrinter.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

class HelixSeedPrinter : public ProductPrinter {
 public:
  HelixSeedPrinter() {}
  HelixSeedPrinter(const Config& conf) : ProductPrinter(conf) {}

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void Print(const art::Handle<HelixSeedCollection>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<HelixSeedCollection>& handle,
             std::ostream& os = std::cout);
  void Print(const HelixSeedCollection& coll, std::ostream& os = std::cout);
  void Print(const art::Ptr<HelixSeed>& ptr, int ind = -1,
             std::ostream& os = std::cout);
  void Print(const mu2e::HelixSeed& obj, int ind = -1,
             std::ostream& os = std::cout);

  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
  void PrintListHeader(std::ostream& os = std::cout);
};

}  // namespace mu2e
#endif

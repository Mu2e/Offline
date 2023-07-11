//
//  Utility class to print CosmicLivetime
//
#ifndef Print_inc_CosmicLivetimePrinter_hh
#define Print_inc_CosmicLivetimePrinter_hh

#include <cstring>
#include <iostream>

#include "Offline/MCDataProducts/inc/CosmicLivetime.hh"
#include "Offline/Print/inc/ProductPrinter.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

class CosmicLivetimePrinter : public ProductPrinter {
 public:
  CosmicLivetimePrinter() {}
  CosmicLivetimePrinter(const Config& conf) : ProductPrinter(conf) {}

  // all the ways to request a printout
  void PrintSubRun(art::SubRun const& subrun,
                   std::ostream& os = std::cout) override;
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void Print(const art::Handle<CosmicLivetime>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<CosmicLivetime>& handle,
             std::ostream& os = std::cout);
  void Print(const mu2e::CosmicLivetime& obj, int ind = -1,
             std::ostream& os = std::cout);
  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
};

}  // namespace mu2e
#endif

//
//  Utility class to print CosmicTrackSeed
//
#ifndef Print_inc_CosmicTrackSeedPrinter_hh
#define Print_inc_CosmicTrackSeedPrinter_hh

#include <cstring>
#include <iostream>

#include "Offline/Print/inc/ProductPrinter.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

class CosmicTrackSeedPrinter : public ProductPrinter {
 public:
  CosmicTrackSeedPrinter() {}
  CosmicTrackSeedPrinter(const Config& conf) : ProductPrinter(conf) {}

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void Print(const art::Handle<CosmicTrackSeedCollection>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<CosmicTrackSeedCollection>& handle,
             std::ostream& os = std::cout);
  void Print(const CosmicTrackSeedCollection& coll, std::ostream& os = std::cout);
  void Print(const art::Ptr<CosmicTrackSeed>& ptr, int ind = -1,
             std::ostream& os = std::cout);
  void Print(const mu2e::CosmicTrackSeed& obj, int ind = -1,
             std::ostream& os = std::cout);

  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
  void PrintListHeader(std::ostream& os = std::cout);
};

}  // namespace mu2e
#endif

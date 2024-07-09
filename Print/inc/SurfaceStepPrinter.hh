//
//  Utility class to print SurfaceStep
//
#ifndef Print_inc_SurfaceStepPrinter_hh
#define Print_inc_SurfaceStepPrinter_hh

#include <cstring>
#include <iostream>

#include "Offline/MCDataProducts/inc/SurfaceStep.hh"
#include "Offline/Print/inc/ProductPrinter.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

class SurfaceStepPrinter : public ProductPrinter {
 public:
  SurfaceStepPrinter(const Config& conf) : ProductPrinter(conf) {}

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void Print(const art::Handle<SurfaceStepCollection>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<SurfaceStepCollection>& handle,
             std::ostream& os = std::cout);
  void Print(const SurfaceStepCollection& coll, std::ostream& os = std::cout);
  void Print(const art::Ptr<SurfaceStep>& ptr, int ind = -1,
             std::ostream& os = std::cout);
  void Print(const mu2e::SurfaceStep& obj, int ind = -1,
             std::ostream& os = std::cout);

  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
  void PrintListHeader(std::ostream& os = std::cout);
};

}  // namespace mu2e
#endif

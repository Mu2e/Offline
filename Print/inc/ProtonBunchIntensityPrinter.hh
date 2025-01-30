//
//  Utility class to print ProtonBunchIntensity
//
#ifndef Print_inc_ProtonBunchIntensityPrinter_hh
#define Print_inc_ProtonBunchIntensityPrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "Offline/Print/inc/ProductPrinter.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

class ProtonBunchIntensityPrinter : public ProductPrinter {
 public:
  ProtonBunchIntensityPrinter() {}
  ProtonBunchIntensityPrinter(const Config& conf) : ProductPrinter(conf) {}

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void Print(const art::Handle<ProtonBunchIntensity>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<ProtonBunchIntensity>& handle,
             std::ostream& os = std::cout);
  void Print(const mu2e::ProtonBunchIntensity& obj, int ind = -1,
             std::ostream& os = std::cout);
  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
  void PrintEndJob(std::ostream& os = std::cout) override;
 private:
  double nevts_ = 0; // number of events processed
  double nPOT_ = 0; // total number of POT for these events
};

}  // namespace mu2e
#endif

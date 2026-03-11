//
//  Utility class to print PrescaleFilterFraction
//
#ifndef Print_inc_PrescaleFilterFractionPrinter_hh
#define Print_inc_PrescaleFilterFractionPrinter_hh

#include <cstring>
#include <iostream>

#include "Offline/DataProducts/inc/PrescaleFilterFraction.hh"
#include "Offline/Print/inc/ProductPrinter.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

class PrescaleFilterFractionPrinter : public ProductPrinter {
 public:
  PrescaleFilterFractionPrinter() {}
  PrescaleFilterFractionPrinter(const Config& conf) : ProductPrinter(conf)  {}

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void PrintSubRun(art::SubRun const& subrun,
                   std::ostream& os = std::cout) override;
  void Print(const art::Handle<PrescaleFilterFraction>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<PrescaleFilterFraction>& handle,
             std::ostream& os = std::cout);
  void Print(const mu2e::PrescaleFilterFraction& obj, int ind = -1,
             std::ostream& os = std::cout);
  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
 private:
};

}  // namespace mu2e
#endif

//
//  Utility class to print FilterFraction
//
#ifndef Print_inc_FilterFractionPrinter_hh
#define Print_inc_FilterFractionPrinter_hh

#include <cstring>
#include <iostream>

#include "Offline/DataProducts/inc/FilterFraction.hh"
#include "Offline/Print/inc/ProductPrinter.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

class FilterFractionPrinter : public ProductPrinter {
 public:
  FilterFractionPrinter() {}
  FilterFractionPrinter(const Config& conf) : ProductPrinter(conf)  {}

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void PrintSubRun(art::SubRun const& subrun,
                   std::ostream& os = std::cout) override;
  void Print(const art::Handle<FilterFraction>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<FilterFraction>& handle,
             std::ostream& os = std::cout);
  void Print(const mu2e::FilterFraction& obj, int ind = -1,
             std::ostream& os = std::cout);
  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
 private:
};

}  // namespace mu2e
#endif
